#Forest seed pilot
#Reported diagnostic: how much forest randomness alone moves delta C. Nothing here sets a main-run option.
#Design: HNSC and KIRP, the full clinical set + RNA, the first 5 repeats of the main run (same folds),
#5 forest seeds per repeat (FOREST_SEED_OFFSET 1..5; 0 is the main run), forest arm only.
#Measures, per tree count (250/500/1000/2000): forest_sd (spread across seeds, same folds), repeat_sd (spread across
#repeats), their ratio (below 0.1 = forest randomness negligible), block-weight stability and fit times.
#Needs the modelling objects for both cohorts. Usage: source("scripts/pilot_forest_seeds.R");
#run_forest_pilot(); analyse_forest_pilot()  (outputs in results/forest_pilot/)

source("paths.R")
source_library()

PILOT_COHORTS <- c("HNSC", "KIRP")
PILOT_SEEDS   <- 1:5          # FOREST_SEED_OFFSET values; 0 is reserved for the main run
PILOT_REPEATS <- 5L
PILOT_WORKERS <- N_WORKERS    # parallel over (fold x seed) tasks; 1 = sequential
PILOT_RATIO   <- 0.1          # criterion: forest_sd < 0.1 x repeat_sd
DIR_PILOT     <- file.path(DIR_RESULTS, "forest_pilot")
dir.create(DIR_PILOT, recursive = TRUE, showWarnings = FALSE)

# Arm switches for the pilot. They are applied ONLY while run_forest_pilot()
# runs and restored afterwards, so the normal pipeline in the same R session is
# unaffected. cv.R on disk is never touched.
.pilot_overrides <- list(RUN_FUSED_ARM = FALSE, RUN_OFFSET_ARM = FALSE,
                         RUN_FOREST_ARM = TRUE,
                         SAVE_PREDICTIONS = FALSE, SAVE_GENE_NAMES = FALSE)

path_pilot_ckpt <- function(code, r)
  file.path(DIR_PILOT, sprintf("%s_pilot_repeat_%02d.rds", code, r))

# settings that change pilot results; checkpoints from other settings are refused
.pilot_fingerprint <- function(code)
  paste(c(sprintf("code=%s seeds=%s repeats=%d seed=%d outer=%d",
                  code, paste(PILOT_SEEDS, collapse = "/"), PILOT_REPEATS, SEED, OUTER_FOLDS),
          sprintf("trees=%s nsets=%d trees_pre=%d splitrule=%s",
                  paste(FOREST_TREES, collapse = "/"), BF_NSETS, BF_NUM_TREES_PRE, RF_SPLITRULE),
          sprintf("n_genes=%s rna_completeness=%s grid=%d min_events=%d",
                  TOP_N_GENES, RNA_MIN_COMPLETENESS, length(GRID_QUANTILES), MIN_TRAIN_EVENTS),
          sprintf("vars=%s", paste(names(CANCER_VARS[[code]]), collapse = ","))),
        collapse = " | ")

#one task = one (fold, seed): runs the full-clinical cell
# ctx holds the fold's objects (built once per fold, shared across seeds).
.pilot_task <- function(task, ctx, subset, sid, code, r, seed_r) {
  assign("FOREST_SEED_OFFSET", task$s, envir = .GlobalEnv)   # read by forest_seed()
  cx  <- ctx[[task$f]]
  out <- run_one_cell(subset, sid, code, r, task$f, seed_r,
                      cx$trc, cx$tec, cx$rn, cx$tau, cx$grid, cx$rna_secs,
                      inner_foldid = NULL, inner_rna = NULL)   # only fused/offset use these
  row <- out$row
  row$pilot_seed <- task$s
  row
}

# Worker entry point. Defined at TOP LEVEL on purpose: a function created
# inside run_forest_pilot_one() would carry that function's whole frame --
# including the ~20,000-gene expression matrix -- to the workers with every
# task. Here it only refers to `ctx`, exported once per repeat.
.pilot_worker <- function(task, subset, sid, code, r, seed_r)
  .pilot_task(task, ctx, subset, sid, code, r, seed_r)

#fit the pilot for one cohort
run_forest_pilot_one <- function(code, cl_workers = NULL) {
  if (!file.exists(path_modelling(code)))
    stop("no modelling object for ", code, " -- run the 'prepare' stage first")
  fp <- .pilot_fingerprint(code)
  fppath <- file.path(DIR_PILOT, paste0(code, "_pilot_settings.txt"))
  if (file.exists(fppath) && !identical(paste(readLines(fppath), collapse = ""), fp))
    stop("Pilot checkpoints for ", code, " were written under different settings. ",
         "Delete ", DIR_PILOT, "/", code, "_pilot_* and re-run.")
  writeLines(fp, fppath)
  
  dat  <- readRDS(path_modelling(code))
  ok   <- !is.na(dat$clinical$OS) & !is.na(dat$clinical$OS.time) & dat$clinical$OS.time > 0
  clin <- dat$clinical[ok, ]; expr <- dat$expr[ok, , drop = FALSE]
  status <- clin$OS
  
  vars   <- names(CANCER_VARS[[code]])
  real   <- all_subsets(vars)
  sid    <- which(vapply(real, function(s) setequal(s, vars), logical(1)))
  if (length(sid) != 1L) stop("Could not locate the full clinical subset for ", code)
  subset <- real[[sid]]
  cat(sprintf("=== %s: %d patients, %d events | cell = full set (subset_id %d, k = %d) ===\n",
              code, nrow(clin), sum(status), sid, length(subset)))
  
  for (r in seq_len(PILOT_REPEATS)) {
    ckpt <- path_pilot_ckpt(code, r)
    if (file.exists(ckpt)) { cat(sprintf("  repeat %d: checkpoint found\n", r)); next }
    seed_r <- SEED + r                                 # identical to the main run
    folds  <- make_outer_folds(status, seed = seed_r)
    t0 <- proc.time()[["elapsed"]]
    
    # fold objects, built once and shared by all seeds (same rules as cv.R)
    ctx <- vector("list", OUTER_FOLDS)
    for (f in seq_len(OUTER_FOLDS)) {
      tei <- which(folds == f); tri <- which(folds != f)
      trc <- clin[tri, ]; tec <- clin[tei, ]
      if (sum(trc$OS) < MIN_TRAIN_EVENTS) { cat(sprintf("  fold %d: SKIPPED (train events)\n", f)); next }
      tau  <- make_tau(trc$OS.time, trc$OS)
      grid <- sort(unique(quantile(trc$OS.time[trc$OS == 1], GRID_QUANTILES)))
      grid <- grid[grid <= tau]
      if (length(grid) < 2) { cat(sprintf("  fold %d: SKIPPED (grid)\n", f)); next }
      trna <- timed(build_rna(expr[tri, , drop = FALSE], expr[tei, , drop = FALSE], TOP_N_GENES))
      ctx[[f]] <- list(trc = trc, tec = tec, rn = trna$value, tau = tau,
                       grid = grid, rna_secs = trna$secs)
    }
    live  <- which(!vapply(ctx, is.null, logical(1)))
    tasks <- expand.grid(f = live, s = PILOT_SEEDS)
    tasks <- split(tasks, seq_len(nrow(tasks)))
    
    rows <- if (is.null(cl_workers)) {
      lapply(tasks, .pilot_task, ctx = ctx, subset = subset, sid = sid,
             code = code, r = r, seed_r = seed_r)
    } else {
      parallel::clusterExport(cl_workers, "ctx", envir = environment())
      parallel::parLapplyLB(cl_workers, tasks, .pilot_worker,
                            subset = subset, sid = sid, code = code, r = r, seed_r = seed_r)
    }
    res <- do.call(rbind, rows)
    saveRDS(res, ckpt)
    cat(sprintf("  repeat %d/%d: %d cells in %.1f min | blockForest errors %d\n",
                r, PILOT_REPEATS, nrow(res), (proc.time()[["elapsed"]] - t0) / 60,
                sum(!is.na(res$bf_error))))
  }
  invisible(TRUE)
}

run_forest_pilot <- function(codes = PILOT_COHORTS) {
  # apply the pilot's switches for the duration of this call, then restore
  keep <- c(names(.pilot_overrides), "FOREST_SEED_OFFSET")
  old  <- mget(keep, envir = .GlobalEnv)
  on.exit(for (nm in keep) assign(nm, old[[nm]], envir = .GlobalEnv), add = TRUE)
  for (nm in names(.pilot_overrides)) assign(nm, .pilot_overrides[[nm]], envir = .GlobalEnv)
  
  cl_workers <- NULL
  if (PILOT_WORKERS > 1) {
    cl_workers <- parallel::makeCluster(PILOT_WORKERS)
    on.exit(tryCatch(parallel::stopCluster(cl_workers), error = function(e) NULL), add = TRUE)
    parallel::clusterEvalQ(cl_workers, { source("paths.R"); source_library(); NULL })
    # the pilot's arm switches must reach the workers, which sourced cv.R afresh
    parallel::clusterExport(cl_workers, c(".pilot_task", ".pilot_worker",
                                          names(.pilot_overrides)),
                            envir = .GlobalEnv)
  }
  for (code in codes) run_forest_pilot_one(code, cl_workers)
  cat("Pilot fitting done. Next: analyse_forest_pilot()\n")
  invisible(TRUE)
}

#analysis
.load_pilot <- function(code) {
  files <- vapply(seq_len(PILOT_REPEATS), function(r) path_pilot_ckpt(code, r), character(1))
  files <- files[file.exists(files)]
  if (!length(files)) return(NULL)
  do.call(rbind, lapply(files, readRDS))
}

analyse_forest_pilot <- function(codes = PILOT_COHORTS) {
  out <- list()
  for (code in codes) {
    res <- .load_pilot(code)
    if (is.null(res)) { cat(code, ": no pilot checkpoints\n"); next }
    n_rep <- length(unique(res$repeat_id))
    
    # block-weight stability: SD across seeds within each (repeat, fold)
    w <- aggregate(bf_w_clin ~ repeat_id + fold, data = res, FUN = sd)
    w_med <- median(res$bf_w_clin, na.rm = TRUE)
    t_bf  <- median(res$secs_bf, na.rm = TRUE)
    t_rf  <- median(res$secs_rfclin, na.rm = TRUE)
    
    for (nt in FOREST_TREES) {
      for (metric in c("delta_c_bf", "delta_ibs_bf")) {
        col <- paste0(metric, "_", nt)
        d <- res[!is.na(res[[col]]), c("repeat_id", "pilot_seed", "fold", col)]
        names(d)[4] <- "v"
        # 1) fold-collapse within (repeat, seed) -- same order as 02_summarise.R
        fm <- aggregate(v ~ repeat_id + pilot_seed, data = d, FUN = mean)
        # 2) forest noise: seed spread within repeat, pooled over repeats
        vr <- aggregate(v ~ repeat_id, data = fm, FUN = var)
        forest_sd <- sqrt(mean(vr$v, na.rm = TRUE))
        # 3) fold-split noise: spread across repeats of the seed-averaged value
        sm <- aggregate(v ~ repeat_id, data = fm, FUN = mean)
        repeat_sd <- sd(sm$v)
        out[[length(out) + 1]] <- data.frame(
          cancer = code, trees = nt, metric = metric,
          n_repeats = n_rep, n_seeds = length(unique(fm$pilot_seed)),
          mean = mean(sm$v), forest_sd = forest_sd, repeat_sd = repeat_sd,
          ratio = forest_sd / repeat_sd,
          passes = (forest_sd / repeat_sd) < PILOT_RATIO,
          bf_w_clin_median = w_med,
          bf_w_clin_seed_sd_median = median(w$bf_w_clin, na.rm = TRUE),
          bf_w_clin_seed_sd_max    = max(w$bf_w_clin, na.rm = TRUE),
          secs_bf_median = t_bf, secs_rfclin_median = t_rf,
          stringsAsFactors = FALSE)
      }
    }
  }
  if (!length(out)) return(invisible(NULL))
  tab <- do.call(rbind, out)
  path <- file.path(DIR_PILOT, "forest_pilot_summary.csv")
  write.csv(tab, path, row.names = FALSE)
  
  cat("\n=== forest seed pilot: delta C (criterion: forest_sd < ",
      PILOT_RATIO, " x repeat_sd) ===\n", sep = "")
  print(tab[tab$metric == "delta_c_bf",
            c("cancer", "trees", "mean", "forest_sd", "repeat_sd", "ratio", "passes")],
        row.names = FALSE, digits = 3)
  cat("\n=== blockForest clinical-block weight (seed stability) ===\n")
  print(unique(tab[, c("cancer", "bf_w_clin_median",
                       "bf_w_clin_seed_sd_median", "bf_w_clin_seed_sd_max")]),
        row.names = FALSE, digits = 3)
  
  # rough main-run cost: every subset x fold x repeat, spread over N_WORKERS
  cat("\n=== rough main-run cost of the forest arm (median cell time; full-set cells are the slowest, so an upper-side guess) ===\n")
  for (code in unique(tab$cancer)) {
    x <- tab[tab$cancer == code, ][1, ]
    n_cells <- 2^length(CANCER_VARS[[code]]) * OUTER_FOLDS * n_repeats_for(code)
    hrs <- n_cells * (x$secs_bf_median + x$secs_rfclin_median) / 3600 / max(1, N_WORKERS)
    cat(sprintf("  %s: %d cells x %.1fs  ->  ~%.1f h on %d workers\n",
                code, n_cells, x$secs_bf_median + x$secs_rfclin_median, hrs, N_WORKERS))
  }
  cat("\nSaved ", path, "\n", sep = "")
  invisible(tab)
}

cat("pilot_forest_seeds.R loaded.  run_forest_pilot()  then  analyse_forest_pilot()\n")