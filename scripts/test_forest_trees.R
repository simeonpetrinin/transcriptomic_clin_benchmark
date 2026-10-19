#Forest tree-count test
#How many trees blockForest needs, in the weight search and in the final forest, on
#repeat 1 fold 1 of each cancer (full main set). Saves only differences, spreads and
#seconds to results/timing/; no C, delta C, OOB error or chosen weight is shown.

TT_CODES        <- c("HNSC", "KIRP")
TT_SEEDS        <- 1:10
TT_SEARCH_TREES <- 1500L  #trees per candidate; an earlier run compared 250, 500 and 1500
TT_FINAL_TREES  <- c(100L, 200L, 250L, 300L, 500L, 1000L, 1500L, 2000L, 3000L)  #compared with the largest
TT_GRID         <- seq(0.05, 1, by = 0.05)  #RNA weight candidates, clinical fixed at 1

#repeat 1 fold 1, built exactly as run_one_repeat() does
.tt_fold <- function(code) {
  dat  <- readRDS(path_modelling(code))
  ok   <- valid_survival(dat$clinical)
  clin <- dat$clinical[ok, ]; expr <- dat$expr[ok, , drop = FALSE]
  vars <- names(CANCER_VARS[[code]])
  RARE_MISSING_VARS <<- rare_missing_vars(clin, vars, MIN_LEVEL_N)
  folds <- make_outer_folds(clin$OS, seed = SEED + 1L)
  trc <- clin[folds != 1, ]; tec <- clin[folds == 1, ]
  tau  <- make_tau(trc$OS.time, trc$OS)
  grid <- sort(unique(quantile(trc$OS.time[trc$OS == 1], GRID_QUANTILES)))
  grid <- grid[grid <= tau]
  rn <- build_rna(expr[folds != 1, , drop = FALSE], expr[folds == 1, , drop = FALSE], TOP_N_GENES)
  cl <- build_clinical(vars, trc, tec)
  fx <- fuse(cl, rn)
  Xtr <- as.matrix(fx$train); Xte <- as.matrix(fx$test)
  safe <- paste0("V", seq_len(ncol(Xtr))); colnames(Xtr) <- safe; colnames(Xte) <- safe
  blocks <- list(clinical = which(fx$block == "clinical"), rna = which(fx$block == "rna"))
  mtry <- as.integer(c(length(blocks$clinical), sqrt(length(blocks$rna))))  #as .blockfor_fixed_clinical
  list(code = code, trc = trc, tec = tec, tau = tau, grid = grid, Xtr = Xtr, Xte = Xte,
       Xcl_tr = as.matrix(cl$train), Xcl_te = as.matrix(cl$test), blocks = blocks, mtry = mtry)
}

#one search: OOB error of every grid weight at nt trees, one seed (runs on a worker)
.tt_search <- function(job, fd) {
  set.seed(SEED + 7000000L + job$seed)
  df <- data.frame(time = fd$trc$OS.time, status = fd$trc$OS, fd$Xtr, check.names = FALSE)
  t0 <- proc.time()[["elapsed"]]
  err <- vapply(TT_GRID, function(w)
    blockForest::blockForest(dependent.variable.name = "time", status.variable.name = "status",
                             data = df, num.trees = job$nt, blocks = fd$blocks,
                             block.weights = c(1, w), mtry = fd$mtry, keep.inbag = TRUE,
                             block.method = "BlockForest", splitrule = RF_SPLITRULE,
                             write.forest = FALSE, num.threads = RF_THREADS)$prediction.error,
    numeric(1))
  list(nt = job$nt, seed = job$seed, err = err,
       secs_per_forest = (proc.time()[["elapsed"]] - t0) / length(TT_GRID))
}

#one seed of the final forests: blockForest at weight w and the clinical RSF (worker)
.tt_final <- function(seed, fd, w) {
  set.seed(SEED + 8000000L + seed)
  df <- data.frame(time = fd$trc$OS.time, status = fd$trc$OS, fd$Xtr, check.names = FALSE)
  t0 <- proc.time()[["elapsed"]]
  bf <- blockForest::blockForest(dependent.variable.name = "time", status.variable.name = "status",
                                 data = df, num.trees = max(TT_FINAL_TREES), blocks = fd$blocks,
                                 block.weights = c(1, w), mtry = fd$mtry, keep.inbag = TRUE,
                                 block.method = "BlockForest", splitrule = RF_SPLITRULE,
                                 num.threads = RF_THREADS)
  secs_bf <- proc.time()[["elapsed"]] - t0
  pb <- .forest_preds(bf, fd$Xte, fd$grid, TT_FINAL_TREES)
  rf <- fit_rsf(fd$Xcl_tr, fd$trc$OS.time, fd$trc$OS, fd$Xcl_te, fd$grid,
                seed = SEED + 9000000L + seed, mtry = ncol(fd$Xcl_tr))
  if (!is.na(rf$error)) stop("clinical RSF failed: ", rf$error)
  sc <- function(p) c(uno_c(fd$tec$OS.time, fd$tec$OS, p$risk, fd$tau),
                      as.numeric(integrated_brier(fd$tec$OS.time, fd$tec$OS, p$surv, fd$grid,
                                                  fd$tau, cens_time = fd$trc$OS.time,
                                                  cens_status = fd$trc$OS)))
  do.call(rbind, lapply(as.character(TT_FINAL_TREES), function(nt) {
    b <- sc(pb[[nt]]); r <- sc(rf$by_trees[[nt]])
    data.frame(seed = seed, trees = as.integer(nt), c_bf = b[1], ibs_bf = b[2],
               c_rf = r[1], ibs_rf = r[2], secs_bf_max = secs_bf)
  }))
}

#run the test on every cancer in codes, in parallel over seeds
test_forest_trees <- function(codes = TT_CODES, workers = N_WORKERS) {
  dir <- file.path(DIR_RESULTS, "timing"); dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  cw <- parallel::makeCluster(workers); on.exit(parallel::stopCluster(cw), add = TRUE)
  parallel::clusterEvalQ(cw, { source("paths.R"); source_library(); NULL })
  parallel::clusterExport(cw, c("TT_GRID", "TT_FINAL_TREES", ".tt_search", ".tt_final"))
  parallel::clusterEvalQ(cw, { FOREST_TREES <- TT_FINAL_TREES; NULL })  #RSF predictions at every test count
  s_rows <- list(); f_rows <- list(); t_rows <- list()
  for (code in codes) {
    t0 <- proc.time()[["elapsed"]]
    fd <- .tt_fold(code)
    cat(sprintf("\n[trees] %s: %d train / %d test patients, %d clinical + %d RNA columns\n",
                code, nrow(fd$trc), nrow(fd$tec), length(fd$blocks$clinical), length(fd$blocks$rna)))
    #1. weight search at each tree count
    jobs <- lapply(seq_len(length(TT_SEARCH_TREES) * length(TT_SEEDS)), function(i)
      list(nt = TT_SEARCH_TREES[(i - 1) %/% length(TT_SEEDS) + 1],
           seed = TT_SEEDS[(i - 1) %% length(TT_SEEDS) + 1]))
    parallel::clusterExport(cw, "fd", envir = environment())
    sr <- parallel::parLapplyLB(cw, jobs, function(j) .tt_search(j, fd))
    E  <- lapply(TT_SEARCH_TREES, function(nt) do.call(rbind, lapply(Filter(function(x) x$nt == nt, sr), `[[`, "err")))
    names(E) <- TT_SEARCH_TREES
    pick <- function(e) TT_GRID[which(e == min(e))[1]]
    ref_w <- median(apply(E[[as.character(max(TT_SEARCH_TREES))]], 1, pick))  #used below, not shown
    for (nt in TT_SEARCH_TREES) {
      e  <- E[[as.character(nt)]]; ch <- apply(e, 1, pick)
      rho <- suppressWarnings(cor(t(e), method = "spearman"))
      s_rows[[length(s_rows) + 1]] <- data.frame(
        cancer = code, search_trees = nt,
        weight_spread = max(ch) - min(ch),  #max minus min chosen RNA weight over seeds
        n_distinct_weights = length(unique(ch)),
        rank_agreement = mean(rho[upper.tri(rho)], na.rm = TRUE),  #mean seed-to-seed Spearman
        signal_to_noise = diff(range(colMeans(e))) / mean(apply(e, 2, sd)),
        secs_per_forest = mean(vapply(Filter(function(x) x$nt == nt, sr), `[[`, numeric(1), "secs_per_forest")))
    }
    #2. final forests at the reference weight
    fr <- do.call(rbind, parallel::parLapplyLB(cw, TT_SEEDS, function(s) .tt_final(s, fd, ref_w)))
    fr$d_bf <- fr$c_bf - fr$c_rf
    top <- fr[fr$trees == max(TT_FINAL_TREES), ]
    for (nt in TT_FINAL_TREES) {
      x <- fr[fr$trees == nt, ]; y <- top[match(x$seed, top$seed), ]
      f_rows[[length(f_rows) + 1]] <- data.frame(
        cancer = code, trees = nt,
        c_bf_diff_mean = mean(abs(x$c_bf - y$c_bf)), c_bf_diff_max = max(abs(x$c_bf - y$c_bf)),
        c_rf_diff_mean = mean(abs(x$c_rf - y$c_rf)), c_rf_diff_max = max(abs(x$c_rf - y$c_rf)),
        delta_c_diff_mean = mean(abs(x$d_bf - y$d_bf)),
        ibs_bf_diff_mean = mean(abs(x$ibs_bf - y$ibs_bf)),
        c_bf_seed_sd = sd(x$c_bf), c_rf_seed_sd = sd(x$c_rf), delta_c_seed_sd = sd(x$d_bf),
        ibs_bf_seed_sd = sd(x$ibs_bf))
    }
    t_rows[[length(t_rows) + 1]] <- data.frame(
      cancer = code, n_train = nrow(fd$trc),
      secs_final_bf_max = median(top$secs_bf_max),
      test_wall_min = round((proc.time()[["elapsed"]] - t0) / 60, 1))
  }
  S <- do.call(rbind, s_rows); Fn <- do.call(rbind, f_rows); Tm <- do.call(rbind, t_rows)
  #projected seconds for one blockForest fit under each search setting
  P <- do.call(rbind, lapply(seq_len(nrow(S)), function(i) {
    s <- S[i, ]  #final forest scaled to the 2000 trees the pipeline uses
    fin <- Tm$secs_final_bf_max[Tm$cancer == s$cancer] * 2000 / max(TT_FINAL_TREES)
    data.frame(cancer = s$cancer, search_trees = s$search_trees,
               min_300_random = round((300 * s$secs_per_forest + fin) / 60, 1),
               min_50_random  = round((50 * s$secs_per_forest + fin) / 60, 1),
               min_20_grid    = round((20 * s$secs_per_forest + fin) / 60, 1))
  }))
  write.csv(S, file.path(dir, "forest_trees_search.csv"), row.names = FALSE)
  write.csv(Fn, file.path(dir, "forest_trees_final.csv"), row.names = FALSE)
  write.csv(Tm, file.path(dir, "forest_trees_timing.csv"), row.names = FALSE)
  write.csv(P, file.path(dir, "forest_trees_projection.csv"), row.names = FALSE)
  cat("\n== weight search: stability by trees per candidate ==\n"); print(S, row.names = FALSE, digits = 3)
  cat("\n== final forest: difference from the largest count and seed spread ==\n"); print(Fn, row.names = FALSE, digits = 3)
  cat("\n== minutes per blockForest fit (single thread) ==\n"); print(P, row.names = FALSE)
  print(Tm, row.names = FALSE)
  invisible(list(search = S, final = Fn, timing = Tm, projection = P))
}