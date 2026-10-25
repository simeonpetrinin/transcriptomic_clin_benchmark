#Secondary runs: IPF-ridge and RNA score on the full lattice; blockForest and CoxBoost on the robustness sample
#Two jobs; the main lattice (IPF-LASSO, offset, clinical Cox) is NOT rerun.
#  A. ridge and RNA-score arms over the full lattice, key <CODE>_rs (IPF-LASSO and offset switched off).
#  B. blockForest and CoxBoost robustness: the 3 fixed specs plus 80 nested sampled subsets per cancer, key <CODE>_robust.
#CoxBoost runs unchanged: cells it cannot fit stay recorded as failed (about 10% in the old robustness run).
#Nothing in R/ is changed: the arm switches are set here and exported to the workers by run_cancer().
#Needs scripts/run_robustness.R with the nested sampler and 80 sampled subsets; stops otherwise.
#Run from the project root in a fresh R session: source("scripts/run_weekend.R"), then the printed steps.

Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1")  #one BLAS thread per worker; workers inherit this
source("master.R")
if (!exists("summarise_cancer", mode = "function")) source(file.path(DIR_SCRIPTS, "02_summarise.R"))

#guard against an old build of run_robustness.R
if (!exists("check_robust_nesting", mode = "function") || !identical(ROBUST_N_SAMPLED, 80L))
  stop("scripts/run_robustness.R is not the version with the nested sampler and 80 sampled subsets")

BACKUP_DIR <- file.path(DIR_RESULTS, "before_weekend")  #old robustness outputs are moved here, never deleted
PREP_MARK  <- file.path(BACKUP_DIR, "robust_prepared.txt")  #written once prepare_robust() has run

WEEKEND_ARMS <- c("forest", "coxboost")  #job B learners

#estimated CPU-hours per cancer for job B (from the earlier 13-spec robustness run, scaled to 83 specs); used only for the run order
FOREST_COST_H   <- c(BRCA = 118.8, LUAD = 99.8, LUSC = 74.8, KIRC = 73.6, SKCM = 73.1, BLCA = 70.4,
                     UCEC = 69.3, LIHC = 48.0, HNSC = 46.4, CESC = 33.5, STAD = 30.4, SARC = 30.4,
                     KIRP = 26.2, COAD = 24.8, LGG = 23.6, OV = 13.3, ESCA = 10.2, READ = 10.2,
                     PAAD = 9.4, GBM = 1.0, LAML = 0.8)  #blockForest + clinical RSF + fused + offset
COXBOOST_COST_H <- c(BRCA = 6.4, LUAD = 12.9, LUSC = 13.2, KIRC = 10.1, SKCM = 12.2, BLCA = 11.0,
                     UCEC = 5.0, LIHC = 7.9, HNSC = 8.2, CESC = 3.1, STAD = 4.6, SARC = 4.0,
                     KIRP = 1.4, COAD = 2.7, LGG = 7.1, OV = 6.8, ESCA = 1.1, READ = 0.9,
                     PAAD = 1.6, GBM = 0.3, LAML = 0.3)  #failed cells stop early, so these are low where many fail
ROBUST_COST_H <- sort(FOREST_COST_H + COXBOOST_COST_H[names(FOREST_COST_H)], decreasing = TRUE)
ROBUST_ORDER  <- names(ROBUST_COST_H)  #most expensive first
if (!setequal(ROBUST_ORDER, CANCERS))
  warning("ROBUST_ORDER and CANCERS differ: ", paste(union(setdiff(ROBUST_ORDER, CANCERS),
                                                           setdiff(CANCERS, ROBUST_ORDER)), collapse = ", "))
ARM_SWITCHES <- c("RUN_FUSED_ARM", "RUN_OFFSET_ARM", "RUN_FOREST_ARM", "RUN_COXBOOST_ARM",
                  "RUN_RIDGE_ARM", "RUN_SCORE_ARM")

rs_key <- function(code) paste0(code, "_rs")

#set the worker count (and optionally arm switches) for one call, restore them on exit
.with_settings <- function(workers, switches, expr) {
  old <- mget(c("N_WORKERS", ARM_SWITCHES), envir = .GlobalEnv)
  on.exit(for (n in names(old)) assign(n, old[[n]], envir = .GlobalEnv), add = TRUE)
  assign("N_WORKERS", as.integer(workers), envir = .GlobalEnv)
  for (n in names(switches)) assign(n, switches[[n]], envir = .GlobalEnv)
  force(expr)
}

#Job A: ridge + RNA score over the full lattice, same seeds and folds as the main run
weekend_ridge_score <- function(codes = intersect(ROBUST_ORDER, setdiff(CANCERS, "LAML")), workers = N_WORKERS) {
  sw <- list(RUN_FUSED_ARM = FALSE, RUN_OFFSET_ARM = FALSE, RUN_FOREST_ARM = FALSE,
             RUN_COXBOOST_ARM = FALSE, RUN_RIDGE_ARM = TRUE, RUN_SCORE_ARM = TRUE)
  failed <- character(0)
  .with_settings(workers, sw, for (code in codes) {
    key <- rs_key(code)
    if (file.exists(path_results(key))) {
      cat(sprintf("[ridge+score] %-5s skipped (%s exists)\n", code, basename(path_results(key)))); next
    }
    ok <- tryCatch({
      timed_stage(code, "ridge_score", run_cancer(readRDS(path_modelling(code)), code, key = key))
      TRUE
    }, error = function(e) { cat(sprintf("[FAILED]      %-5s ridge+score: %s\n", code, conditionMessage(e))); FALSE })
    if (ok) tryCatch(summarise_cancer(key), error = function(e)
      cat(sprintf("[summary]     %-5s not written: %s\n", code, conditionMessage(e))))
    if (!ok) failed <- c(failed, code)
    invisible(gc())
  })
  if (length(failed)) cat("[ridge+score] failed: ", paste(failed, collapse = ", "), "\n", sep = "")
  invisible(list(failed = failed))
}

#the old <CODE>_robust files (results, summaries, timing, checkpoints, figure); the _g100/_gall gene runs are left alone
.robust_files <- function(code) {
  k <- paste0("^", code, "_robust_")
  c(list.files(DIR_RESULTS, full.names = TRUE, pattern = paste0(k,
                                                                "(results|genes|preds|manifest|specs|subset_summary|per_repeat|timing_summary)(_norescue)?\\.")),
    list.files(file.path(DIR_RESULTS, "timing"), full.names = TRUE,
               pattern = paste0(k, "(timing|timing_repeats|fold_timing)\\.csv$")),
    list.files(DIR_CHECKPOINTS, full.names = TRUE, pattern = paste0(k, "(repeat_[0-9]+\\.rds|settings\\.txt)$")),
    list.files(DIR_FIGURES, full.names = TRUE, pattern = paste0(k, "summary_check")))
}

#Job B, once before it starts: move the earlier 13-spec robustness outputs to results/before_weekend, write the 83-spec file
prepare_robust <- function(codes = CANCERS) {
  if (file.exists(PREP_MARK))
    stop("prepare_robust() already ran (", PREP_MARK, "); running it again would move the new results away")
  n_moved <- 0L
  for (code in codes) {
    f <- .robust_files(code)
    dest <- file.path(BACKUP_DIR, basename(dirname(f)), basename(f))  #results/, timing/, checkpoints/, figures/
    if (any(file.exists(dest))) stop("backup already holds ", dest[file.exists(dest)][1])
    for (d in unique(dirname(dest))) dir.create(d, recursive = TRUE, showWarnings = FALSE)
    ok <- file.copy(f, dest, copy.date = TRUE)
    if (!all(ok)) stop("could not copy ", f[!ok][1], "; nothing deleted for ", code)
    unlink(f); n_moved <- n_moved + length(f)
    cat(sprintf("[prepare]   %-5s %d old robustness file(s) moved\n", code, length(f)))
  }
  old_specs <- file.path(DIR_RESULTS, "robustness_specs.csv")
  if (file.exists(old_specs)) file.rename(old_specs, file.path(BACKUP_DIR, "robustness_specs.csv"))
  write_robust_specs(codes)
  writeLines(sprintf("prepared %s, %d files moved", format(Sys.time()), n_moved), PREP_MARK)
  invisible(n_moved)
}

#Job B: blockForest + CoxBoost robustness (fused and offset run beside them on the same cells, as before)
weekend_robust <- function(codes = intersect(ROBUST_ORDER, setdiff(CANCERS, "LAML")), workers = N_WORKERS,
                           arms = WEEKEND_ARMS) {
  if (!file.exists(PREP_MARK)) stop("run prepare_robust() first")
  .with_settings(workers, list(), full_run_robustness(codes, arms = arms, force = FALSE))
}

#everything left after the LAML checks: job A, then job B; safe to stop and call again (checkpoints)
weekend_rest <- function(workers = N_WORKERS) {
  weekend_ridge_score(workers = workers)
  weekend_robust(workers = workers)
}

#matched cells (same variables, repeat, fold) of two results files; descriptive only, you judge
compare_runs <- function(new_file, old_file, cols) {
  new <- read.csv(new_file, stringsAsFactors = FALSE); old <- read.csv(old_file, stringsAsFactors = FALSE)
  k <- c("variables", "repeat_id", "fold")
  cols <- intersect(cols, intersect(names(new), names(old)))
  m <- merge(old[c(k, cols)], new[c(k, cols)], by = k, suffixes = c(".old", ".new"))
  cat(sprintf("%s vs %s: %d matched cells\n", basename(new_file), basename(old_file), nrow(m)))
  for (v in cols) {
    a <- m[[paste0(v, ".old")]]; b <- m[[paste0(v, ".new")]]; both <- !is.na(a) & !is.na(b)
    cat(sprintf("  %-16s both %5d | equal (1e-9) %5d | max |diff| %.3g | NA old only %d, new only %d\n",
                v, sum(both), sum(abs(a[both] - b[both]) < 1e-9),
                if (any(both)) max(abs(a[both] - b[both])) else NA_real_,
                sum(is.na(a) & !is.na(b)), sum(is.na(b) & !is.na(a))))
  }
  invisible(m)
}

#check of job A: clinical side against the main lattice, and how many ridge/score cells fitted
check_ridge_score <- function(code = "LAML") {
  compare_runs(path_results(rs_key(code)), path_results(code), c("c_clin", "ibs_clin"))
  new <- read.csv(path_results(rs_key(code)), stringsAsFactors = FALSE)
  for (a in c("ridge", "score"))
    cat(sprintf("  %-6s fitted %d of %d cells\n", a, sum(!is.na(new[[paste0("c_", a)]])), nrow(new)))
}

#check of job B: clinical side against the main lattice; forest, CoxBoost, fused, offset against the old robustness run
check_robust <- function(code = "LAML") {
  new_f <- path_results(robust_key(code))
  compare_runs(new_f, path_results(code), c("c_clin", "ibs_clin"))
  old_f <- file.path(BACKUP_DIR, "results", basename(new_f))
  if (file.exists(old_f))
    compare_runs(new_f, old_f, c("c_fused", "c_offset", "c_rfclin_2000", "c_bf_2000", "ibs_bf_2000", "c_cb"))
  new <- read.csv(new_f, stringsAsFactors = FALSE)
  cat(sprintf("  CoxBoost fitted %d of %d cells\n", sum(!is.na(new$c_cb)), nrow(new)))
}

#this machine's forest speed against the earlier run (cells both runs fitted) and the projected wall time for job B
forest_speed <- function(code = "LAML", workers = N_WORKERS) {
  m <- compare_runs(path_results(robust_key(code)),
                    file.path(BACKUP_DIR, "results", paste0(robust_key(code), "_results.csv")),
                    c("secs_bf", "secs_rfclin"))
  t_old <- rowSums(m[c("secs_bf.old", "secs_rfclin.old")], na.rm = TRUE)
  t_new <- rowSums(m[c("secs_bf.new", "secs_rfclin.new")], na.rm = TRUE)
  ratio <- median(t_new / t_old, na.rm = TRUE)
  todo  <- intersect(ROBUST_ORDER, CANCERS)
  left  <- todo[!file.exists(path_results(robust_key(todo)))]
  cat(sprintf("forest cell time here / earlier run = %.2f (median of %d cells)\n", ratio, nrow(m)))
  cat(sprintf("job B left: %d cancers, about %.0f h at %d workers (rough: earlier estimate x ratio / workers)\n",
              length(left), sum(ROBUST_COST_H[left]) * ratio / workers, workers))
  invisible(ratio)
}

cat("\nrun_weekend.R loaded (workers =", N_WORKERS, "from settings.txt)\n")
cat("  1. prepare_robust()                # once: moves the old robustness files, writes the 80-subset specs\n")
cat("  2. weekend_ridge_score(\"LAML\"); check_ridge_score(\"LAML\")   # minutes\n")
cat("  3. weekend_robust(\"LAML\");      check_robust(\"LAML\"); forest_speed(\"LAML\")\n")
cat("  4. weekend_rest()                  # job A for the other 20, then job B; call again to resume\n")