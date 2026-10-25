#Example run
#Reruns one cancer (LAML or GBM) from data/processed/<CODE>_modelling.rds and compares it with the paper:
#  1. the main lattice (IPF-LASSO and offset arms) against example/<CODE>_expected_subset_summary.csv;
#  2. the IPF-ridge and RNA-score arms over the same lattice, as in scripts/run_weekend.R (job A),
#     against example/<CODE>_expected_rs_subset_summary.csv.
#Outputs use the keys <CODE>_example and <CODE>_example_rs. Run from the project root: Rscript scripts/run_example.R LAML
#(or GBM); source("scripts/run_example.R") runs LAML.

source("master.R")
source(file.path(DIR_SCRIPTS, "02_summarise.R"))

EXAMPLE_CODE <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(EXAMPLE_CODE)) EXAMPLE_CODE <- "LAML"
if (!EXAMPLE_CODE %in% c("LAML", "GBM")) stop("the example runs LAML or GBM")
EXAMPLE_TOL <- 1e-3  #largest allowed absolute difference; results are reported to three decimals

#remove outputs and checkpoints of an earlier example run
clear_example <- function(key) {
  unlink(c(path_results(key), path_summary(key), path_timing(key), path_manifest(key),
           sub("_results\\.csv$", "_genes.csv.gz", path_results(key)),
           sub("_results\\.csv$", "_preds.rds", path_results(key)),
           list.files(DIR_CHECKPOINTS, pattern = paste0("^", key, "_"), full.names = TRUE)))
}

#compare a rerun summary with the paper, per model and clinical subset; TRUE if every quantity is within EXAMPLE_TOL
compare_example <- function(key, expected) {
  new  <- read.csv(path_summary(key), stringsAsFactors = FALSE)
  old  <- read.csv(expected, stringsAsFactors = FALSE)
  cols <- c("c_ref", "c_model", "delta_c", "ibs_ref", "ibs_model", "delta_ibs")
  m <- merge(old[, c("arm", "subset_id", cols)], new[, c("arm", "subset_id", cols)],
             by = c("arm", "subset_id"), suffixes = c("_paper", "_rerun"))
  gap <- vapply(cols, function(k) max(abs(m[[paste0(k, "_paper")]] - m[[paste0(k, "_rerun")]]), na.rm = TRUE), numeric(1))
  cat(sprintf("\n%s (arms: %s): largest absolute difference from the paper, per quantity:\n",
              key, paste(sort(unique(old$arm)), collapse = ", ")))
  print(signif(gap, 3))
  nrow(m) == nrow(old) && all(gap < EXAMPLE_TOL)
}

dat <- readRDS(path_modelling(EXAMPLE_CODE))

#1. main lattice: IPF-LASSO and offset arms (the cv.R defaults)
main_key <- paste0(EXAMPLE_CODE, "_example")
clear_example(main_key)
run_cancer(dat, EXAMPLE_CODE, key = main_key)
summarise_cancer(main_key)
ok_main <- compare_example(main_key, file.path("example", paste0(EXAMPLE_CODE, "_expected_subset_summary.csv")))

#2. IPF-ridge and RNA-score arms, IPF-LASSO and offset switched off, settings restored afterwards
rs_key <- paste0(EXAMPLE_CODE, "_example_rs")
clear_example(rs_key)
switches <- c(RUN_FUSED_ARM = FALSE, RUN_OFFSET_ARM = FALSE, RUN_FOREST_ARM = FALSE,
              RUN_COXBOOST_ARM = FALSE, RUN_RIDGE_ARM = TRUE, RUN_SCORE_ARM = TRUE)
old_switches <- mget(names(switches), envir = .GlobalEnv)
for (n in names(switches)) assign(n, switches[[n]], envir = .GlobalEnv)
run_cancer(dat, EXAMPLE_CODE, key = rs_key)
for (n in names(old_switches)) assign(n, old_switches[[n]], envir = .GlobalEnv)
summarise_cancer(rs_key)
ok_rs <- compare_example(rs_key, file.path("example", paste0(EXAMPLE_CODE, "_expected_rs_subset_summary.csv")))

if (ok_main && ok_rs) {
  cat("PASS: the rerun reproduces the paper's", EXAMPLE_CODE, "results (IPF-LASSO, offset, IPF-ridge, RNA score).\n")
} else {
  stop("DIFFERENT from the paper's ", EXAMPLE_CODE, " results; check that the package versions match renv.lock (renv::status()).")
}
