#Example run
#Reruns one cancer (LAML or GBM) from data/processed/<CODE>_modelling.rds and compares it with example/<CODE>_expected_subset_summary.csv.
#Outputs use the key <CODE>_example. Run from the project root: Rscript scripts/run_example.R LAML (or GBM);
#source("scripts/run_example.R") runs LAML.

source("master.R")
source(file.path(DIR_SCRIPTS, "02_summarise.R"))

EXAMPLE_CODE <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(EXAMPLE_CODE)) EXAMPLE_CODE <- "LAML"
if (!EXAMPLE_CODE %in% c("LAML", "GBM")) stop("the example runs LAML or GBM")
EXAMPLE_KEY  <- paste0(EXAMPLE_CODE, "_example")
EXPECTED     <- file.path("example", paste0(EXAMPLE_CODE, "_expected_subset_summary.csv"))
EXAMPLE_TOL  <- 1e-3  #largest allowed absolute difference; results are reported to three decimals

#remove outputs and checkpoints of an earlier example run
unlink(c(path_results(EXAMPLE_KEY), path_summary(EXAMPLE_KEY), path_timing(EXAMPLE_KEY), path_manifest(EXAMPLE_KEY),
         sub("_results\\.csv$", "_genes.csv.gz", path_results(EXAMPLE_KEY)),
         sub("_results\\.csv$", "_preds.rds", path_results(EXAMPLE_KEY)),
         list.files(DIR_CHECKPOINTS, pattern = paste0("^", EXAMPLE_KEY, "_"), full.names = TRUE)))

dat <- readRDS(path_modelling(EXAMPLE_CODE))
run_cancer(dat, EXAMPLE_CODE, key = EXAMPLE_KEY)
summarise_cancer(EXAMPLE_KEY)

#compare with the paper, per model and clinical subset
new  <- read.csv(path_summary(EXAMPLE_KEY), stringsAsFactors = FALSE)
old  <- read.csv(EXPECTED, stringsAsFactors = FALSE)
cols <- c("c_ref", "c_model", "delta_c", "ibs_ref", "ibs_model", "delta_ibs")
m <- merge(old[, c("arm", "subset_id", cols)], new[, c("arm", "subset_id", cols)],
           by = c("arm", "subset_id"), suffixes = c("_paper", "_rerun"))
gap <- vapply(cols, function(k) max(abs(m[[paste0(k, "_paper")]] - m[[paste0(k, "_rerun")]]), na.rm = TRUE), numeric(1))
cat("\nlargest absolute difference from the paper, per quantity:\n"); print(signif(gap, 3))
if (nrow(m) == nrow(old) && all(gap < EXAMPLE_TOL)) {
  cat("PASS: the rerun reproduces the paper's", EXAMPLE_CODE, "results.\n")
} else {
  stop("DIFFERENT from the paper's ", EXAMPLE_CODE, " results; check that the package versions match renv.lock (renv::status()).")
}
