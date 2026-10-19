#Example run
#Reruns TCGA-LAML from data/processed/LAML_modelling.rds and compares it with example/LAML_expected_subset_summary.csv.
#Outputs use the key LAML_example. Run from the project root: source("scripts/run_example.R")

source("master.R")
source(file.path(DIR_SCRIPTS, "02_summarise.R"))

EXAMPLE_CODE <- "LAML"
EXAMPLE_KEY  <- "LAML_example"
EXPECTED     <- file.path("example", "LAML_expected_subset_summary.csv")
EXAMPLE_TOL  <- 1e-6  #largest allowed absolute difference

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
  cat("PASS: the rerun reproduces the paper's LAML results.\n")
} else {
  stop("DIFFERENT from the paper's LAML results; check that the package versions match renv.lock (renv::status()).")
}
