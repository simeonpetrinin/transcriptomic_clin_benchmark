#External cohorts: loads the main pipeline and the external-cohort functions. Source from the project root.
source("master.R")
source(file.path(DIR_R, "external_data.R"))
source(file.path(DIR_R, "external_validation.R"))

cat("Parallel independent-validation extension loaded.\n")
cat("  run_external_cohort(\"GSE68465\")\n")
cat("  run_external_panel()\n")
cat("  compare_parallel_validation()\n")
