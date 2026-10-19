#Run sweeps
#Entry point for the three sensitivity sweeps (gene count, PFI endpoint, missing data).
#Set CODES and WHICH below; source from the project root.
source("paths.R")
source_library()

CODES <- names(CANCER_VARS)  #e.g. c("HNSC", "ESCA") on one machine
WHICH <- c("genes", "pfi", "missing")  #any subset, e.g. "pfi"

run_sweeps(codes = CODES, which = WHICH)