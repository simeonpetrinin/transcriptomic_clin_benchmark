#Sweep 2: PFI endpoint
#Re-runs the main lattice with progression-free interval instead of overall survival
#(LAML skipped, no PFI). Output key <CODE>_pfi; defines only.
SWEEP_PFI_EXCLUDE <- c("LAML")

build_pfi_data <- function(code, force = FALSE) {
  key <- paste0(code, "_pfi")
  if (!force && file.exists(path_modelling(key))) return(invisible(readRDS(path_modelling(key))))
  if (!file.exists(path_modelling(code))) build_modelling_data(code)
  dat <- readRDS(path_modelling(code))
  cdr <- as.data.frame(readxl::read_xlsx(CDR_PATH, sheet = "TCGA-CDR", guess_max = 11000,
                                         .name_repair = "unique_quiet"))
  cdr <- cdr[cdr$type == code, c("bcr_patient_barcode", "PFI", "PFI.time")]
  i <- match(dat$clinical$patient_id, cdr$bcr_patient_barcode)
  n_os <- sum(dat$clinical$OS == 1, na.rm = TRUE)
  dat$clinical$OS      <- suppressWarnings(as.integer(cdr$PFI[i]))
  dat$clinical$OS.time <- suppressWarnings(as.numeric(cdr$PFI.time[i]))
  dat$endpoint <- "PFI"
  cat(sprintf("[sweep pfi] %s: %d patients | OS events %d -> PFI events %d | PFI missing %d\n",
              code, nrow(dat$clinical), n_os, sum(dat$clinical$OS == 1, na.rm = TRUE),
              sum(is.na(dat$clinical$OS) | is.na(dat$clinical$OS.time))))
  saveRDS(dat, path_modelling(key))
  invisible(dat)
}

sweep_pfi <- function(codes = names(CANCER_VARS)) {
  out <- list()
  for (code in setdiff(codes, SWEEP_PFI_EXCLUDE)) {
    key <- paste0(code, "_pfi")
    cat(sprintf("\n[sweep pfi] %s  (key %s)\n", code, key))
    out[[key]] <- tryCatch({
      dat <- build_pfi_data(code)
      run_cancer(dat, code, key = key); "ok"
    }, error = function(e) paste("FAILED:", conditionMessage(e)))
  }
  skipped <- intersect(codes, SWEEP_PFI_EXCLUDE)
  if (length(skipped)) cat("[sweep pfi] skipped (no PFI data):", paste(skipped, collapse = ", "), "\n")
  invisible(unlist(out))
}