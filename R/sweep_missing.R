#Sweep 3: missing-data handling
#Re-runs the extended specification with CLINICAL_MISSING = "impute" (training-fold
#median/mode, no indicators) instead of the primary missing-as-category, with the same
#repeats as the extended run. Output key <CODE>_ext_impute; defines only.

sweep_missing <- function(codes = names(EXT_VARS)) {
  old <- CLINICAL_MISSING
  on.exit(CLINICAL_MISSING <<- old, add = TRUE)
  CLINICAL_MISSING <<- "impute"
  out <- list()
  for (code in intersect(codes, names(EXT_VARS))) {
    key <- paste0(code, "_ext_impute")
    cat(sprintf("\n[sweep missing] %s  CLINICAL_MISSING = impute  (key %s)\n", code, key))
    out[[key]] <- tryCatch({
      #rebuild the extended object if it is missing or lacks a variable EXT_VARS now asks for
      f <- path_modelling(ext_key(code))
      dat <- if (file.exists(f)) readRDS(f) else NULL
      if (is.null(dat) || !all(extended_vars(code) %in% names(dat$clinical)))
        dat <- build_modelling_data(code, vars = extended_vars(code), key = ext_key(code))
      run_cancer(dat, code, subsets = extended_specs(code), key = key,
                 n_repeats = EXT_REPEATS); "ok"
    }, error = function(e) paste("FAILED:", conditionMessage(e)))
  }
  invisible(unlist(out))
}