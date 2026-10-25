#Sweep 1: RNA gene count
#Re-runs the main lattice with other top-variance gene counts (SWEEP_GENE_COUNTS).
#Inf = every gene that passes the training completeness check (build_rna takes the minimum).
#Output key <CODE>_g<N>, or <CODE>_gall for Inf; defines only.
SWEEP_GENE_COUNTS <- c(100, Inf)

#file-key label for a gene count: 100 -> "g100", Inf -> "gall"
gene_count_label <- function(n) if (is.infinite(n)) "gall" else sprintf("g%d", as.integer(n))

#full lattice for the chosen cancers
sweep_genes <- function(codes = names(CANCER_VARS), counts = SWEEP_GENE_COUNTS) {
  out <- list()
  for (code in codes) for (n in counts) {
    key <- paste0(code, "_", gene_count_label(n))
    cat(sprintf("\n[sweep genes] %s  n_genes = %s  (key %s)\n", code, format(n), key))
    if (!file.exists(path_modelling(code))) build_modelling_data(code)
    dat <- readRDS(path_modelling(code))
    out[[key]] <- tryCatch({ run_cancer(dat, code, n_genes = n, key = key); "ok" },
                           error = function(e) paste("FAILED:", conditionMessage(e)))
  }
  invisible(unlist(out))
}

#the robustness sample the gene sweep was run on (omics only, age, main set, 10 sampled), drawn before
#the robustness run grew to 80 sampled subsets (results/robustness_specs.csv)
GENE_SWEEP_SPECS_FILE <- file.path(DIR_RESULTS, "robustness_specs_gene_sweep.csv")

#robustness sample only (GENE_SWEEP_SPECS_FILE); IPF-LASSO and offset arms
#Output key <CODE>_robust_g<N> / <CODE>_robust_gall
sweep_genes_robust <- function(codes = names(CANCER_VARS), counts = SWEEP_GENE_COUNTS) {
  if (!exists(".saved_robust_specs", mode = "function"))
    stop(".saved_robust_specs() is not loaded: source master.R (it sources scripts/run_robustness.R)")
  #sensitivity learners off for the whole sweep, previous settings restored on exit
  flags <- c("RUN_FOREST_ARM", "RUN_COXBOOST_ARM", "RUN_RIDGE_ARM", "RUN_SCORE_ARM")
  old <- mget(flags, envir = .GlobalEnv)
  on.exit(for (f in flags) assign(f, old[[f]], envir = .GlobalEnv), add = TRUE)
  for (f in flags) assign(f, FALSE, envir = .GlobalEnv)
  out <- list()
  for (code in codes) for (n in counts) {
    key <- paste0(code, "_robust_", gene_count_label(n))
    cat(sprintf("\n[sweep genes, robustness sample] %s  n_genes = %s  (key %s)\n", code, format(n), key))
    if (!file.exists(path_modelling(code))) build_modelling_data(code)
    dat <- readRDS(path_modelling(code))
    out[[key]] <- tryCatch({
      specs <- .saved_robust_specs(code, GENE_SWEEP_SPECS_FILE)
      if (is.null(specs)) stop("no gene-sweep specs for ", code, " in ", GENE_SWEEP_SPECS_FILE)
      run_cancer(dat, code, n_genes = n, subsets = specs, key = key); "ok" },
                           error = function(e) paste("FAILED:", conditionMessage(e)))
  }
  invisible(unlist(out))
}