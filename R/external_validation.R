# PARALLEL EXTERNAL VALIDATION
# Runs each cohort with the same design as the TCGA cohorts, through cv.R, etc

EXTERNAL_MIN_COMPLETENESS <- 0.50  #as for TCGA
EXTERNAL_VALIDATION_DESIGN <- "parallel_independent_repeated_nested_cv_v1"
AUTOMATIC_EXTERNAL_COHORTS <- c("GSE68465", "GSE39582", "GSE65858", "GSE13507",
                                "GSE84437", "GSE14520")
# Manual download first (see external_data.R): CGGA_LGG, CGGA_GBM, BEATAML
#the panel reported in the paper, run by default (SCANB, METABRIC and EMTAB1980 remain callable by name)
EXTERNAL_PANEL <- c("GSE68465", "GSE39582", "GSE65858", "GSE13507", "GSE14520", "GSE84437",
                    "CGGA_LGG", "CGGA_GBM", "BEATAML")
EXTERNAL_MAX_TOP_SHARE <- 0.95  #largest category (MISSING counted as a level) at most 95%, as for TCGA

#variable type as stored in the external cohort (TCGA's types are not imposed)
.external_var_type <- function(v, x, source_types) if (is.numeric(x)) "numeric" else "categorical"

#rare-level rule: a category level held by fewer than MIN_LEVEL_N patients -> missing (never merged)
external_rare_levels_to_na <- function(clinical, cohort_id, n_min = MIN_LEVEL_N) {
  sp <- EXTERNAL_COHORTS[[cohort_id]]
  source_types <- CANCER_VARS[[sp$cancer]]
  log <- list()
  for (v in intersect(sp$variables, names(clinical))) {
    x <- clinical[[v]]
    if (.external_var_type(v, x, source_types) != "categorical") next
    tb <- table(x, useNA = "no")
    rare <- names(tb)[tb < n_min]
    if (!length(rare)) next
    clinical[[v]][as.character(x) %in% rare] <- NA
    log[[v]] <- data.frame(cohort = cohort_id, variable = v, level = rare,
                           n = as.integer(tb[rare]), action = "set to missing",
                           stringsAsFactors = FALSE)
    cat(sprintf("  [rare -> NA] %s$%s: %s\n", cohort_id, v,
                paste(sprintf("%s(%d)", rare, as.integer(tb[rare])), collapse = ", ")))
  }
  log <- if (length(log)) do.call(rbind, log) else
    data.frame(cohort = character(0), variable = character(0), level = character(0),
               n = integer(0), action = character(0))
  list(clinical = clinical, log = log)
}

.normalise_external_clinical <- function(d) {
  for (v in names(d)) {
    if (v == "patient_id") next
    if (is.numeric(d[[v]])) {
      d[[v]][!is.finite(d[[v]])] <- NA_real_
    } else {
      z <- trimws(as.character(mark_na(d[[v]])))
      z <- tolower(z)
      z[z %in% c("", ".", "na", "n/a", "unknown", "not available",
                 "not reported", "not applicable", "--")] <- NA_character_  #"." = blank in GSE14520
      d[[v]] <- z
    }
  }
  d
}

audit_external_variables <- function(external_clinical, cohort_id,
                                     min_complete = EXTERNAL_MIN_COMPLETENESS,
                                     max_top_share = EXTERNAL_MAX_TOP_SHARE) {
  sp <- EXTERNAL_COHORTS[[cohort_id]]
  if (is.null(sp)) stop("Unknown external cohort: ", cohort_id)
  source_types <- CANCER_VARS[[sp$cancer]]
  variable_min <- sp$variable_min_complete
  variable_top <- sp$variable_max_top_share  #per-variable exception to the 95% rule
  as_cat <- identical(CLINICAL_MISSING, "category")
  
  rows <- lapply(sp$variables, function(v) {
    required_complete <- min_complete
    if (!is.null(variable_min) && v %in% names(variable_min))
      required_complete <- unname(variable_min[[v]])
    
    required_top <- max_top_share
    if (!is.null(variable_top) && v %in% names(variable_top)) required_top <- unname(variable_top[[v]])
    present <- v %in% names(external_clinical)
    x <- if (present) external_clinical[[v]] else NULL
    complete <- if (present) mean(!is.na(x)) else NA_real_
    type <- .external_var_type(v, x, source_types)
    n_na <- if (present) sum(is.na(x)) else NA_integer_
    #levels as the model sees them: a missing group of MIN_LEVEL_N or more is its own level,
    #a smaller one is imputed with the most common level
    counts <- if (present && type == "categorical") {
      tb <- table(x, useNA = "no")
      if (as_cat && n_na >= MIN_LEVEL_N) c(tb, MISSING = n_na) else {
        if (length(tb)) tb[which.max(tb)] <- tb[which.max(tb)] + n_na
        tb
      }
    } else NULL
    levels <- if (!present) NA_integer_ else if (type == "categorical") length(counts) else
      length(unique(x[!is.na(x)]))
    top_share <- if (length(counts)) max(counts) / sum(counts) else NA_real_
    reason <- if (!present) {
      "absent from external cohort"
    } else if (complete < required_complete) {
      sprintf("external completeness %.1f%% < %.0f%%", 100 * complete,
              100 * required_complete)
    } else if (levels < 2L) {
      "constant in external cohort"
    } else if (is.finite(top_share) && top_share > required_top) {
      sprintf("largest level %.1f%% > %.0f%%", 100 * top_share, 100 * required_top)
    } else {
      "included"
    }
    
    data.frame(
      variable = v, type = type, external_complete = complete,
      required_external_complete = required_complete,
      external_levels = levels, top_share = top_share,
      included = identical(reason, "included"), reason = reason,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.external_variable_types <- function(audit) {
  z <- audit[audit$included, c("variable", "type"), drop = FALSE]
  if (!nrow(z)) return(list())
  out <- as.list(z$type)
  names(out) <- z$variable
  out
}

.external_subset_index <- function(vars) {
  real <- all_subsets(vars)
  rbind(
    data.frame(subset_id = 0L, n_vars = 0L, variables = ANCHOR_LABEL,
               stringsAsFactors = FALSE),
    data.frame(
      subset_id = seq_along(real),
      n_vars = vapply(real, length, integer(1)),
      variables = vapply(real, paste, character(1), collapse = "+"),
      stringsAsFactors = FALSE
    )
  )
}

.external_design_fingerprint <- function(cohort_id, vars, n_genes) {
  sp <- EXTERNAL_COHORTS[[cohort_id]]
  paste(c(
    paste0("design=", EXTERNAL_VALIDATION_DESIGN),
    paste0("cohort=", cohort_id), paste0("tcga_comparator=", sp$cancer),
    paste0("variables=", paste(vars, collapse = ",")),
    paste0("n_genes=", n_genes),
    paste0("repeats_rule=subsets>=", REPEATS_SUBSET_THRESHOLD, "->", N_REPEATS_LARGE,
           " else ", N_REPEATS_SMALL, " (own lattice: ", 2^length(vars) - 1, " subsets)"),
    paste0("outer_folds=", OUTER_FOLDS), paste0("inner_folds=", INNER_FOLDS),
    paste0("seed=", SEED), paste0("fused=", RUN_FUSED_ARM),
    paste0("offset=", RUN_OFFSET_ARM),
    paste0("fused_fallback=", FUSED_FALLBACK, " short=", FALLBACK_SHORT, " start_tol=", FALLBACK_START_TOL),
    paste0("rare_levels=<", MIN_LEVEL_N, "->missing; top_share<=", EXTERNAL_MAX_TOP_SHARE)
  ), collapse = " | ")
}

.clear_external_cv_run <- function(cohort_id) {
  repeat_files <- list.files(
    DIR_EXTERNAL_CHECKPOINTS,
    pattern = paste0("^", cohort_id, "_repeat_[0-9]+\\.rds$"),
    full.names = TRUE
  )
  targets <- c(
    repeat_files, path_external_settings(cohort_id),
    path_external_design_guard(cohort_id), path_external_results(cohort_id),
    path_external_predictions(cohort_id), path_external_genes(cohort_id),
    path_external_manifest(cohort_id), path_external_variable_audit(cohort_id),
    path_external_level_fixes(cohort_id),
    path_external_subsets(cohort_id),
    # Remove files produced by the superseded transport-validation extension.
    file.path(DIR_EXTERNAL_RESULTS, paste0(cohort_id, "_checkpoint.rds")),
    file.path(DIR_EXTERNAL_RESULTS, paste0(cohort_id, "_gene_audit.csv")),
    file.path(DIR_EXTERNAL_RESULTS, paste0(cohort_id, "_selected_genes.csv.gz")),
    file.path(DIR_EXTERNAL_RESULTS, paste0(cohort_id, "_predictions.rds"))
  )
  targets <- unique(targets[file.exists(targets)])
  if (length(targets)) unlink(targets)
  invisible(targets)
}

.with_external_cv_context <- function(cohort_id, variable_types, dat, n_genes) {
  sp <- EXTERNAL_COHORTS[[cohort_id]]
  engine_env <- environment(run_cancer)
  # Repeats follow the cohort's OWN lattice size (n_repeats_for reads
  # CANCER_VARS[[cohort_id]], set below) -- the same compute rule as TCGA.
  keys <- c("DIR_RESULTS", "DIR_CHECKPOINTS", "CANCER_VARS")
  old <- mget(keys, envir = engine_env, inherits = FALSE)
  on.exit({
    for (k in keys) assign(k, old[[k]], envir = engine_env)
  }, add = TRUE)
  
  cv_vars <- old$CANCER_VARS
  cv_vars[[cohort_id]] <- variable_types
  
  assign("DIR_RESULTS", DIR_EXTERNAL_RESULTS, envir = engine_env)
  assign("DIR_CHECKPOINTS", DIR_EXTERNAL_CHECKPOINTS, envir = engine_env)
  assign("CANCER_VARS", cv_vars, envir = engine_env)
  
  run_cancer(dat, cohort_id, n_genes = n_genes)
}

.is_complete_saved_external_run <- function(res, subset_index) {
  required <- c("subset_id", "repeat_id", "fold", "c_clin", "c_fused",
                "delta_c", "c_offset", "delta_c_offset")
  if (!is.data.frame(res) || !all(required %in% names(res)) || !nrow(res))
    return(FALSE)
  if (anyNA(res[c("subset_id", "repeat_id", "fold")])) return(FALSE)
  if (anyDuplicated(res[c("subset_id", "repeat_id", "fold")])) return(FALSE)
  
  expected_subsets <- sort(unique(subset_index$subset_id))
  observed_subsets <- sort(unique(res$subset_id))
  repeats <- sort(unique(res$repeat_id))
  folds <- sort(unique(res$fold))
  if (!identical(observed_subsets, expected_subsets) ||
      !identical(folds, seq_len(OUTER_FOLDS)) ||
      !identical(repeats, seq_len(length(repeats)))) return(FALSE)
  
  nrow(res) == length(expected_subsets) * length(repeats) * OUTER_FOLDS
}

.run_external_with_saved_result_recovery <- function(cohort_id, variable_types,
                                                     dat, n_genes, subset_index) {
  tryCatch(
    .with_external_cv_context(cohort_id, variable_types, dat, n_genes),
    error = function(e) {
      # run_cancer() writes its CSV before printing final diagnostics and before
      # its cluster-cleanup on.exit handler runs.  Recover only the specific
      # teardown failure and only when the saved result grid is demonstrably
      # complete; all genuine fitting/data errors still propagate unchanged.
      if (!identical(conditionMessage(e), "subscript out of bounds")) stop(e)
      saved_path <- path_external_results(cohort_id)
      if (!file.exists(saved_path)) stop(e)
      saved <- tryCatch(
        read.csv(saved_path, stringsAsFactors = FALSE, check.names = FALSE),
        error = function(...) NULL
      )
      if (is.null(saved) || !.is_complete_saved_external_run(saved, subset_index))
        stop(e)
      warning(cohort_id,
              ": recovered a complete saved result after worker-cluster cleanup failed.",
              call. = FALSE)
      saved
    }
  )
}

run_external_cohort <- function(cohort_id, n_genes = TOP_N_GENES,
                                force = FALSE, prepare_missing = TRUE) {
  sp <- EXTERNAL_COHORTS[[cohort_id]]
  if (is.null(sp)) stop("Unknown external cohort: ", cohort_id)
  modelling_path <- path_external_modelling(cohort_id)
  if (!file.exists(modelling_path) && prepare_missing)
    prepare_external_cohort(cohort_id)
  if (!file.exists(modelling_path))
    stop("Prepare ", cohort_id, " first: prepare_external_cohort(\"", cohort_id, "\")")
  if (force) .clear_external_cv_run(cohort_id)
  
  dat <- readRDS(modelling_path)
  validate_external_object(dat, cohort_id)
  dat$clinical <- .normalise_external_clinical(dat$clinical)
  fx <- external_rare_levels_to_na(dat$clinical, cohort_id)  #modelling cohort: survival filter applied at prepare
  dat$clinical <- fx$clinical
  write.csv(fx$log, path_external_level_fixes(cohort_id), row.names = FALSE)
  
  audit <- audit_external_variables(dat$clinical, cohort_id)
  write.csv(audit, path_external_variable_audit(cohort_id), row.names = FALSE)
  variable_types <- .external_variable_types(audit)
  vars <- names(variable_types)
  if (!length(vars)) stop(cohort_id, ": no clinical variables passed the external audit")
  
  idx <- .external_subset_index(vars)
  idx$cohort <- cohort_id
  idx$tcga_comparator <- sp$cancer
  write.csv(idx, path_external_subsets(cohort_id), row.names = FALSE)
  
  design_fp <- .external_design_fingerprint(cohort_id, vars, n_genes)
  guard <- path_external_design_guard(cohort_id)
  if (file.exists(guard)) {
    old <- paste(readLines(guard, warn = FALSE), collapse = "")
    if (!identical(old, design_fp))
      stop("Parallel-validation settings changed for ", cohort_id,
           ". Re-run with force = TRUE.\nwas: ", old, "\nnow: ", design_fp)
  } else {
    writeLines(design_fp, guard)
  }
  
  cat(sprintf(
    "=== parallel external validation: %s (%s comparator) | %d patients | %d variables | %d specifications ===\n",
    cohort_id, sp$cancer, nrow(dat$clinical), length(vars), nrow(idx)
  ))
  res <- .run_external_with_saved_result_recovery(
    cohort_id, variable_types, dat, n_genes, idx
  )
  
  # Preserve the run_cancer schema and add explicit replication metadata.
  res <- cbind(
    validation_design = EXTERNAL_VALIDATION_DESIGN,
    cohort = cohort_id, tcga_comparator = sp$cancer,
    res
  )
  write.csv(res, path_external_results(cohort_id), row.names = FALSE)
  
  write(c(
    "", paste0("validation_design  : ", EXTERNAL_VALIDATION_DESIGN),
    paste0("cohort_citation   : ", sp$citation),
    paste0("tcga_comparator   : ", sp$cancer),
    "replication_rule  : all preprocessing, feature selection, tuning and fitting use this cohort's training folds only",
    "cross_dataset_fit : none; compare independently estimated specification curves"
  ), file = path_external_manifest(cohort_id), append = TRUE)
  res
}

run_external_panel <- function(cohort_ids = EXTERNAL_PANEL,
                               n_genes = TOP_N_GENES, force = FALSE,
                               prepare_missing = TRUE) {
  failed <- character(0)
  for (id in cohort_ids) {
    ok <- tryCatch({
      run_external_cohort(id, n_genes = n_genes, force = force,
                          prepare_missing = prepare_missing)
      TRUE
    }, error = function(e) {
      message("[parallel external validation FAILED] ", id, ": ", conditionMessage(e))
      FALSE
    })
    if (!ok) failed <- c(failed, id)
  }
  completed <- setdiff(cohort_ids, failed)
  cat(sprintf("Parallel external panel complete: %d succeeded, %d failed.\n",
              length(completed), length(failed)))
  invisible(list(failed = failed, completed = completed))
}

.mean_or_na <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
.median_or_na <- function(x) if (any(is.finite(x))) median(x[is.finite(x)]) else NA_real_
.slope_or_na <- function(y, x) {
  ok <- is.finite(y) & is.finite(x)
  if (sum(ok) < 2L || length(unique(x[ok])) < 2L) return(NA_real_)
  unname(coef(lm(y[ok] ~ x[ok]))[2L])
}

.summarise_independent_curve <- function(d, dataset, cohort, cancer) {
  z <- d[!d$is_anchor, , drop = FALSE]
  if (!nrow(z)) stop("No non-anchor specifications in ", cohort)
  by_spec <- split(z, z$subset_id)
  spec <- do.call(rbind, lapply(by_spec, function(q) data.frame(
    subset_id = q$subset_id[1L], n_vars = q$n_vars[1L],
    delta_c = .mean_or_na(q$delta_c), delta_ibs = .mean_or_na(q$delta_ibs),
    delta_c_offset = .mean_or_na(q$delta_c_offset),
    delta_ibs_offset = .mean_or_na(q$delta_ibs_offset),
    stringsAsFactors = FALSE
  )))
  data.frame(
    dataset = dataset, cohort = cohort, cancer = cancer,
    patients = max(d$n_train + d$n_test, na.rm = TRUE),
    specifications = nrow(spec), fold_rows = nrow(z),
    median_delta_c = .median_or_na(spec$delta_c),
    delta_c_slope_per_variable = .slope_or_na(spec$delta_c, spec$n_vars),
    median_delta_c_offset = .median_or_na(spec$delta_c_offset),
    delta_c_offset_slope_per_variable = .slope_or_na(spec$delta_c_offset, spec$n_vars),
    median_delta_ibs = .median_or_na(spec$delta_ibs),
    median_delta_ibs_offset = .median_or_na(spec$delta_ibs_offset),
    stringsAsFactors = FALSE
  )
}

summarise_external_panel <- function(cohort_ids = EXTERNAL_PANEL) {
  out <- list()
  for (id in cohort_ids) {
    p <- path_external_results(id)
    if (!file.exists(p)) next
    d <- read.csv(p, stringsAsFactors = FALSE)
    sp <- EXTERNAL_COHORTS[[id]]
    out[[id]] <- .summarise_independent_curve(d, "external", id, sp$cancer)
  }
  if (!length(out)) stop("No parallel external result files found")
  ans <- do.call(rbind, out)
  rownames(ans) <- NULL
  write.csv(ans, file.path(DIR_EXTERNAL_RESULTS, "external_panel_summary.csv"),
            row.names = FALSE)
  ans
}

compare_parallel_validation <- function(cohort_ids = EXTERNAL_PANEL) {
  ext <- summarise_external_panel(cohort_ids)
  tcga <- list()
  for (code in unique(ext$cancer)) {
    p <- path_results(code)
    if (!file.exists(p)) next
    d <- read.csv(p, stringsAsFactors = FALSE)
    tcga[[code]] <- .summarise_independent_curve(d, "TCGA", code, code)
  }
  ans <- rbind(if (length(tcga)) do.call(rbind, tcga) else NULL, ext)
  rownames(ans) <- NULL
  write.csv(ans, file.path(DIR_EXTERNAL_RESULTS, "parallel_comparison_summary.csv"),
            row.names = FALSE)
  ans
}