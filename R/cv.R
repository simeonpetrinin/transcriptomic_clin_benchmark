#CV
#Nested cross-validation for one cancer: the clinical Cox model and every switched-on
#clinical + RNA model over each clinical subset, on shared outer folds, scored by
#Uno's C, the IBS and D-calibration per fold. Defines only; run_cancer() is the entry.


suppressPackageStartupMessages(library(survival))

OUTER_FOLDS      <- 5
N_REPEATS        <- 5  #NOT the count used: runs use n_repeats_for(code), below
GRID_QUANTILES   <- seq(0.1, 0.9, length.out = 25)  #IBS grid: q10-q90 of TRAIN event times

#OFFSET ARM (Volkmann et al. 2019)
RUN_OFFSET_ARM   <- TRUE  #sensitivity analysis (the fused arm is primary)
SAVE_OFFSET_SURV <- FALSE  #survival matrices for the offset arm as well

#arm switches (the clinical Cox always runs; an arm switched off writes NA columns)
RUN_FUSED_ARM    <- TRUE

#forest arm switch
RUN_FOREST_ARM     <- FALSE  #OFF for the lattice and extended runs (too slow)
FOREST_SEED_OFFSET <- 0L
RUN_RIDGE_ARM      <- FALSE  #ridge fused arm (robustness runs only)
RUN_SCORE_ARM      <- FALSE  #RNA-score arm (robustness runs only)
RUN_COXBOOST_ARM   <- FALSE  #CoxBoost favoring clinical (robustness runs only)
SAVE_FOREST_PREDS  <- FALSE  #forest risk + survival per tree count in _preds.rds

#repeats: compute rule after Herrmann et al. (2021), keyed on lattice size here
N_REPEATS_LARGE  <- 5L
N_REPEATS_SMALL  <- 10L
REPEATS_SUBSET_THRESHOLD <- 1023L  #2^10 - 1
n_subsets_for <- function(code) {
  v <- CANCER_VARS[[code]]
  if (is.null(v)) stop("no CANCER_VARS entry for ", code)
  2^length(v) - 1
}
n_repeats_for <- function(code)
  if (n_subsets_for(code) >= REPEATS_SUBSET_THRESHOLD) N_REPEATS_LARGE else N_REPEATS_SMALL
SEED             <- 42  #master seed
MIN_TRAIN_EVENTS <- 5  #skip a fold whose TRAIN side is event-starved
if (!exists("N_WORKERS")) N_WORKERS <- 1L  #set in settings.txt (read by paths.R); 1 = sequential
SAVE_GENE_NAMES  <- TRUE  #write <CODE>_genes.csv.gz
SAVE_PREDICTIONS <- TRUE  #write <CODE>_preds.rds (see below)

#predictions are saved so a metric can be fixed or added later without refitting

HOST <- tryCatch(as.character(Sys.info()[["nodename"]]), error = function(e) "unknown")

#stratified outer folds on the event indicator
make_outer_folds <- function(status, k = OUTER_FOLDS, seed = SEED) {
  set.seed(seed)
  fold <- integer(length(status))
  for (lev in unique(status)) {  #stratify: split each class separately
    idx <- which(status == lev)
    fold[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
  }
  fold
}

#deterministic seed for ONE outer-training set
inner_fold_seed <- function(r, f) SEED + 1000L * as.integer(r) +
  100L * as.integer(f)

#deterministic seed for ONE cell's forests
forest_seed <- function(r, f, sid, offset = FOREST_SEED_OFFSET)
  as.integer(SEED + 1e8 * as.integer(offset) + 1e6 * as.integer(r) +
               1e5 * as.integer(f) + as.integer(sid))

timed <- function(expr) {
  t0  <- proc.time()[["elapsed"]]
  val <- force(expr)
  list(value = val, secs = round(proc.time()[["elapsed"]] - t0, 4))
}

#null-coalesce for the result row.
nz <- function(x, alt) if (is.null(x) || length(x) == 0) alt else x

#D-calibration needs survival on a fine grid (every training event time)
.split_surv <- function(m, grid, dgrid) {
  if (is.null(m) || !is.list(m)) return(m)
  pg  <- m$times
  cut <- function(p) {
    if (!is.matrix(p$surv) || is.null(pg) || ncol(p$surv) != length(pg)) return(p)
    p$surv_d <- p$surv[, match(dgrid, pg), drop = FALSE]
    p$surv   <- p$surv[, match(grid,  pg), drop = FALSE]
    p
  }
  m <- cut(m)
  if (!is.null(m$by_trees)) m$by_trees <- lapply(m$by_trees, cut)
  m$times <- grid; m$dtimes <- dgrid
  m
}

#D-calibration of one model's test-fold curves (NA if the fit failed)
.dcal <- function(surv_d, dgrid, tec) {
  if (!is.matrix(surv_d) || nrow(surv_d) != nrow(tec)) return(NA_real_)
  as.numeric(d_calibration(tec$OS.time, tec$OS, surv_d, dgrid))
}

#score one forest fit at one tree count -> c(C, IBS)
.score_forest <- function(m, nt, tec, trc, tau) {
  if (is.null(m) || is.null(m$by_trees)) return(c(NA_real_, NA_real_, NA_real_))
  p <- m$by_trees[[as.character(nt)]]
  if (is.null(p)) return(c(NA_real_, NA_real_, NA_real_))
  c(uno_c(tec$OS.time, tec$OS, p$risk, tau),
    as.numeric(integrated_brier(tec$OS.time, tec$OS, p$surv, m$times, tau,
                                cens_time = trc$OS.time, cens_status = trc$OS)),
    .dcal(p$surv_d, m$dtimes, tec))
}

#all forest columns for one row, identical columns on every row (normal, anchor, arm off)
forest_row_cols <- function(m4, m5, s1, is_anchor, ran, tec, trc, tau, fseed,
                            secs4, secs5, secs4_all = NA_real_, secs5_all = NA_real_,
                            not_fit_msg = NA_character_) {
  out <- list(); t_sc <- 0; t_sb <- 0
  for (nt in FOREST_TREES) {
    t0 <- .now()
    sc <- if (!ran) c(NA_real_, NA_real_, NA_real_) else
      if (is_anchor) s1 else .score_forest(m4, nt, tec, trc, tau)
    t_sc <- t_sc + (.now() - t0); t0 <- .now()
    sb <- if (!ran) c(NA_real_, NA_real_, NA_real_) else .score_forest(m5, nt, tec, trc, tau)
    t_sb <- t_sb + (.now() - t0)
    out[[paste0("c_rfclin_",     nt)]] <- sc[1]
    out[[paste0("ibs_rfclin_",   nt)]] <- sc[2]
    out[[paste0("c_bf_",         nt)]] <- sb[1]
    out[[paste0("ibs_bf_",       nt)]] <- sb[2]
    out[[paste0("delta_c_bf_",   nt)]] <- sb[1] - sc[1]
    out[[paste0("delta_ibs_bf_", nt)]] <- sb[2] - sc[2]
    out[[paste0("dcal_rfclin_",  nt)]] <- sc[3]
    out[[paste0("dcal_bf_",      nt)]] <- sb[3]
    out[[paste0("delta_dcal_bf_", nt)]] <- sb[3] - sc[3]
  }
  out <- c(out, list(
    forest_seed    = if (ran) fseed else NA_integer_,
    rfclin_mtry    = nz(m4$mtry,   NA_integer_),
    bf_w_clin      = nz(m5$w_clin, NA_real_),  #tuned block weights
    bf_w_rna       = nz(m5$w_rna,  NA_real_),
    rfclin_warning = nz(m4$warning, NA_character_),
    bf_warning     = nz(m5$warning, NA_character_),
    rfclin_error   = if (!is.na(not_fit_msg)) not_fit_msg else nz(m4$error, NA_character_),
    bf_error       = if (!is.na(not_fit_msg)) not_fit_msg else nz(m5$error, NA_character_),
    secs_rfclin    = secs4,
    secs_bf        = secs5,
    secs_rfclin_total = round(secs4_all + t_sc, 4),
    secs_bf_total     = round(secs5_all + t_sb, 4)))
  as.data.frame(out, stringsAsFactors = FALSE)
}

#boosting columns for one row (same rule: identical columns on every row)
coxboost_row_cols <- function(m6, s1, ran, tec, trc, tau, secs6, secs6_all = NA_real_,
                              not_fit_msg = NA_character_) {
  t0 <- .now()
  sb <- c(NA_real_, NA_real_, NA_real_)
  if (ran && !is.null(m6) && !(length(m6$risk) == 1 && is.na(m6$risk))) {
    sb <- c(uno_c(tec$OS.time, tec$OS, m6$risk, tau),
            as.numeric(integrated_brier(tec$OS.time, tec$OS, m6$surv, m6$times, tau,
                                        cens_time = trc$OS.time, cens_status = trc$OS)),
            .dcal(m6$surv_d, m6$dtimes, tec))
  }
  data.frame(
    c_cb = sb[1], ibs_cb = sb[2],
    delta_c_cb = if (ran) sb[1] - s1[1] else NA_real_,
    delta_ibs_cb = if (ran) sb[2] - s1[2] else NA_real_,
    dcal_cb = sb[3], delta_dcal_cb = if (ran) sb[3] - s1[3] else NA_real_,
    cb_stepno = nz(m6$stepno, NA_integer_), cb_at_boundary = nz(m6$at_boundary, NA),
    cb_n_genes = nz(m6$n_genes, NA_integer_),
    cb_warning = nz(m6$warning, NA_character_),
    cb_error = if (!is.na(not_fit_msg)) not_fit_msg else nz(m6$error, NA_character_),
    secs_cb = secs6, secs_cb_total = round(secs6_all + (.now() - t0), 4),
    stringsAsFactors = FALSE)
}

#ridge / RNA-score columns (same rule: identical columns on every row)
.simple_arm_cols <- function(prefix, m, s1, ran, tec, trc, tau, secs, extra,
                             secs_all = NA_real_, not_fit_msg = NA_character_) {
  t0 <- .now()
  sb <- c(NA_real_, NA_real_, NA_real_)
  if (ran && !is.null(m) && !(length(m$risk) == 1 && is.na(m$risk)))
    sb <- c(uno_c(tec$OS.time, tec$OS, m$risk, tau),
            as.numeric(integrated_brier(tec$OS.time, tec$OS, m$surv, m$times, tau,
                                        cens_time = trc$OS.time, cens_status = trc$OS)),
            .dcal(m$surv_d, m$dtimes, tec))
  out <- list(sb[1], sb[2], if (ran) sb[1] - s1[1] else NA_real_, if (ran) sb[2] - s1[2] else NA_real_,
              sb[3], if (ran) sb[3] - s1[3] else NA_real_)
  names(out) <- paste0(c("c_", "ibs_", "delta_c_", "delta_ibs_", "dcal_", "delta_dcal_"), prefix)
  for (nm in names(extra)) {
    col <- if (startsWith(nm, paste0(prefix, "_"))) nm else paste0(prefix, "_", nm)
    out[[col]] <- nz(m[[nm]], extra[[nm]])
  }
  out[[paste0(prefix, "_error")]] <- if (!is.na(not_fit_msg)) not_fit_msg else nz(m$error, NA_character_)
  out[[paste0("secs_", prefix)]] <- secs
  out[[paste0("secs_", prefix, "_total")]] <- round(secs_all + (.now() - t0), 4)
  as.data.frame(out, stringsAsFactors = FALSE)
}
#at_boundary = lightest end of the path chosen; at_top = heaviest (lambda_max)
RIDGE_EXTRA <- list(sel_idx = NA_integer_, at_boundary = NA, at_top = NA,
                    lambda_fraction = NA_real_)
SCORE_EXTRA <- list(score_used = NA, score_n_genes = NA_integer_, separated = NA,
                    at_boundary = NA)

#fingerprint of the settings that change results
settings_fingerprint <- function(code, n_genes, spec_desc = NULL,
                                 n_repeats = n_repeats_for(code)) {
  paste(c(sprintf("code=%s", code),
          sprintf("n_genes=%s", n_genes),
          sprintf("outer=%d repeats=%d seed=%d", OUTER_FOLDS, n_repeats, SEED),
          sprintf("min_events=%d", MIN_TRAIN_EVENTS),
          sprintf("lambda_min_ratio=%s nlambda=%s", LAMBDA_MIN_RATIO, NLAMBDA),
          sprintf("inner_folds=%s", INNER_FOLDS),
          "inner_preprocessing=nested_v1 lambda_alignment=lambda_over_lambda_max grouped_deviance",
          sprintf("cox_ties=%s rna_completeness=%s", COX_TIES, RNA_MIN_COMPLETENESS),
          #These two change EVERY IBS number.
          sprintf("ibs_grid_points=%d g_floor=%s", length(GRID_QUANTILES), G_FLOOR),
          sprintf("dcal_bins=%d dcal_grid=train_event_times", D_CAL_BINS),
          sprintf("coef_separation=%s", COEF_SEPARATION),
          sprintf("fused_fallback=%s short=%s start_tol=%s maxiter=%d", FUSED_FALLBACK, FALLBACK_SHORT,
                  FALLBACK_START_TOL, FALLBACK_MAXITER),
          sprintf("fused_arm=%s offset_arm=%s offset_nlambda=%d offset_lmr=%s",
                  RUN_FUSED_ARM, RUN_OFFSET_ARM, OFFSET_NLAMBDA, OFFSET_LAMBDA_MIN_RATIO),
          sprintf("forest_arm=%s forest_trees=%s bf_nsets=%d bf_grid=%s bf_trees_pre=%d splitrule=%s forest_seed_offset=%d bf_fix_clinical_weight=%s",
                  RUN_FOREST_ARM, paste(FOREST_TREES, collapse = "/"), BF_NSETS,
                  paste(BF_WEIGHT_GRID, collapse = "/"),
                  BF_NUM_TREES_PRE, RF_SPLITRULE, FOREST_SEED_OFFSET, BF_FIX_CLINICAL_WEIGHT),
          sprintf("coxboost_arm=%s coxboost_maxstepno=%d", RUN_COXBOOST_ARM, COXBOOST_MAXSTEPNO),
          sprintf("ridge_arm=%s ridge_lmr=%s score_arm=%s score_alpha=%s score_lmr=%s score=prevalidated",
                  RUN_RIDGE_ARM, RIDGE_LAMBDA_MIN_RATIO, RUN_SCORE_ARM, SCORE_ALPHA, SCORE_LAMBDA_MIN_RATIO),
          sprintf("clinical_missing=%s missing_level=%s reference=most_frequent", CLINICAL_MISSING, MISSING_LEVEL),
          sprintf("rare_missing=%s", paste(sort(RARE_MISSING_VARS), collapse = ",")),
          #lattice runs are defined by CANCER_VARS
          if (is.null(spec_desc))
            sprintf("vars=%s", paste(names(CANCER_VARS[[code]]), collapse = ","))
          else sprintf("specs=%s", spec_desc)),
        collapse = " | ")
}

write_manifest <- function(code, seeds, n_subsets, n_rows, total_secs, n_genes,
                           key = code) {
  gv <- tryCatch(as.character(packageVersion("glmnet")),   error = function(e) "NA")
  sv <- tryCatch(as.character(packageVersion("survival")), error = function(e) "NA")
  rv <- tryCatch(as.character(packageVersion("ranger")),      error = function(e) "NA")
  bv <- tryCatch(as.character(packageVersion("blockForest")), error = function(e) "NA")
  writeLines(c(
    sprintf("cancer            : %s%s", code, if (!identical(key, code)) sprintf("  (run: %s)", key) else ""),
    sprintf("clinical_missing  : %s (%s)", CLINICAL_MISSING,
            if (identical(CLINICAL_MISSING, "category"))
              "missing as its own level + numeric indicators"
            else "training-fold median/mode imputation, no indicators"),
    sprintf("rare_missing      : %s (missing for under %d patients: imputed, no missing column)",
            if (length(RARE_MISSING_VARS)) paste(RARE_MISSING_VARS, collapse = ", ") else "none",
            MIN_LEVEL_N),
    sprintf("timestamp         : %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    sprintf("host              : %s", HOST),
    sprintf("outer_folds       : %d", OUTER_FOLDS),
    sprintf("n_repeats         : %d  (compute rule: %d if >= %d clinical subsets, else %d; this lattice: %d%s)",
            length(seeds), N_REPEATS_LARGE, REPEATS_SUBSET_THRESHOLD, N_REPEATS_SMALL,
            as.integer(n_subsets_for(code)),
            if (length(seeds) != n_repeats_for(code)) "; set by the caller for this run" else ""),
    sprintf("offset_arm        : %s (Volkmann 2019; frozen-coefficient coxph, iter.max = 0)",
            RUN_OFFSET_ARM),
    sprintf("forest_arm        : %s (blockForest vs clinical RSF; always.select.block = clinical; clinical RSF mtry = all clinical columns)",
            RUN_FOREST_ARM),
    sprintf("forest_trees      : %s (one forest grown at max; smaller counts = first n trees)",
            paste(FOREST_TREES, collapse = ", ")),
    sprintf("blockfor_settings : %s, num.trees.pre = %d, splitrule = %s (package defaults apart from the grid)",
            if (is.null(BF_WEIGHT_GRID)) sprintf("nsets = %d random", BF_NSETS)
            else sprintf("RNA weight grid of %d (%s to %s)", length(BF_WEIGHT_GRID),
                         min(BF_WEIGHT_GRID), max(BF_WEIGHT_GRID)),
            BF_NUM_TREES_PRE, RF_SPLITRULE),
    sprintf("coxboost_arm      : %s (clinical mandatory/unpenalised; steps by Verweij CV on the fused arm's inner folds, max %d)",
            RUN_COXBOOST_ARM, COXBOOST_MAXSTEPNO),
    sprintf("bf_clinical_weight: %s", if (isTRUE(BF_FIX_CLINICAL_WEIGHT))
      "FIXED at 1 (only the RNA block weight tuned; models.R .blockfor_fixed_clinical)"
      else "tuned by blockfor() (package default)"),
    sprintf("forest_seed_offset: %d (0 = main run)", FOREST_SEED_OFFSET),
    sprintf("master_seed       : %d", SEED),
    sprintf("per_repeat_seeds  : %s", paste(seeds, collapse = ", ")),
    sprintf("inner_folds       : %d (stratified, shared within repeat x outer fold)", INNER_FOLDS),
    sprintf("top_n_genes       : %s", n_genes),
    "inner_preprocessing: fit independently within each inner-training split",
    "lambda_alignment  : lambda/lambda_max; each training split defines its own path",
    "inner_score       : grouped Cox deviance; validation-patient weighted mean",
    "inner_path_policy : candidate must be scored in every inner fold; no tail extrapolation",
    sprintf("cox_ties          : %s", COX_TIES),
    sprintf("ibs_interval      : q10-q90 of training-fold event times"),
    sprintf("ibs_grid_points   : %d requested quantiles (duplicates removed if tied)", length(GRID_QUANTILES)),
    sprintf("d_calibration     : %d bins, survival read at every training-fold event time (SurvBoard)", D_CAL_BINS),
    sprintf("uno_tau_convention: max event time per training fold (make_tau)"),
    sprintf("min_train_events  : %d", MIN_TRAIN_EVENTS),
    sprintf("n_workers         : %d", N_WORKERS),
    sprintf("omics_only_anchor : subset_id 0 (%s)", ANCHOR_LABEL),
    sprintf("n_subsets         : %d (incl. anchor)", n_subsets),
    sprintf("n_result_rows     : %d", n_rows),
    sprintf("total_cell_seconds: %.1f  (%.2f h; every model, training + prediction + scoring)", total_secs, total_secs / 3600),
    sprintf("glmnet_version    : %s", gv),
    sprintf("survival_version  : %s", sv),
    sprintf("ranger_version    : %s", rv),
    sprintf("blockForest_version: %s", bv),
    sprintf("R_version         : %s", R.version.string)
  ), path_manifest(key))
  cat("Saved ", path_manifest(key), "\n", sep = "")
}

#ONE CELL: one clinical subset within one fold.
run_one_cell <- function(subset, sid, code, r, f, seed_r,
                         trc, tec, rn, tau, grid, dgrid, rna_secs, inner_foldid, inner_rna) {
  t_cell <- .now()  #whole-cell time (secs_total)
  is_anchor <- (length(subset) == 0)
  pgrid <- sort(unique(c(grid, dgrid)))  #models predict on both grids
  label     <- if (is_anchor) ANCHOR_LABEL else paste(subset, collapse = "+")
  
  cl <- build_clinical(subset, trc, tec)  #0-col matrix when is_anchor
  
  #Do not silently drop a real subset that becomes degenerate in this fold.
  if (!is_anchor && ncol(cl$train) == 0) {
    msg <- "degenerate clinical subset: no usable columns after fold-specific preprocessing"
    row <- data.frame(
      cancer = code, subset_id = sid, n_vars = length(subset), variables = label,
      is_anchor = FALSE, repeat_id = r, seed = seed_r, fold = f,
      c_clin = NA_real_, ibs_clin = NA_real_,
      c_fused = NA_real_, ibs_fused = NA_real_,
      delta_c = NA_real_, delta_ibs = NA_real_,
      dcal_clin = NA_real_, dcal_fused = NA_real_, delta_dcal = NA_real_,
      dcal_offset = NA_real_, delta_dcal_offset = NA_real_,
      fused_n_path = NA_integer_, fused_sel_idx = NA_integer_,
      fused_path_short = NA, fused_solver_error = NA, fused_at_boundary = NA,
      fused_n_genes = NA_integer_,
      fused_lambda_fraction = NA_real_, fused_inner_n_path = NA_integer_,
      fused_inner_at_boundary = NA, fused_inner_solver_error = NA,
      fused_outer_clin_cols = NA_integer_,
      fused_inner_clin_cols_min = NA_integer_, fused_inner_clin_cols_max = NA_integer_,
      fused_inner_clin_mismatch = NA, fused_inner_vars_dropped = NA_character_,
      fused_fallback = NA_character_, fused_fallback_why = NA_character_,
      clin_n_coef = NA_integer_, clin_max_abs_coef = NA_real_,
      clin_n_bad_coef = NA_integer_, clin_n_iter = NA_integer_,
      clin_separated = NA, clin_warning = NA_character_,
      c_offset = NA_real_, ibs_offset = NA_real_,
      delta_c_offset = NA_real_, delta_ibs_offset = NA_real_,
      offset_n_path = NA_integer_, offset_sel_idx = NA_integer_,
      offset_path_short = NA, offset_solver_error = NA,
      offset_at_boundary = NA, offset_n_genes = NA_integer_,
      offset_lambda_fraction = NA_real_, offset_inner_n_path = NA_integer_,
      offset_inner_at_boundary = NA, offset_inner_solver_error = NA,
      offset_inner_clin_mismatch = NA, offset_error = paste0("not fit: ", msg),
      ibs_offset_floored = NA_integer_, secs_offset = NA_real_,
      fused_error = paste0("not fit: ", msg), clin_error = msg,
      ibs_clin_floored = NA_integer_, ibs_fused_floored = NA_integer_,
      secs_clin = 0, secs_fused = 0, secs_total = round(.now() - t_cell, 4),
      secs_clin_total = NA_real_, secs_fused_total = NA_real_, secs_offset_total = NA_real_,
      secs_rna_block = rna_secs,
      n_train = nrow(trc), n_test = nrow(tec),
      n_events_train = sum(trc$OS), n_cols_fused = NA_integer_,
      host = HOST,
      stringsAsFactors = FALSE
    )
    row <- cbind(row, forest_row_cols(
      NULL, NULL, c(NA_real_, NA_real_, NA_real_), is_anchor = FALSE, ran = FALSE,
      tec, trc, NA_real_, NA_integer_, NA_real_, NA_real_,
      not_fit_msg = paste0("not fit: ", msg)))
    row <- cbind(row, coxboost_row_cols(NULL, c(NA_real_, NA_real_, NA_real_), FALSE, tec, trc,
                                        NA_real_, NA_real_, not_fit_msg = paste0("not fit: ", msg)),
                 .simple_arm_cols("ridge", NULL, c(NA_real_, NA_real_, NA_real_), FALSE, tec, trc, NA_real_, NA_real_,
                                  RIDGE_EXTRA, not_fit_msg = paste0("not fit: ", msg)),
                 .simple_arm_cols("score", NULL, c(NA_real_, NA_real_, NA_real_), FALSE, tec, trc, NA_real_, NA_real_,
                                  SCORE_EXTRA, not_fit_msg = paste0("not fit: ", msg)))
    return(list(row = row, genes = character(0), preds = NULL))
  }
  
  #clinical-only side (timed)
  t1 <- .now()
  if (is_anchor) {
    tm1 <- timed(fit_km(trc$OS.time, trc$OS, nrow(tec), pgrid))
    m1  <- .split_surv(tm1$value, grid, dgrid)
    
    ibs1 <- integrated_brier(
      tec$OS.time, tec$OS, m1$surv, m1$times, tau,
      cens_time = trc$OS.time,
      cens_status = trc$OS
    )
    ibs1_floored <- attr(ibs1, "n_floored")
    if (is.null(ibs1_floored)) ibs1_floored <- NA_integer_
    
    s1 <- c(0.5, as.numeric(ibs1), .dcal(m1$surv_d, dgrid, tec))
    
  } else {
    tm1 <- timed(fit_clinical(cl$train, trc$OS.time, trc$OS, cl$test, pgrid))
    m1  <- .split_surv(tm1$value, grid, dgrid)
    
    if (length(m1$risk) == 1 && is.na(m1$risk)) {
      s1 <- c(NA_real_, NA_real_, NA_real_)
      ibs1_floored <- NA_integer_
    } else {
      c1 <- uno_c(tec$OS.time, tec$OS, m1$risk, tau)
      
      ibs1 <- integrated_brier(
        tec$OS.time, tec$OS, m1$surv, m1$times, tau,
        cens_time = trc$OS.time,
        cens_status = trc$OS
      )
      ibs1_floored <- attr(ibs1, "n_floored")
      if (is.null(ibs1_floored)) ibs1_floored <- NA_integer_
      
      s1 <- c(c1, as.numeric(ibs1), .dcal(m1$surv_d, dgrid, tec))
    }
  }
  
  secs1_total <- .now() - t1  #fit + prediction + scoring
  
  #clinical + RNA (block penalty).
  t2 <- .now()
  m2 <- NULL; s2 <- c(NA_real_, NA_real_, NA_real_); ibs2_floored <- NA_integer_
  tm2 <- list(secs = NA_real_)
  if (isTRUE(RUN_FUSED_ARM)) {
    fx  <- fuse(cl, rn)
    tm2 <- timed(fit_clin_rna(
      fx$train, trc$OS.time, trc$OS, fx$test,
      fx$penalty, pgrid,
      foldid = inner_foldid, clinical_subset = subset, clinical_train = trc,
      inner_rna = inner_rna
    ))
    m2 <- .split_surv(tm2$value, grid, dgrid)
    
    if (length(m2$risk) == 1 && is.na(m2$risk)) {
      s2 <- c(NA_real_, NA_real_, NA_real_)
      ibs2_floored <- NA_integer_
    }
    else {
      c2 <- uno_c(tec$OS.time, tec$OS, m2$risk, tau)
      
      ibs2 <- integrated_brier(
        tec$OS.time, tec$OS, m2$surv, m2$times, tau,
        cens_time = trc$OS.time,
        cens_status = trc$OS
      )
      ibs2_floored <- attr(ibs2, "n_floored")
      if (is.null(ibs2_floored)) ibs2_floored <- NA_integer_
      
      s2 <- c(c2, as.numeric(ibs2), .dcal(m2$surv_d, dgrid, tec))
    }
  }
  
  secs2_total <- if (isTRUE(RUN_FUSED_ARM)) .now() - t2 else NA_real_
  
  #offset arm (Volkmann): clinical LP frozen as an offset
  t3 <- .now()
  m3 <- NULL; s3 <- c(NA_real_, NA_real_, NA_real_); ibs3_floored <- NA_integer_
  tm3 <- list(secs = NA_real_)
  if (isTRUE(RUN_OFFSET_ARM)) {
    tm3 <- timed(fit_offset_rna(
      cl$train, cl$test, rn$train, rn$test, trc$OS.time, trc$OS,
      pgrid, foldid = inner_foldid,
      clinical_subset = subset, clinical_train = trc, inner_rna = inner_rna
    ))
    m3 <- .split_surv(tm3$value, grid, dgrid)
    if (!(length(m3$risk) == 1 && is.na(m3$risk))) {
      c3   <- uno_c(tec$OS.time, tec$OS, m3$risk, tau)
      ibs3 <- integrated_brier(
        tec$OS.time, tec$OS, m3$surv, m3$times, tau,
        cens_time = trc$OS.time, cens_status = trc$OS
      )
      ibs3_floored <- attr(ibs3, "n_floored")
      if (is.null(ibs3_floored)) ibs3_floored <- NA_integer_
      s3 <- c(c3, as.numeric(ibs3), .dcal(m3$surv_d, dgrid, tec))
    }
  }
  
  secs3_total <- if (isTRUE(RUN_OFFSET_ARM)) .now() - t3 else NA_real_
  
  #forest arm: clinical RSF + blockForest
  m4 <- NULL; m5 <- NULL
  tm4 <- list(secs = NA_real_); tm5 <- list(secs = NA_real_)
  fseed <- forest_seed(r, f, sid)
  if (isTRUE(RUN_FOREST_ARM)) {
    if (!is_anchor) {  #anchor clinical side = KM (m1)
      tm4 <- timed(fit_rsf(cl$train, trc$OS.time, trc$OS, cl$test, pgrid, fseed,
                           mtry = ncol(cl$train)))
      m4  <- .split_surv(tm4$value, grid, dgrid)
    }
    fxf <- fuse(cl, rn)
    tm5 <- timed(fit_blockforest(fxf$train, trc$OS.time, trc$OS, fxf$test,
                                 fxf$block, pgrid, fseed))
    m5  <- .split_surv(tm5$value, grid, dgrid)
  }
  fcols <- forest_row_cols(m4, m5, s1, is_anchor, isTRUE(RUN_FOREST_ARM),
                           tec, trc, tau, fseed, nz(m4$train_secs, NA_real_),
                           nz(m5$train_secs, NA_real_), tm4$secs, tm5$secs)
  
  #boosting arm: CoxBoost favoring clinical
  m6 <- NULL; tm6 <- list(secs = NA_real_)
  if (isTRUE(RUN_COXBOOST_ARM)) {
    fxb <- fuse(cl, rn)
    tm6 <- timed(fit_coxboost(fxb$train, trc$OS.time, trc$OS, fxb$test,
                              fxb$block, pgrid, inner_foldid))
    m6  <- .split_surv(tm6$value, grid, dgrid)
  }
  bcols <- coxboost_row_cols(m6, s1, isTRUE(RUN_COXBOOST_ARM), tec, trc, tau,
                             nz(m6$train_secs, NA_real_), tm6$secs)
  
  #ridge fused arm (alpha = 0) and RNA-score arm
  m7 <- NULL; tm7 <- list(secs = NA_real_); m8 <- NULL; tm8 <- list(secs = NA_real_)
  if (isTRUE(RUN_RIDGE_ARM)) {
    fxr <- fuse(cl, rn)
    tm7 <- timed(fit_clin_rna(fxr$train, trc$OS.time, trc$OS, fxr$test, fxr$penalty, pgrid,
                              foldid = inner_foldid, lambda_min_ratio = RIDGE_LAMBDA_MIN_RATIO,
                              clinical_subset = subset, clinical_train = trc,
                              inner_rna = inner_rna, alpha = 0))
    m7 <- .split_surv(tm7$value, grid, dgrid)
    if (is.numeric(m7$sel_idx) && length(m7$sel_idx) == 1 && !is.na(m7$sel_idx))
      m7$at_top <- m7$sel_idx == 1L
  }
  if (isTRUE(RUN_SCORE_ARM)) {
    tm8 <- timed(fit_rna_score(cl$train, cl$test, rn$train, rn$test, trc$OS.time, trc$OS,
                               pgrid, inner_foldid, inner_rna))
    m8 <- .split_surv(tm8$value, grid, dgrid)
  }
  rcols <- .simple_arm_cols("ridge", m7, s1, isTRUE(RUN_RIDGE_ARM), tec, trc, tau,
                            nz(m7$train_secs, NA_real_), RIDGE_EXTRA, tm7$secs)
  scols <- .simple_arm_cols("score", m8, s1, isTRUE(RUN_SCORE_ARM), tec, trc, tau,
                            nz(m8$train_secs, NA_real_), SCORE_EXTRA, tm8$secs)
  
  row <- data.frame(
    cancer = code, subset_id = sid, n_vars = length(subset), variables = label,
    is_anchor = is_anchor, repeat_id = r, seed = seed_r, fold = f,
    c_clin = s1[1], ibs_clin = s1[2],
    c_fused = s2[1], ibs_fused = s2[2],
    delta_c = s2[1] - s1[1], delta_ibs = s2[2] - s1[2],
    #D-calibration statistic: lower = better calibrated
    dcal_clin = s1[3], dcal_fused = s2[3], delta_dcal = s2[3] - s1[3],
    dcal_offset = s3[3], delta_dcal_offset = s3[3] - s1[3],
    fused_n_path      = if (is.null(m2$n_path))      NA else m2$n_path,
    fused_sel_idx     = if (is.null(m2$sel_idx))     NA else m2$sel_idx,
    #three separate fused diagnostics: path_short, solver_error, at_boundary
    fused_path_short   = if (is.null(m2$n_path)) NA else m2$n_path < NLAMBDA,
    fused_solver_error = if (is.null(m2$jerr))   NA else m2$jerr != 0,
    fused_at_boundary  = if (is.null(m2$at_boundary)) NA else m2$at_boundary,
    fused_n_genes      = if (is.null(m2$n_genes))     NA else m2$n_genes,
    fused_lambda_fraction    = nz(m2$lambda_fraction,    NA_real_),
    fused_inner_n_path       = nz(m2$inner_n_path,       NA_integer_),
    fused_inner_at_boundary  = nz(m2$inner_at_boundary,  NA),
    fused_inner_solver_error = nz(m2$inner_solver_error, NA),
    #INNER-SPLIT CLINICAL WIDTH
    fused_outer_clin_cols     = nz(m2$outer_clin_cols,     NA_integer_),
    fused_inner_clin_cols_min = nz(m2$inner_clin_cols_min, NA_integer_),
    fused_inner_clin_cols_max = nz(m2$inner_clin_cols_max, NA_integer_),
    fused_inner_clin_mismatch = nz(m2$inner_clin_mismatch, NA),
    fused_inner_vars_dropped  = nz(m2$inner_vars_dropped,  NA_character_),
    #penalized fallback: where it was used (inner1..5, final) and why (error, short, start, beyond, predict)
    fused_fallback     = nz(m2$fallback,     NA_character_),
    fused_fallback_why = nz(m2$fallback_why, NA_character_),
    #CLINICAL-ARM CONVERGENCE, the other half of delta_c
    clin_n_coef       = nz(m1$n_coef,       NA_integer_),
    clin_max_abs_coef = nz(m1$max_abs_coef, NA_real_),
    clin_n_bad_coef   = nz(m1$n_bad_coef,   NA_integer_),
    clin_n_iter       = nz(m1$n_iter,       NA_integer_),
    clin_separated    = nz(m1$separated,    NA),
    clin_warning      = nz(m1$warning,      NA_character_),
    #OFFSET ARM (Volkmann)
    c_offset           = s3[1],
    ibs_offset         = s3[2],
    delta_c_offset     = s3[1] - s1[1],
    delta_ibs_offset   = s3[2] - s1[2],
    offset_n_path            = nz(m3$n_path,            NA_integer_),
    offset_sel_idx           = nz(m3$sel_idx,           NA_integer_),
    offset_path_short        = if (is.null(m3$n_path)) NA else m3$n_path < OFFSET_NLAMBDA,
    offset_solver_error      = if (is.null(m3$jerr))   NA else m3$jerr != 0,
    offset_at_boundary       = nz(m3$at_boundary,       NA),
    offset_n_genes           = nz(m3$n_genes,           NA_integer_),
    offset_lambda_fraction   = nz(m3$lambda_fraction,   NA_real_),
    offset_inner_n_path      = nz(m3$inner_n_path,      NA_integer_),
    offset_inner_at_boundary = nz(m3$inner_at_boundary, NA),
    offset_inner_solver_error= nz(m3$inner_solver_error,NA),
    offset_inner_clin_mismatch = nz(m3$inner_clin_mismatch, NA),
    offset_error             = nz(m3$error,             NA_character_),
    ibs_offset_floored       = ibs3_floored,
    secs_offset              = nz(m3$train_secs, NA_real_),
    fused_error        = nz(m2$error, NA_character_),
    clin_error         = nz(m1$error, NA_character_),
    ibs_clin_floored   = ibs1_floored,
    ibs_fused_floored  = ibs2_floored,
    secs_clin  = nz(m1$train_secs, NA_real_),  #training time only (see header)
    secs_fused = nz(m2$train_secs, NA_real_),
    secs_clin_total   = round(secs1_total, 4),  #training + prediction + scoring
    secs_fused_total  = round(secs2_total, 4),
    secs_offset_total = round(secs3_total, 4),
    secs_total = round(.now() - t_cell, 4),  #the whole cell, every model
    secs_rna_block = rna_secs,  #per-FOLD cost, repeated on each row
    n_train = nrow(trc), n_test = nrow(tec),
    n_events_train = sum(trc$OS),
    n_cols_fused   = ncol(cl$train) + ncol(rn$train),  #not fx: fx is absent when RUN_FUSED_ARM = FALSE
    host = HOST,
    stringsAsFactors = FALSE
  )
  row <- cbind(row, fcols, bcols, rcols, scols)
  list(row   = row,
       genes = if (is.null(m2$genes)) character(0) else m2$genes,
       preds = if (!SAVE_PREDICTIONS) NULL else list(
         key       = sprintf("%s|r%d|f%d|s%d", code, r, f, sid),
         test_time = tec$OS.time, test_status = tec$OS,
         tau = tau, times = grid,
         risk_clin = m1$risk, surv_clin = m1$surv,
         risk_fused = m2$risk, surv_fused = m2$surv,
         risk_offset = if (is.null(m3)) NULL else m3$risk,
         surv_offset = if (is.null(m3) || !isTRUE(SAVE_OFFSET_SURV)) NULL else m3$surv,
         #forest predictions per tree count
         forest = if (!isTRUE(SAVE_FOREST_PREDS) || !isTRUE(RUN_FOREST_ARM)) NULL else
           lapply(setNames(FOREST_TREES, FOREST_TREES), function(nt) {
             k  <- as.character(nt)
             pc <- if (is_anchor) list(risk = m1$risk, surv = m1$surv) else
               if (is.null(m4$by_trees)) NULL else m4$by_trees[[k]]
             pb <- if (is.null(m5$by_trees)) NULL else m5$by_trees[[k]]
             list(risk_rfclin = pc$risk, surv_rfclin = pc$surv,
                  risk_bf = pb$risk, surv_bf = pb$surv)
           })))
}

#compute-time tracking: results/timing/<key>_timing.csv, _timing_repeats.csv and timing_log.csv
hardware_info <- function() {
  sys  <- tryCatch(Sys.info(), error = function(e) c(sysname = NA, release = NA))
  cpu  <- NA_character_; ram_gb <- NA_real_
  if (.Platform$OS.type == "windows") {
    cpu <- Sys.getenv("PROCESSOR_IDENTIFIER", NA_character_)
    ram <- tryCatch(suppressWarnings(system2("powershell",
                                             c("-NoProfile", "-Command",
                                               "(Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory"),
                                             stdout = TRUE, stderr = FALSE)), error = function(e) character(0))
    ram <- suppressWarnings(as.numeric(trimws(ram[nzchar(trimws(ram))][1])))
    if (length(ram) && !is.na(ram)) ram_gb <- round(ram / 1024^3, 1)
    nm <- tryCatch(suppressWarnings(system2("powershell",
                                            c("-NoProfile", "-Command", "(Get-CimInstance Win32_Processor).Name"),
                                            stdout = TRUE, stderr = FALSE)), error = function(e) character(0))
    nm <- trimws(nm[nzchar(trimws(nm))])
    if (length(nm)) cpu <- nm[1]
  } else if (file.exists("/proc/cpuinfo")) {
    l <- grep("^model name", readLines("/proc/cpuinfo"), value = TRUE)
    if (length(l)) cpu <- trimws(sub(".*:", "", l[1]))
    m <- grep("^MemTotal", readLines("/proc/meminfo"), value = TRUE)
    if (length(m)) ram_gb <- round(as.numeric(gsub("\\D", "", m[1])) / 1024^2, 1)
  }
  data.frame(host = HOST, os = paste(sys[["sysname"]], sys[["release"]]),
             cpu = cpu, cores_logical = parallel::detectCores(logical = TRUE),
             cores_physical = parallel::detectCores(logical = FALSE), ram_gb = ram_gb,
             r_version = paste(R.version$major, R.version$minor, sep = "."),
             glmnet = as.character(utils::packageVersion("glmnet")),
             survival = as.character(utils::packageVersion("survival")),
             stringsAsFactors = FALSE)
}

.arm_total <- function(x) if (all(is.na(x))) NA_real_ else round(sum(x, na.rm = TRUE), 1)
.arm_median <- function(x) if (all(is.na(x))) NA_real_ else round(median(x, na.rm = TRUE), 3)
.arm_max <- function(x) if (all(is.na(x))) NA_real_ else round(max(x, na.rm = TRUE), 1)

write_timing <- function(code, key, res, rep_info, wall_secs, n_specs, n_patients, fold_times = NULL) {
  dir <- file.path(DIR_RESULTS, "timing"); dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  g <- function(col) if (col %in% names(res)) res[[col]] else rep(NA_real_, nrow(res))
  rna <- unique(res[, intersect(c("repeat_id", "fold", "secs_rna_block"), names(res)), drop = FALSE])
  row <- data.frame(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), cancer = code, run_key = key,
    n_patients = n_patients, n_specs = n_specs, n_cells = nrow(res),
    n_repeats = length(rep_info$repeat_id), folds = OUTER_FOLDS, workers = N_WORKERS,
    repeats_computed_now = sum(rep_info$source == "computed"),
    repeats_from_checkpoint = sum(rep_info$source == "checkpoint"),
    wall_secs_this_call = round(wall_secs, 1),
    wall_secs_all_repeats = .arm_total(rep_info$wall_secs),
    compute_secs_rna_block = .arm_total(rna$secs_rna_block),
    compute_secs_clin = .arm_total(g("secs_clin")),     median_cell_clin = .arm_median(g("secs_clin")),
    compute_secs_fused = .arm_total(g("secs_fused")),   median_cell_fused = .arm_median(g("secs_fused")),
    compute_secs_offset = .arm_total(g("secs_offset")), median_cell_offset = .arm_median(g("secs_offset")),
    compute_secs_rsf = .arm_total(g("secs_rfclin")),    median_cell_rsf = .arm_median(g("secs_rfclin")),
    compute_secs_blockforest = .arm_total(g("secs_bf")), median_cell_blockforest = .arm_median(g("secs_bf")),
    compute_secs_coxboost = .arm_total(g("secs_cb")), median_cell_coxboost = .arm_median(g("secs_cb")),
    compute_secs_ridge = .arm_total(g("secs_ridge")), compute_secs_score = .arm_total(g("secs_score")),
    #totals: training + prediction + scoring per model
    total_secs_clin = .arm_total(g("secs_clin_total")), total_secs_fused = .arm_total(g("secs_fused_total")),
    total_secs_offset = .arm_total(g("secs_offset_total")), total_secs_rsf = .arm_total(g("secs_rfclin_total")),
    total_secs_blockforest = .arm_total(g("secs_bf_total")), total_secs_coxboost = .arm_total(g("secs_cb_total")),
    total_secs_ridge = .arm_total(g("secs_ridge_total")), total_secs_score = .arm_total(g("secs_score_total")),
    total_secs_all_cells = .arm_total(g("secs_total")),
    #peak memory over the repeats (MB)
    peak_ram_main_mb   = .arm_max(rep_info$ram_main_mb),
    peak_ram_worker_mb = .arm_max(rep_info$ram_worker_max_mb),
    peak_ram_total_mb  = .arm_max(rowSums(cbind(rep_info$ram_main_mb,
                                                rep_info$ram_workers_sum_mb), na.rm = TRUE)),
    stringsAsFactors = FALSE)
  arms <- c("compute_secs_rna_block", "compute_secs_clin", "compute_secs_fused",
            "compute_secs_offset", "compute_secs_rsf", "compute_secs_blockforest",
            "compute_secs_coxboost", "compute_secs_ridge", "compute_secs_score")
  row$compute_secs_all_arms <- round(sum(unlist(row[arms]), na.rm = TRUE), 1)
  row <- cbind(row, hardware_info())
  write.csv(row, file.path(dir, paste0(key, "_timing.csv")), row.names = FALSE)
  #per-repeat breakdown: which machine computed it and what it cost
  per <- lapply(rep_info$repeat_id, function(r) {
    rr <- res[res$repeat_id == r, , drop = FALSE]
    if (!nrow(rr)) return(data.frame(host = NA_character_, compute_secs_cox = NA_real_,
                                     compute_secs_forest = NA_real_))
    gr  <- function(col) if (col %in% names(rr)) rr[[col]] else NA_real_
    rna <- unique(rr[, intersect(c("fold", "secs_rna_block"), names(rr)), drop = FALSE])
    hs  <- unique(stats::na.omit(as.character(gr("host"))))
    data.frame(host = if (length(hs)) paste(hs, collapse = "+") else NA_character_,
               compute_secs_cox = round(sum(c(rna$secs_rna_block, gr("secs_clin"),
                                              gr("secs_fused"), gr("secs_offset"),
                                              gr("secs_cb"), gr("secs_ridge"), gr("secs_score")),
                                            na.rm = TRUE), 1),
               compute_secs_forest = .arm_total(c(gr("secs_rfclin"), gr("secs_bf"))))
  })
  write.csv(cbind(as.data.frame(rep_info), do.call(rbind, per)),
            file.path(dir, paste0(key, "_timing_repeats.csv")), row.names = FALSE)
  #per-fold wall time (repeats from checkpoints written before this change have none)
  if (!is.null(fold_times) && nrow(fold_times)) {
    fold_times$source <- rep_info$source[match(fold_times$repeat_id, rep_info$repeat_id)]
    write.csv(fold_times, file.path(dir, paste0(key, "_fold_timing.csv")), row.names = FALSE)
  }
  logp <- file.path(dir, "timing_log.csv")
  write.table(row, logp, sep = ",", row.names = FALSE, col.names = !file.exists(logp),
              append = file.exists(logp))
  cat(sprintf("Saved %s  (compute %.1f min all arms | wall %.1f min this call)\n",
              file.path(dir, paste0(key, "_timing.csv")),
              row$compute_secs_all_arms / 60, wall_secs / 60))
  invisible(row)
}

#one parallel task: runs subset `si` with the fold objects exported to the worker
.run_cell_on_worker <- function(si)
  with(fold_env, run_one_cell(subset_list[[si]], subset_ids[si], code, r, f, seed_r,
                              trc, tec, rn, tau, grid, dgrid, rna_secs,
                              inner_foldid, inner_rna))

#peak memory of this R process since the last gc(reset = TRUE), in MB.
.peak_mb <- function() {
  g <- gc()
  round(sum(g[, which(colnames(g) == "max used") + 1]), 1)
}

run_one_repeat <- function(clin, expr, status, subset_list, subset_ids,
                           code, n_genes, r, seed_r, cl_workers = NULL,
                           n_repeats = n_repeats_for(code)) {
  folds <- make_outer_folds(status, seed = seed_r)  #shared across subsets this repeat
  cat(sprintf("\n===== repeat %d/%d (seed %d) =====\n", r, n_repeats, seed_r))
  rep_rows <- list()
  gene_rows <- list()
  pred_rows <- list()
  fold_rows <- list()  #one timing row per outer fold, saved with the checkpoint
  fold_row <- function(f, status, n_cells)
    data.frame(repeat_id = r, fold = f, status = status, n_cells = n_cells,
               workers = if (is.null(cl_workers)) 1L else N_WORKERS,
               wall_secs = round(proc.time()[["elapsed"]] - fold_t0, 2),
               host = Sys.info()[["nodename"]], finished = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
               stringsAsFactors = FALSE)
  
  for (f in seq_len(OUTER_FOLDS)) {
    fold_t0 <- proc.time()[["elapsed"]]
    
    tei <- which(folds == f)
    tri <- which(folds != f)
    trc <- clin[tri, ]
    tec <- clin[tei, ]
    tre <- expr[tri, , drop = FALSE]
    tee <- expr[tei, , drop = FALSE]
    
    if (sum(trc$OS) < MIN_TRAIN_EVENTS) {  #low-event guard
      cat(sprintf("  fold %d: SKIPPED (only %d train events)\n", f, sum(trc$OS)))
      fold_rows[[length(fold_rows) + 1]] <- fold_row(f, "skipped_low_events", 0L)
      next
    }
    
    #BUILT ONCE PER FOLD, REUSED BY EVERY SUBSET
    tau  <- make_tau(trc$OS.time, trc$OS)  #Uno C horizon: max TRAIN event time
    grid <- sort(unique(quantile(trc$OS.time[trc$OS == 1], GRID_QUANTILES)))  #IBS: TRAIN q10-q90
    grid <- grid[grid <= tau]
    if (length(grid) < 2) {
      cat(sprintf("  fold %d: SKIPPED (grid too short)\n", f))
      fold_rows[[length(fold_rows) + 1]] <- fold_row(f, "skipped_grid", 0L)
      next
    }
    dgrid <- sort(unique(trc$OS.time[trc$OS == 1]))  #D-calibration: every TRAIN event time
    
    #ONE stratified inner-CV partition for this outer-training set.
    inner_foldid <- make_inner_folds(
      trc$OS, k = INNER_FOLDS, seed = inner_fold_seed(r, f)
    )
    
    trna <- timed(list(
      outer = build_rna(tre, tee, n_genes),
      inner = build_inner_rna(tre, inner_foldid, n_genes)
    ))
    rn <- trna$value$outer
    inner_rna <- trna$value$inner
    cat(sprintf("  fold %d: outer + inner RNA blocks built in %.2fs (%d outer genes), %d subsets%s ...\n",
                f, trna$secs, ncol(rn$train), length(subset_list),
                if (!is.null(cl_workers)) sprintf(" on %d workers", N_WORKERS) else ""))
    
    if (is.null(cl_workers)) {
      out <- lapply(seq_along(subset_list), function(si)
        run_one_cell(subset_list[[si]], subset_ids[si], code, r, f, seed_r,
                     trc, tec, rn, tau, grid, dgrid, trna$secs, inner_foldid, inner_rna))
    } else {
      #send this fold's objects to every worker once, then only the subset number per task
      fold_env <- list(subset_list = subset_list, subset_ids = subset_ids, code = code,
                       r = r, f = f, seed_r = seed_r, trc = trc, tec = tec, rn = rn,
                       tau = tau, grid = grid, dgrid = dgrid, rna_secs = trna$secs,
                       inner_foldid = inner_foldid, inner_rna = inner_rna)
      parallel::clusterExport(cl_workers, "fold_env", envir = environment())
      out <- parallel::parLapplyLB(cl_workers, seq_along(subset_list), .run_cell_on_worker)
    }
    
    for (o in out) {
      if (is.null(o)) next
      rep_rows[[length(rep_rows) + 1]] <- o$row
      if (SAVE_PREDICTIONS && !is.null(o$preds))
        pred_rows[[length(pred_rows) + 1]] <- o$preds
      if (SAVE_GENE_NAMES && length(o$genes))
        gene_rows[[length(gene_rows) + 1]] <- data.frame(
          cancer = code, subset_id = o$row$subset_id, n_vars = o$row$n_vars,
          repeat_id = r, fold = f, gene = o$genes, stringsAsFactors = FALSE)
    }
    fold_rows[[length(fold_rows) + 1]] <- fold_row(f, "computed", sum(!vapply(out, is.null, logical(1))))
    cat(sprintf("  repeat %d fold %d done in %.1fs\n", r, f, fold_rows[[length(fold_rows)]]$wall_secs))
  }
  
  list(results = do.call(rbind, rep_rows),
       genes   = if (length(gene_rows)) do.call(rbind, gene_rows) else NULL,
       preds   = if (length(pred_rows)) pred_rows else NULL,
       fold_times = do.call(rbind, fold_rows))
}

#run one cancer end to end with checkpoints.
run_cancer <- function(dat, code, n_genes = TOP_N_GENES,
                       subsets = NULL, key = code, n_repeats = n_repeats_for(code)) {
  wall0 <- proc.time()[["elapsed"]]
  
  ok   <- !is.na(dat$clinical$OS) & !is.na(dat$clinical$OS.time) & dat$clinical$OS.time > 0
  clin <- dat$clinical[ok, ]; expr <- dat$expr[ok, , drop = FALSE]
  status <- clin$OS
  
  if (is.null(subsets)) {
    real_subsets <- all_subsets(names(CANCER_VARS[[code]]))
    subset_list  <- c(list(character(0)), real_subsets)  #anchor first (subset_id 0)
    subset_ids   <- c(0L, seq_along(real_subsets))
    spec_desc    <- NULL
  } else {
    if (is.null(names(subsets)) || anyDuplicated(names(subsets)))
      stop("`subsets` must be a list with unique names")
    miss <- setdiff(unique(unlist(subsets)), names(dat$clinical))
    if (length(miss))
      stop("variables not in the modelling data for ", key, ": ", paste(miss, collapse = ", "))
    subset_list <- unname(subsets)
    subset_ids  <- seq_along(subsets)
    spec_desc   <- paste(sprintf("%s:%s", names(subsets),
                                 vapply(subsets, paste, character(1), collapse = "+")),
                         collapse = ";")
    dir.create(DIR_RESULTS, showWarnings = FALSE, recursive = TRUE)
    write.csv(data.frame(subset_id = subset_ids, spec = names(subsets),
                         variables = vapply(subsets, paste, character(1), collapse = "+"),
                         n_vars = lengths(subsets)),
              file.path(DIR_RESULTS, paste0(key, "_specs.csv")), row.names = FALSE)
  }
  N_REPEATS    <- as.integer(n_repeats)  #default: the compute rule, 5 for large lattices, 10 otherwise
  RARE_MISSING_VARS <<- rare_missing_vars(clin, unique(unlist(subset_list)), MIN_LEVEL_N)
  if (length(RARE_MISSING_VARS))
    cat(sprintf("rare missing groups (under %d patients, imputed, no missing column): %s\n",
                MIN_LEVEL_N, paste(RARE_MISSING_VARS, collapse = ", ")))
  seeds        <- SEED + seq_len(N_REPEATS)
  
  #stale-checkpoint guard
  fp   <- settings_fingerprint(code, n_genes, spec_desc, n_repeats = N_REPEATS)
  fppath <- file.path(DIR_CHECKPOINTS, paste0(key, "_settings.txt"))
  if (file.exists(fppath)) {
    old <- paste(readLines(fppath), collapse = "")
    if (!identical(old, fp))
      stop("Checkpoints for ", key, " were written under different settings.\n",
           "  was: ", old, "\n  now: ", fp, "\n",
           "  Delete the checkpoints for this cancer and re-run.")
  } else {
    old_checkpoints <- vapply(seq_len(N_REPEATS), function(r)
      file.exists(path_checkpoint(key, r)), logical(1))
    if (any(old_checkpoints))
      stop("Existing checkpoints have no settings fingerprint. Re-run with force = TRUE.")
    dir.create(dirname(fppath), recursive = TRUE, showWarnings = FALSE)
    writeLines(fp, fppath)
  }
  
  cat(sprintf("=== %s: %d patients, %d %s x %d folds x %d repeats ===\n",
              key, nrow(clin), length(subset_list),
              if (is.null(subsets)) "subsets (incl anchor)" else "specifications",
              OUTER_FOLDS, N_REPEATS))
  
  #worker pool, created once and reused across every fold
  cl_workers <- NULL
  if (N_WORKERS > 1) {
    cl_workers <- parallel::makeCluster(N_WORKERS)
    #Cluster teardown is housekeeping, not part of model estimation.
    on.exit({
      tryCatch(
        parallel::stopCluster(cl_workers),
        error = function(e)
          warning("Worker cluster cleanup failed after model fitting: ",
                  conditionMessage(e), call. = FALSE)
      )
    }, add = TRUE)
    parallel::clusterEvalQ(cl_workers, {
      source("paths.R"); source_library(); NULL
    })
    #export every setting that can change a result (add new ones here)
    SETTINGS_CONSTANTS <- c("NLAMBDA", "LAMBDA_MIN_RATIO", "COX_TIES",
                            "COEF_SEPARATION", "INNER_FOLDS",
                            "TOP_N_GENES", "RNA_MIN_COMPLETENESS", "G_FLOOR",
                            "GRID_QUANTILES", "MIN_TRAIN_EVENTS",
                            "RUN_OFFSET_ARM", "SAVE_OFFSET_SURV", "RUN_FUSED_ARM",
                            "OFFSET_NLAMBDA", "OFFSET_LAMBDA_MIN_RATIO",
                            "OUTER_FOLDS", "N_REPEATS", "SEED",
                            "N_REPEATS_LARGE", "N_REPEATS_SMALL", "REPEATS_SUBSET_THRESHOLD",
                            "SAVE_PREDICTIONS", "SAVE_GENE_NAMES",
                            "FOREST_TREES", "BF_NUM_TREES_PRE", "BF_NSETS", "BF_WEIGHT_GRID",
                            "RF_SPLITRULE", "RF_THREADS", "RUN_FOREST_ARM", "BF_FIX_CLINICAL_WEIGHT",
                            "RUN_COXBOOST_ARM", "COXBOOST_MAXSTEPNO",
                            "RUN_RIDGE_ARM", "RUN_SCORE_ARM", "RIDGE_LAMBDA_MIN_RATIO", "D_CAL_BINS",
                            "SCORE_ALPHA", "SCORE_LAMBDA_MIN_RATIO", "RARE_MISSING_VARS",
                            "FOREST_SEED_OFFSET", "SAVE_FOREST_PREDS",
                            "CLINICAL_MISSING", "MISSING_LEVEL")
    missing_const <- SETTINGS_CONSTANTS[!vapply(SETTINGS_CONSTANTS, exists, logical(1))]
    if (length(missing_const))
      stop("Settings constant(s) not found before clusterExport: ",
           paste(missing_const, collapse = ", "))
    parallel::clusterExport(cl_workers,
                            c("run_one_cell", "timed", "nz", ".simple_arm_cols", "RIDGE_EXTRA", "SCORE_EXTRA",
                              "ANCHOR_LABEL", "HOST", SETTINGS_CONSTANTS),
                            envir = environment())
  }
  
  results <- vector("list", N_REPEATS)
  rep_info <- data.frame(repeat_id = seq_len(N_REPEATS), source = NA_character_,
                         wall_secs = NA_real_, ram_main_mb = NA_real_,
                         ram_worker_max_mb = NA_real_, ram_workers_sum_mb = NA_real_,
                         stringsAsFactors = FALSE)
  ram_cols <- c("ram_main_mb", "ram_worker_max_mb", "ram_workers_sum_mb")
  for (r in seq_len(N_REPEATS)) {
    ckpt <- path_checkpoint(key, r)
    if (file.exists(ckpt)) {
      cat(sprintf("===== repeat %d/%d: checkpoint found, loading =====\n", r, N_REPEATS))
      results[[r]] <- readRDS(ckpt)
      rep_info$source[r] <- "checkpoint"
      if (is.list(results[[r]]) && !is.data.frame(results[[r]]) && !is.null(results[[r]]$wall_secs))
        rep_info$wall_secs[r] <- results[[r]]$wall_secs
      for (k in ram_cols) if (!is.null(results[[r]][[k]])) rep_info[[k]][r] <- results[[r]][[k]]
      next
    }
    #peak memory is measured per repeat: reset the counters, run, read them
    invisible(gc(reset = TRUE))
    if (!is.null(cl_workers)) invisible(parallel::clusterEvalQ(cl_workers, gc(reset = TRUE)))
    rt0 <- proc.time()[["elapsed"]]
    results[[r]] <- run_one_repeat(clin, expr, status, subset_list, subset_ids,
                                   code, n_genes, r, seeds[r], cl_workers,
                                   n_repeats = N_REPEATS)
    results[[r]]$wall_secs <- round(proc.time()[["elapsed"]] - rt0, 2)
    results[[r]]$ram_main_mb <- .peak_mb()
    wk <- if (is.null(cl_workers)) numeric(0) else
      unlist(parallel::clusterEvalQ(cl_workers, .peak_mb()))
    results[[r]]$ram_worker_max_mb  <- if (length(wk)) max(wk) else NA_real_
    results[[r]]$ram_workers_sum_mb <- if (length(wk)) sum(wk) else NA_real_
    rep_info$source[r] <- "computed"; rep_info$wall_secs[r] <- results[[r]]$wall_secs
    for (k in ram_cols) rep_info[[k]][r] <- results[[r]][[k]]
    saveRDS(results[[r]], ckpt)
    cat(sprintf("saved checkpoint: %s\n", ckpt))
  }
  
  #tolerate checkpoints written before the list-return change
  results <- lapply(results, function(x)
    if (is.data.frame(x)) list(results = x, genes = NULL) else x)
  
  res  <- do.call(rbind, lapply(results, `[[`, "results"))
  gsel <- do.call(rbind, lapply(results, `[[`, "genes"))
  prd  <- unlist(lapply(results, `[[`, "preds"), recursive = FALSE)
  ftm  <- do.call(rbind, lapply(results, `[[`, "fold_times"))  #NULL for old checkpoints
  
  write.csv(res, path_results(key), row.names = FALSE)
  cat("\nSaved ", path_results(key), "  (", nrow(res), " rows )\n", sep = "")
  
  if (!is.null(gsel)) {
    #derive from path_results so it follows any redirect (e.g. the gene sweep)
    gpath <- sub("_results\\.csv$", "_genes.csv.gz", path_results(key))
    gz <- gzfile(gpath, "w"); write.csv(gsel, gz, row.names = FALSE); close(gz)
    cat(sprintf("Saved %s  (%d gene selections, %d distinct genes)\n",
                gpath, nrow(gsel), length(unique(gsel$gene))))
  }
  
  if (!is.null(prd)) {
    ppath <- sub("_results\\.csv$", "_preds.rds", path_results(key))
    saveRDS(prd, ppath, compress = "xz")
    cat(sprintf("Saved %s  (%d cells; metrics can now be recomputed without refitting)\n",
                ppath, length(prd)))
  }
  
  fit_secs  <- sum(res$secs_total, na.rm = TRUE)
  wall_secs <- proc.time()[["elapsed"]] - wall0
  write_manifest(code, seeds, length(subset_list), nrow(res), fit_secs, n_genes, key)
  tryCatch(write_timing(code, key, res, rep_info, wall_secs, length(subset_list), nrow(clin),
                        fold_times = ftm),
           error = function(e) warning("timing file not written: ", conditionMessage(e), call. = FALSE))
  
  cat(sprintf("fused diagnostics: path_short %d | solver_error %d | at_boundary %d | errored %d (of %d)\n",
              sum(res$fused_path_short,   na.rm = TRUE),
              sum(res$fused_solver_error, na.rm = TRUE),
              sum(res$fused_at_boundary,  na.rm = TRUE),
              sum(!is.na(res$fused_error)), nrow(res)))
  cat(sprintf("clin  diagnostics: separated %d | non-finite coef %d | errored %d | max |beta| %.1f\n",
              sum(res$clin_separated, na.rm = TRUE),
              sum(res$clin_n_bad_coef > 0, na.rm = TRUE),
              sum(!is.na(res$clin_error)),
              suppressWarnings(max(res$clin_max_abs_coef, na.rm = TRUE))))
  n_mis <- sum(res$fused_inner_clin_mismatch, na.rm = TRUE)
  if (n_mis > 0) {
    dropped <- res$fused_inner_vars_dropped[!is.na(res$fused_inner_vars_dropped)]
    tb <- sort(table(dropped), decreasing = TRUE)
    cat(sprintf("inner-split clinical width mismatch: %d cells (%.1f%%)\n",
                n_mis, 100 * n_mis / nrow(res)))
    cat(sprintf("  variables lost in >=1 inner split: %s\n",
                paste(sprintf("%s (%d)", names(tb)[seq_len(min(5, length(tb)))],
                              as.integer(tb)[seq_len(min(5, length(tb)))]),
                      collapse = ", ")))
  }
  if (any(!is.na(res$fused_error))) {
    tb <- sort(table(substr(res$fused_error[!is.na(res$fused_error)], 1, 70)),
               decreasing = TRUE)
    cat("  most common fused error(s):\n")
    for (i in seq_len(min(3, length(tb))))
      cat(sprintf("    [%d x] %s\n", tb[i], names(tb)[i]))
  }
  nfl <- sum(res$ibs_fused_floored, na.rm = TRUE) + sum(res$ibs_clin_floored, na.rm = TRUE)
  if (nfl > 0)
    cat(sprintf("  IBS censoring floor bound %d times -- see metrics.R G_FLOOR\n", nfl))
  if (isTRUE(RUN_FOREST_ARM)) {
    cat(sprintf("forest diagnostics: clinical RSF errored %d | blockForest errored %d (of %d)\n",
                sum(!is.na(res$rfclin_error)), sum(!is.na(res$bf_error)), nrow(res)))
    cat(sprintf("  blockForest block weights: median clinical %.3f | median RNA %.3f\n",
                suppressWarnings(median(res$bf_w_clin, na.rm = TRUE)),
                suppressWarnings(median(res$bf_w_rna,  na.rm = TRUE))))
    cat(sprintf("  forest timing: median clinical RSF %.2fs | median blockForest %.2fs\n",
                suppressWarnings(median(res$secs_rfclin, na.rm = TRUE)),
                suppressWarnings(median(res$secs_bf,     na.rm = TRUE))))
    if (any(!is.na(res$bf_error))) {
      tb <- sort(table(substr(res$bf_error[!is.na(res$bf_error)], 1, 70)), decreasing = TRUE)
      cat("  most common blockForest error(s):\n")
      for (i in seq_len(min(3, length(tb))))
        cat(sprintf("    [%d x] %s\n", tb[i], names(tb)[i]))
    }
  }
  if (isTRUE(RUN_COXBOOST_ARM)) {
    cat(sprintf("coxboost diagnostics: errored %d | at max steps %d | median steps %s | median genes %s | zero-gene fits %.1f%% (of %d)\n",
                sum(!is.na(res$cb_error)), sum(res$cb_at_boundary, na.rm = TRUE),
                suppressWarnings(median(res$cb_stepno, na.rm = TRUE)),
                suppressWarnings(median(res$cb_n_genes, na.rm = TRUE)),
                100 * mean(res$cb_n_genes == 0, na.rm = TRUE), nrow(res)))
    if (any(!is.na(res$cb_error))) {
      tb <- sort(table(substr(res$cb_error[!is.na(res$cb_error)], 1, 70)), decreasing = TRUE)
      for (i in seq_len(min(3, length(tb)))) cat(sprintf("    [%d x] %s\n", tb[i], names(tb)[i]))
    }
  }
  cat(sprintf("genes selected: median %s, zero-gene fits %.1f%%\n",
              median(res$fused_n_genes, na.rm = TRUE),
              100 * mean(res$fused_n_genes == 0, na.rm = TRUE)))
  cat(sprintf("timing: %.1f min of cell compute, all models (wall %.1f min, %dx workers) | median cell %.2fs\n",
              fit_secs / 60, wall_secs / 60, N_WORKERS,
              median(res$secs_total, na.rm = TRUE)))
  res
}

#NOTE: delete a cancer's checkpoints and _settings.txt before re-running with changed settings