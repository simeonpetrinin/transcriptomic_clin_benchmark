#Summarise
#Averages one run's raw CV results into one row per clinical subset for every model
#run (Uno's C, IBS, D-calibration); the only place averaging happens. Also writes a
#timing table and a check figure. Run: summarise_cancer("LUAD")
#Rescue: fused cells refitted by scripts/run_rescue.R replace the failed glmnet cells,
#in the main lattice and in the robustness run (same cells). Rescued cells get no D-calibration
#(not recomputed). summarise_cancer(key, rescue = FALSE) writes the old drop-the-cells version
#with a "_norescue" suffix, as the sensitivity analysis.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(patchwork)
})
#master.R's summarise() masks dplyr's, so dplyr verbs are written dplyr::verb()

#a missing flag is not a passing flag
flag     <- function(x) if (is.null(x)) FALSE else (!is.na(x) & as.logical(x))
has_text <- function(x) if (is.null(x)) FALSE else (!is.na(x) & nzchar(as.character(x)))

#bad-cell rules for the fused arm
FUSED_RULES <- list(
  spec   = function(d) flag(d$fused_path_short) | flag(d$fused_solver_error),
  solver = function(d) flag(d$fused_solver_error) | flag(d$fused_at_boundary),
  both   = function(d) flag(d$fused_path_short) | flag(d$fused_solver_error) |
    flag(d$clin_separated)
)

#the arms: column names and the rule that marks a failed fit
arm_table <- function(res, fused_rule) {
  arms <- list(
    fused    = list(c = "c_fused",  c_ref = "c_clin", ibs = "ibs_fused",  ibs_ref = "ibs_clin",
                    dc = "delta_c", dibs = "delta_ibs",
                    dcal = "dcal_fused", dcal_ref = "dcal_clin", ddcal = "delta_dcal",
                    bad = FUSED_RULES[[fused_rule]]),
    offset   = list(c = "c_offset", c_ref = "c_clin", ibs = "ibs_offset", ibs_ref = "ibs_clin",
                    dc = "delta_c_offset", dibs = "delta_ibs_offset",
                    dcal = "dcal_offset", dcal_ref = "dcal_clin", ddcal = "delta_dcal_offset",
                    bad = function(d) flag(d$offset_path_short) | flag(d$offset_solver_error)),
    ridge    = list(c = "c_ridge",  c_ref = "c_clin", ibs = "ibs_ridge",  ibs_ref = "ibs_clin",
                    dc = "delta_c_ridge", dibs = "delta_ibs_ridge",
                    dcal = "dcal_ridge", dcal_ref = "dcal_clin", ddcal = "delta_dcal_ridge",
                    bad = function(d) has_text(d$ridge_error)),
    coxboost = list(c = "c_cb",     c_ref = "c_clin", ibs = "ibs_cb",     ibs_ref = "ibs_clin",
                    dc = "delta_c_cb", dibs = "delta_ibs_cb",
                    dcal = "dcal_cb", dcal_ref = "dcal_clin", ddcal = "delta_dcal_cb",
                    bad = function(d) has_text(d$cb_error)),
    score    = list(c = "c_score",  c_ref = "c_clin", ibs = "ibs_score",  ibs_ref = "ibs_clin",
                    dc = "delta_c_score", dibs = "delta_ibs_score",
                    dcal = "dcal_score", dcal_ref = "dcal_clin", ddcal = "delta_dcal_score",
                    bad = function(d) has_text(d$score_error)))
  #one blockForest arm per tree count found in the file
  trees <- sort(as.integer(sub("^c_bf_", "", grep("^c_bf_[0-9]+$", names(res), value = TRUE))))
  for (nt in trees)
    arms[[paste0("blockforest_", nt)]] <- list(
      c = paste0("c_bf_", nt), c_ref = paste0("c_rfclin_", nt),
      ibs = paste0("ibs_bf_", nt), ibs_ref = paste0("ibs_rfclin_", nt),
      dc = paste0("delta_c_bf_", nt), dibs = paste0("delta_ibs_bf_", nt),
      dcal = paste0("dcal_bf_", nt), dcal_ref = paste0("dcal_rfclin_", nt),
      ddcal = paste0("delta_dcal_bf_", nt),
      bad = function(d) has_text(d$bf_error) | has_text(d$rfclin_error))
  #keep the arms that were run: columns present and at least one delta_c
  ran <- vapply(arms, function(a) a$dc %in% names(res) && any(!is.na(res[[a$dc]])), logical(1))
  arms[ran]
}

#file paths: paths.R helpers where they exist
out_path <- function(helper, key, suffix, dir = DIR_RESULTS) {
  if (exists(helper, mode = "function")) get(helper)(key) else
    file.path(dir, paste0(key, suffix))
}

#keys the rescue applies to: the main lattice and the robustness run (the gene sweep ran with the fallback)
rescue_code <- function(key) {
  code <- sub("_robust$", "", key)
  if (key %in% c(code, paste0(code, "_robust")) && !grepl("_", code)) code else NA_character_
}

#overlay the rescued fused cells, matched on variables, repeat and fold
apply_rescue <- function(res, key) {
  res$rescued <- FALSE
  code <- rescue_code(key)
  f <- if (is.na(code)) "" else file.path(DIR_RESULTS, paste0(code, "_rescue.csv"))
  if (!nzchar(f) || !file.exists(f)) {
    cat(sprintf("[rescue] %s: no rescue file, nothing overlaid\n", key)); return(res)
  }
  r <- read.csv(f, stringsAsFactors = FALSE)
  id_res <- paste(res$variables, res$repeat_id, res$fold, sep = "|")
  id_r   <- paste(r$variables, r$repeat_id, r$fold, sep = "|")
  if (anyDuplicated(id_r)) stop("duplicate cells in ", basename(f))
  i <- match(id_r, id_res)
  #main lattice: every rescue cell must be there; robustness run: only its own subsets are
  if (anyNA(i) && key == code) stop(sum(is.na(i)), " rescue cells not found in ", key)
  r <- r[!is.na(i), ]; i <- i[!is.na(i)]
  #same split, same clinical model: the clinical C must agree exactly
  if (max(abs(res$c_clin[i] - r$c_clin), na.rm = TRUE) > 1e-9) stop("clinical C differs between ", key, " and the rescue")
  ok <- is.na(r$error) | !nzchar(as.character(r$error))  #cells still failing stay bad
  i <- i[ok]; r <- r[ok, ]
  res$c_fused[i] <- r$c_fused; res$ibs_fused[i] <- r$ibs_fused
  res$delta_c[i] <- r$delta_c; res$delta_ibs[i] <- r$delta_ibs
  if ("fused_n_genes" %in% names(res)) res$fused_n_genes[i] <- r$n_genes
  for (col in intersect(c("dcal_fused", "delta_dcal"), names(res))) res[[col]][i] <- NA_real_  #not recomputed
  for (col in intersect(c("fused_path_short", "fused_solver_error"), names(res))) res[[col]][i] <- FALSE
  res$rescued[i] <- TRUE
  cat(sprintf("[rescue] %s: %d cells overlaid, %d still failing (left out)\n", key, sum(ok), sum(!ok)))
  res
}

#"_norescue" inserted before the extension, for the sensitivity files
norescue_path <- function(f) sub("(\\.[a-z]+)$", "_norescue\\1", f)

#average one arm to one row per subset (folds within a repeat, then repeats)
average_arm <- function(res, a, drop) {
  num <- function(col) if (col %in% names(res)) as.numeric(res[[col]]) else NA_real_
  d <- data.frame(subset_id = res$subset_id, n_vars = res$n_vars,
                  variables = res$variables, is_anchor = as.logical(res$is_anchor),
                  repeat_id = res$repeat_id,
                  c_model = num(a$c), c_ref = num(a$c_ref), delta_c = num(a$dc),
                  ibs_model = num(a$ibs), ibs_ref = num(a$ibs_ref), delta_ibs = num(a$dibs),
                  dcal_model = num(a$dcal), dcal_ref = num(a$dcal_ref), delta_dcal = num(a$ddcal),
                  stringsAsFactors = FALSE)
  d[drop, c("c_model", "c_ref", "delta_c", "ibs_model", "ibs_ref", "delta_ibs",
            "dcal_model", "dcal_ref", "delta_dcal")] <- NA_real_
  d$dcal_ref[is.na(d$dcal_model)] <- NA_real_  #comparator D-calibration on the same cells as the model's (rescued cells have none)
  nan_to_na <- function(x) ifelse(is.nan(x), NA_real_, x)
  #1. folds -> one value per repeat
  per_rep <- d %>%
    dplyr::group_by(subset_id, n_vars, variables, is_anchor, repeat_id) %>%
    dplyr::summarise(
      c_model = nan_to_na(mean(c_model, na.rm = TRUE)),
      c_ref = nan_to_na(mean(c_ref, na.rm = TRUE)),
      delta_c = nan_to_na(mean(delta_c, na.rm = TRUE)),
      ibs_model = nan_to_na(mean(ibs_model, na.rm = TRUE)),
      ibs_ref = nan_to_na(mean(ibs_ref, na.rm = TRUE)),
      delta_ibs = nan_to_na(mean(delta_ibs, na.rm = TRUE)),
      dcal_model = nan_to_na(mean(dcal_model, na.rm = TRUE)),
      dcal_ref = nan_to_na(mean(dcal_ref, na.rm = TRUE)),
      delta_dcal = nan_to_na(mean(delta_dcal, na.rm = TRUE)),
      .groups = "drop")
  #2-3. repeats -> estimate (mean) and spread (SD); the per-repeat table rides along
  out <- per_rep %>%
    dplyr::group_by(subset_id, n_vars, variables, is_anchor) %>%
    dplyr::summarise(
      n_repeats      = dplyr::n(),
      n_repeats_used = sum(!is.na(delta_c)),
      delta_c_sd     = sd(delta_c, na.rm = TRUE),
      delta_ibs_sd   = sd(delta_ibs, na.rm = TRUE),
      delta_dcal_sd  = sd(delta_dcal, na.rm = TRUE),
      c_model   = nan_to_na(mean(c_model, na.rm = TRUE)),
      c_ref     = nan_to_na(mean(c_ref, na.rm = TRUE)),
      delta_c   = nan_to_na(mean(delta_c, na.rm = TRUE)),
      ibs_model = nan_to_na(mean(ibs_model, na.rm = TRUE)),
      ibs_ref   = nan_to_na(mean(ibs_ref, na.rm = TRUE)),
      delta_ibs = nan_to_na(mean(delta_ibs, na.rm = TRUE)),
      dcal_model = nan_to_na(mean(dcal_model, na.rm = TRUE)),
      dcal_ref   = nan_to_na(mean(dcal_ref, na.rm = TRUE)),
      delta_dcal = nan_to_na(mean(delta_dcal, na.rm = TRUE)),
      .groups = "drop")
  attr(out, "per_repeat") <- per_rep
  out
}

summarise_cancer <- function(key = "LUAD", fused_rule = "spec", write_xlsx = FALSE, rescue = TRUE) {
  if (!fused_rule %in% names(FUSED_RULES))
    stop("fused_rule must be one of: ", paste(names(FUSED_RULES), collapse = ", "))
  
  #1. read the raw results, then overlay the rescue (or not, for the sensitivity files)
  res <- read.csv(out_path("path_results", key, "_results.csv"), stringsAsFactors = FALSE)
  res$is_anchor <- as.logical(res$is_anchor)
  if (rescue) res <- apply_rescue(res, key) else res$rescued <- FALSE
  sfx <- if (rescue) identity else norescue_path  #output names: unchanged, or with "_norescue"
  n_rep   <- length(unique(res$repeat_id))
  n_folds <- length(unique(res$fold))
  arms <- arm_table(res, fused_rule)
  if (!length(arms)) stop("no arm has any results in ", key)
  cat(sprintf("=== %s === %d cells | %d repeats x %d folds | arms: %s\n", key, nrow(res),
              n_rep, n_folds, paste(names(arms), collapse = ", ")))
  if ("fused" %in% names(arms)) {
    alt <- vapply(FUSED_RULES, function(f) sum(f(res) | is.na(res$delta_c)), integer(1))
    cat(sprintf("fused bad cells under each rule: %s (using '%s')\n",
                paste(sprintf("%s %d", names(alt), alt), collapse = " | "), fused_rule))
  }
  
  #2. average every arm, with and without the bad-cell rule
  out <- list(); reps <- list()
  for (arm in names(arms)) {
    a    <- arms[[arm]]
    bad  <- a$bad(res) | is.na(res[[a$dc]])
    clean <- average_arm(res, a, drop = bad)
    all   <- average_arm(res, a, drop = rep(FALSE, nrow(res)))
    resc <- if (arm == "fused") res$rescued & !bad else rep(FALSE, nrow(res))  #only the fused arm was rescued
    cells <- data.frame(subset_id = res$subset_id, bad = bad, resc = resc) %>%
      dplyr::group_by(subset_id) %>%
      dplyr::summarise(n_cells = dplyr::n(), n_bad_cells = sum(bad), n_rescued_cells = sum(resc),
                       .groups = "drop")
    s <- clean %>%
      dplyr::left_join(dplyr::select(all, subset_id, delta_c_all = delta_c,
                                     delta_ibs_all = delta_ibs, delta_dcal_all = delta_dcal),
                       by = "subset_id") %>%
      dplyr::left_join(cells, by = "subset_id")
    s$key <- key; s$arm <- arm; s$primary <- arm == "fused"
    out[[arm]] <- s
    pr <- attr(clean, "per_repeat"); pr$key <- key; pr$arm <- arm
    reps[[arm]] <- pr
    
    #report: how many cells were dropped, and where
    cat(sprintf("\n[%s] dropped %d of %d cells (%.1f%%)%s\n", arm, sum(bad), length(bad),
                100 * mean(bad), if (any(resc)) sprintf(", %d rescued cells kept", sum(resc)) else ""))
    if (any(bad)) {
      gr <- tapply(bad, res$n_vars, mean)
      cat("  dropped share by n_vars: ",
          paste(sprintf("%s:%.0f%%", names(gr), 100 * as.numeric(gr)), collapse = " "), "\n", sep = "")
    }
    lat <- s[!s$is_anchor, ]; anc <- s[s$is_anchor, ]
    if (nrow(anc))
      cat(sprintf("  anchor: C = %.3f vs Kaplan-Meier %.3f (delta %+.3f, SD %.3f)\n",
                  anc$c_model, anc$c_ref, anc$delta_c, anc$delta_c_sd))
    if (nrow(lat)) {
      cat(sprintf("  subsets %d | comparator C %.3f-%.3f | mean delta_c %+.4f | mean delta_ibs %+.4f\n",
                  nrow(lat), min(lat$c_ref, na.rm = TRUE), max(lat$c_ref, na.rm = TRUE),
                  mean(lat$delta_c, na.rm = TRUE), mean(lat$delta_ibs, na.rm = TRUE)))
      cat(sprintf("  RNA raises C in %d of %d subsets; lowers IBS in %d of %d; lowers D-calibration in %d of %d (mean delta %+.3f)\n",
                  sum(lat$delta_c > 0, na.rm = TRUE), sum(!is.na(lat$delta_c)),
                  sum(lat$delta_ibs < 0, na.rm = TRUE), sum(!is.na(lat$delta_ibs)),
                  sum(lat$delta_dcal < 0, na.rm = TRUE), sum(!is.na(lat$delta_dcal)),
                  mean(lat$delta_dcal, na.rm = TRUE)))
      n_empty <- sum(is.na(lat$delta_c))
      if (n_empty) cat(sprintf("  WARNING: %d subsets have no usable estimate\n", n_empty))
      thin <- lat$n_repeats_used < lat$n_repeats
      if (any(thin, na.rm = TRUE))
        cat(sprintf("  %d subsets rest on fewer than %d repeats (min %d)\n",
                    sum(thin, na.rm = TRUE), n_rep, min(lat$n_repeats_used)))
      if (sum(!is.na(lat$delta_c)) > 2)
        cat(sprintf("  corr(delta_c, comparator C) = %+.3f | corr(delta_c, n_vars) = %+.3f (descriptive: subsets share patients and variables)\n",
                    cor(lat$c_ref, lat$delta_c, use = "complete.obs"),
                    cor(lat$n_vars, lat$delta_c, use = "complete.obs")))
    }
  }
  #delta = model - comparator: up = RNA helps for C
  summ <- dplyr::bind_rows(out) %>%
    dplyr::select(key, arm, primary, subset_id, n_vars, variables, is_anchor,
                  n_cells, n_bad_cells, n_rescued_cells, n_repeats, n_repeats_used,
                  c_ref, c_model, delta_c, delta_c_sd, delta_c_all,
                  ibs_ref, ibs_model, delta_ibs, delta_ibs_sd, delta_ibs_all,
                  dcal_ref, dcal_model, delta_dcal, delta_dcal_sd, delta_dcal_all) %>%
    dplyr::arrange(arm, subset_id)
  
  #save to .csv
  f_summ <- sfx(out_path("path_summary", key, "_subset_summary.csv"))
  write.csv(summ, f_summ, row.names = FALSE)
  #per-repeat values (bad cells dropped), for analyses that split the repeats
  f_rep <- sfx(out_path("path_per_repeat", key, "_per_repeat.csv"))
  write.csv(dplyr::bind_rows(reps) %>% dplyr::arrange(arm, subset_id, repeat_id),
            f_rep, row.names = FALSE)
  
  #3. timing: median seconds per fit, by number of clinical variables
  timing <- NULL
  secs <- setdiff(grep("^secs_", names(res), value = TRUE), "secs_rna_block")
  secs <- secs[vapply(secs, function(s) any(!is.na(res[[s]]) & res[[s]] > 0), logical(1))]
  if (length(secs)) {
    timing <- res %>%
      dplyr::group_by(n_vars) %>%
      dplyr::summarise(cells = dplyr::n(),
                       dplyr::across(dplyr::all_of(secs), \(x) median(x, na.rm = TRUE)),
                       .groups = "drop")
    names(timing)[-(1:2)] <- paste0("median_", names(timing)[-(1:2)])
    f_time <- sfx(out_path("path_timing", key, "_timing_summary.csv"))
    write.csv(timing, f_time, row.names = FALSE)
    cat(sprintf("\ntiming: %.1f min of fitting in total\n",
                sum(unlist(res[secs]), na.rm = TRUE) / 60))
    print(as.data.frame(timing), row.names = FALSE, digits = 3)
  }
  
  #4. check figure: C on top, then IBS, then D-calibration
  prim <- summ[summ$arm == names(arms)[1], ]
  lat  <- prim[!prim$is_anchor & !is.na(prim$delta_c), ]
  anc  <- prim[prim$is_anchor, ]
  byk <- summ %>%
    dplyr::group_by(arm, n_vars) %>%
    dplyr::summarise(delta_c = mean(delta_c, na.rm = TRUE),
                     ibs_gain = -mean(delta_ibs, na.rm = TRUE),
                     dcal_gain = -mean(delta_dcal, na.rm = TRUE), .groups = "drop") %>%
    dplyr::filter(is.finite(delta_c))
  ANCHOR_COL <- "#e8871a"; TREND_COL <- "#DC143C"
  scatter <- function(y, ysd, x, xlab, ylab, anc_x, anc_y) {
    p <- ggplot(lat, aes(.data[[x]], .data[[y]])) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
      geom_errorbar(aes(ymin = .data[[y]] - .data[[ysd]], ymax = .data[[y]] + .data[[ysd]]),
                    width = 0, colour = "grey80", linewidth = 0.4)
    if (nrow(lat) > 2)
      p <- p + geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                           colour = TREND_COL, linewidth = 0.9)
    p <- p + geom_point(aes(fill = n_vars), shape = 21, size = 2.4, colour = "white", stroke = 0.3)
    if (nrow(anc) && is.finite(anc_y))
      p <- p + annotate("point", x = anc_x, y = anc_y, shape = 23, size = 4,
                        fill = ANCHOR_COL, colour = "black")
    p + scale_fill_viridis_c(name = "# clinical\nvariables") +
      labs(x = xlab, y = ylab) +
      theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank())
  }
  lines <- function(y, ylab) {
    ggplot(byk, aes(n_vars, .data[[y]], colour = arm)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
      geom_line(linewidth = 0.8) + geom_point(size = 1.8) +
      scale_x_continuous(breaks = 0:max(summ$n_vars)) +
      labs(x = "Clinical variables in the baseline (0 = omics-only anchor)", y = ylab,
           colour = "arm") +
      theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank())
  }
  lat$ibs_gain <- -lat$delta_ibs; lat$dcal_gain <- -lat$delta_dcal
  pA <- scatter("delta_c", "delta_c_sd", "c_ref", "Comparator Uno's C",
                expression(paste(Delta, "C (up = RNA helps)")),
                if (nrow(anc)) anc$c_ref else NA, if (nrow(anc)) anc$delta_c else NA)
  pB <- lines("delta_c", expression(paste("Mean ", Delta, "C")))
  pC <- scatter("ibs_gain", "delta_ibs_sd", "ibs_ref", "Comparator IBS (lower = better)",
                "IBS gain (up = RNA helps)",
                if (nrow(anc)) anc$ibs_ref else NA, if (nrow(anc)) -anc$delta_ibs else NA)
  pD <- lines("ibs_gain", "Mean IBS gain")
  pE <- scatter("dcal_gain", "delta_dcal_sd", "dcal_ref", "Comparator D-calibration (lower = better)",
                "D-calibration gain (up = RNA helps)",
                if (nrow(anc)) anc$dcal_ref else NA, if (nrow(anc)) -anc$delta_dcal else NA)
  pF <- lines("dcal_gain", "Mean D-calibration gain")
  fig <- (pA + pB) / (pC + pD) / (pE + pF) + plot_annotation(
    title = sprintf("%s | %d subsets + omics-only anchor | %d repeats x %d-fold CV | left: %s arm, bars = SD over repeats",
                    key, nrow(lat), n_rep, n_folds, names(arms)[1]),
    theme = theme(plot.title = element_text(face = "bold", size = 11)))
  f_fig <- if (exists("path_figure", mode = "function")) path_figure(key, "summary_check") else
    file.path(if (exists("DIR_FIGURES")) DIR_FIGURES else DIR_RESULTS,
              paste0(key, "_summary_check.png"))
  f_fig <- sfx(f_fig)
  ggsave(f_fig, fig, width = 13, height = 13, dpi = 200)
  
  #5. optional .xlsx
  if (write_xlsx && requireNamespace("writexl", quietly = TRUE)) {
    sheets <- list(subset_summary = summ)
    if (!is.null(timing)) sheets$timing <- timing
    writexl::write_xlsx(sheets, sub("\\.csv$", ".xlsx", f_summ))
  }
  
  cat("\nSaved ", paste(c(f_summ, f_rep, if (!is.null(timing)) f_time, f_fig), collapse = ", "),
      "\n", sep = "")
  invisible(list(summary = summ, timing = timing))
}