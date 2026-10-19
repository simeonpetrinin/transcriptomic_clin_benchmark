#Baseline selection
#Selection-corrected best and worst clinical subsets per cancer (Varma & Simon 2006).
#Inside each outer training fold, every subset is scored by inner CV on the fused arm's
#inner folds by partial likelihood, C, IBS and ICI; the best and worst subset under each
#metric are then scored on the outer test fold (one figure panel per metric).

SEL_METRICS       <- c(pl = TRUE, c = FALSE, ibs = TRUE, ici = TRUE)  #TRUE = lower is better
SEL_ICI_MIN_FOLDS <- 3L  #inner folds with a usable ICI needed to rank a subset by ICI

path_sel       <- function(code) file.path(DIR_RESULTS, paste0(code, "_baseline_selection.csv"))
path_sel_inner <- function(code) file.path(DIR_RESULTS, paste0(code, "_baseline_inner_scores.csv.gz"))
path_sel_ckpt  <- function(code, r) file.path(DIR_CHECKPOINTS, sprintf("%s_selection_repeat_%02d.rds", code, r))

#clinical Cox as in fit_clinical(), without the unused standard errors
.sel_fit <- function(Xtr, time_tr, status_tr, Xte, times) {
  tryCatch(withCallingHandlers({
    dtr <- data.frame(.t = time_tr, .s = status_tr, Xtr, check.names = TRUE)
    dte <- data.frame(Xte, check.names = TRUE)
    fit <- coxph(Surv(.t, .s) ~ ., data = dtr, ties = COX_TIES, model = TRUE)
    eta_tr <- as.numeric(predict(fit, type = "lp"))
    eta_te <- as.numeric(predict(fit, newdata = dte, type = "lp"))
    S <- t(summary(survfit(fit, newdata = dte, se.fit = FALSE), times = times, extend = TRUE)$surv)
    if (is.null(dim(S))) S <- matrix(S, nrow = nrow(Xte))
    if (any(!is.finite(c(eta_tr, eta_te)))) stop("non-finite risk")
    list(eta_tr = eta_tr, eta_te = eta_te, surv = S)
  }, warning = function(w) invokeRestart("muffleWarning")),
  error = function(e) NULL)
}

#times a fold needs: IBS grid (as in cv.R) and the ICI time (median training event time)
.sel_times <- function(tr) {
  tau  <- make_tau(tr$OS.time, tr$OS)
  grid <- sort(unique(quantile(tr$OS.time[tr$OS == 1], GRID_QUANTILES)))
  grid <- grid[grid <= tau]
  t0   <- unname(quantile(tr$OS.time[tr$OS == 1], 0.5))
  list(tau = tau, grid = grid, t0 = t0, pg = sort(unique(c(grid, t0))))
}

#score one fitted model on held-out rows -> c(C, IBS, ICI)
.sel_scores <- function(m, te, tr, tm) {
  S <- m$surv
  c(c   = uno_c(te$OS.time, te$OS, m$eta_te, tm$tau),
    ibs = as.numeric(integrated_brier(te$OS.time, te$OS, S[, match(tm$grid, tm$pg), drop = FALSE],
                                      tm$grid, tm$tau, cens_time = tr$OS.time, cens_status = tr$OS)),
    ici = ici(te$OS.time, te$OS, S[, match(tm$t0, tm$pg)], tm$t0))
}

#inner scores of one subset: partial likelihood (lower better), C, IBS, ICI
sel_inner_score <- function(subset, trc, foldid) {
  fail <- c(pl = NA_real_, c = NA_real_, ibs = NA_real_, ici = NA_real_)
  ids <- sort(unique(foldid))
  sc <- matrix(NA_real_, length(ids), 3, dimnames = list(NULL, c("c", "ibs", "ici")))
  loss <- numeric(length(ids)); nv <- integer(length(ids))
  for (k in ids) {
    tri <- which(foldid != k); vai <- which(foldid == k)
    itr <- trc[tri, , drop = FALSE]; iva <- trc[vai, , drop = FALSE]
    cl <- build_clinical(subset, itr, iva)
    if (ncol(cl$train) == 0) return(fail)
    tm <- .sel_times(itr)
    if (length(tm$grid) < 2) return(fail)
    m <- .sel_fit(cl$train, itr$OS.time, itr$OS, cl$test, tm$pg)
    if (is.null(m)) return(fail)  #a failed inner fit removes the subset from every ranking
    y   <- Surv(trc$OS.time[c(tri, vai)], trc$OS[c(tri, vai)])
    eta <- matrix(c(m$eta_tr, m$eta_te), ncol = 1)
    d_all <- glmnet::coxnet.deviance(pred = eta, y = y, std.weights = FALSE, cox.ties = COX_TIES)
    d_tr  <- glmnet::coxnet.deviance(pred = eta[seq_along(tri), , drop = FALSE], y = y[seq_along(tri), ],
                                     std.weights = FALSE, cox.ties = COX_TIES)
    loss[k] <- (d_all - d_tr) / length(vai); nv[k] <- length(vai)  #as in tune_clin_rna()
    sc[k, ] <- .sel_scores(m, iva, itr, tm)
  }
  pl <- sum(loss * nv) / sum(nv)
  n_ici <- sum(is.finite(sc[, "ici"]))
  c(pl  = if (is.finite(pl)) pl else Inf,  #separated fits rank last, not excluded
    c   = if (all(is.finite(sc[, "c"])))   mean(sc[, "c"])   else NA_real_,
    ibs = if (all(is.finite(sc[, "ibs"]))) mean(sc[, "ibs"]) else NA_real_,
    ici = if (n_ici >= SEL_ICI_MIN_FOLDS) mean(sc[is.finite(sc[, "ici"]), "ici"]) else NA_real_)
}

#worker task: score subset number si with this fold's objects
.sel_on_worker <- function(si) with(sel_env, sel_inner_score(subset_list[[si]], trc, foldid))

#pick one subset; ties go to the lower subset_id (fewer variables)
.sel_pick <- function(x, ids, lower_is_best) {
  ok <- which(!is.na(x))
  if (!length(ok)) return(NA_integer_)
  ids[ok][if (lower_is_best) which.min(x[ok]) else which.max(x[ok])]
}

#outer test-fold scores of one pick: C and IBS from the results file, PL and ICI from a refit
.sel_outer <- function(sid, prefix, subset_list, res_rf, trc, tec) {
  cols <- c("c_clin", "ibs_clin", "c_fused", "ibs_fused", "delta_c", "delta_ibs")
  out <- setNames(as.list(rep(NA_real_, length(cols) + 2)), c(cols, "pl_clin", "ici_clin"))
  if (!is.na(sid)) {
    row <- res_rf[res_rf$subset_id == sid, , drop = FALSE]
    if (nrow(row) == 1) out[cols] <- as.list(row[1, cols])
    cl <- build_clinical(subset_list[[sid]], trc, tec); tm <- .sel_times(trc)
    m  <- .sel_fit(cl$train, trc$OS.time, trc$OS, cl$test, tm$pg)
    if (!is.null(m)) {
      s <- .sel_scores(m, tec, trc, tm)
      if (is.finite(out$c_clin) && abs(s[["c"]] - out$c_clin) > 1e-6)  #refit must match the lattice
        stop(sprintf("refit C %.6f differs from the results %.6f (subset %d)", s[["c"]], out$c_clin, sid))
      out$ici_clin <- s[["ici"]]
      out$pl_clin  <- cox_deviance(tec$OS.time, tec$OS, m$eta_te)  #test-fold deviance per patient
    }
  }
  names(out) <- paste0(prefix, "_outer_", names(out))
  as.data.frame(out)
}

#one repeat: rebuild the outer and inner folds exactly as cv.R does, rank, score
sel_one_repeat <- function(code, clin, res, subset_list, subset_ids, labels, r, seed_r, cl_workers) {
  folds <- make_outer_folds(clin$OS, seed = seed_r)
  rows <- list(); inner <- list()
  for (f in seq_len(OUTER_FOLDS)) {
    res_rf <- res[res$repeat_id == r & res$fold == f, , drop = FALSE]
    if (!nrow(res_rf)) { cat(sprintf("  r%d f%d: not in results (skipped in the lattice)\n", r, f)); next }
    trc <- clin[folds != f, , drop = FALSE]; tec <- clin[folds == f, , drop = FALSE]
    if (!all(res_rf$n_test == nrow(tec)))  #fold rebuild must match the lattice
      stop(sprintf("%s r%d f%d: rebuilt test fold size differs from the results", code, r, f))
    foldid <- make_inner_folds(trc$OS, k = INNER_FOLDS, seed = inner_fold_seed(r, f))
    t0 <- .now()
    sc <- if (is.null(cl_workers)) {
      lapply(seq_along(subset_list), function(si) sel_inner_score(subset_list[[si]], trc, foldid))
    } else {
      sel_env <- list(subset_list = subset_list, trc = trc, foldid = foldid)
      parallel::clusterExport(cl_workers, "sel_env", envir = environment())
      parallel::parLapplyLB(cl_workers, seq_along(subset_list), .sel_on_worker)
    }
    sc <- do.call(rbind, sc)
    picks <- integer(0); crit <- character(0)
    for (mt in names(SEL_METRICS)) for (end in c("best", "worst")) {
      lower <- if (end == "best") SEL_METRICS[[mt]] else !SEL_METRICS[[mt]]
      picks[paste0(end, "_", mt)] <- .sel_pick(sc[, mt], subset_ids, lower)
      crit[paste0(end, "_", mt)]  <- mt
    }
    base <- data.frame(cancer = code, repeat_id = r, fold = f,
                       n_ranked_pl = sum(!is.na(sc[, "pl"])), n_ranked_ici = sum(!is.na(sc[, "ici"])),
                       outer_n_test = nrow(tec), outer_n_deaths = sum(tec$OS),  #for the per-death ratio
                       outer_pl_null = cox_deviance(tec$OS.time, tec$OS, rep(0, nrow(tec))),  #no-predictor deviance
                       secs = NA_real_, stringsAsFactors = FALSE)
    for (p in names(picks)) {
      base[[paste0(p, "_subset_id")]] <- picks[[p]]
      base[[paste0(p, "_variables")]] <- labels[match(picks[[p]], subset_ids)]
      base[[paste0(p, "_inner")]]     <- sc[match(picks[[p]], subset_ids), crit[[p]]]
    }
    outer <- do.call(cbind, lapply(names(picks), function(p)
      .sel_outer(picks[[p]], p, subset_list, res_rf, trc, tec)))
    base$secs <- round(.now() - t0, 1)
    rows[[length(rows) + 1]] <- cbind(base, outer)
    inner[[length(inner) + 1]] <- data.frame(repeat_id = r, fold = f, subset_id = subset_ids,
                                             inner_pl = sc[, "pl"], inner_c = sc[, "c"],
                                             inner_ibs = sc[, "ibs"], inner_ici = sc[, "ici"])
    cat(sprintf("  r%d f%d: best by PL %s | best by C %s | %.0fs\n", r, f,
                base$best_pl_variables, base$best_c_variables, base$secs))
  }
  list(rows = do.call(rbind, rows), inner = do.call(rbind, inner))
}

#one cancer end to end, with a checkpoint per repeat
select_baselines_cancer <- function(code, force = FALSE) {
  if (!file.exists(path_results(code))) stop("no lattice results for ", code, " -- run the lattice first")
  if (file.exists(path_sel(code)) && !force) {
    cat(sprintf("[selection]  %-5s skipped (%s exists)\n", code, basename(path_sel(code)))); return(invisible(FALSE))
  }
  dat  <- readRDS(path_modelling(code))
  ok   <- !is.na(dat$clinical$OS) & !is.na(dat$clinical$OS.time) & dat$clinical$OS.time > 0
  clin <- dat$clinical[ok, ]  #same rows, same order as run_cancer()
  res  <- read.csv(path_results(code), stringsAsFactors = FALSE)
  res  <- res[!as.logical(res$is_anchor), , drop = FALSE]
  
  subset_list <- all_subsets(names(CANCER_VARS[[code]]))  #anchor left out on purpose
  subset_ids  <- seq_along(subset_list)
  labels      <- vapply(subset_list, paste, character(1), collapse = "+")
  RARE_MISSING_VARS <<- rare_missing_vars(clin, unique(unlist(subset_list)), MIN_LEVEL_N)
  n_rep <- n_repeats_for(code); seeds <- SEED + seq_len(n_rep)
  cat(sprintf("=== selection %s: %d subsets x %d folds x %d repeats ===\n",
              code, length(subset_list), OUTER_FOLDS, n_rep))
  
  cl_workers <- NULL
  if (N_WORKERS > 1) {  #core count from settings.txt
    cl_workers <- parallel::makeCluster(N_WORKERS)
    on.exit(try(parallel::stopCluster(cl_workers), silent = TRUE), add = TRUE)
    parallel::clusterEvalQ(cl_workers, { source("paths.R"); source_library(); NULL })
    parallel::clusterExport(cl_workers, c("sel_inner_score", ".sel_on_worker", ".sel_fit", ".sel_times",
                                          ".sel_scores", "SEL_ICI_MIN_FOLDS", "RARE_MISSING_VARS",
                                          "CLINICAL_MISSING", "MISSING_LEVEL", "COX_TIES",
                                          "GRID_QUANTILES", "G_FLOOR"),
                            envir = environment())
  }
  
  out <- vector("list", n_rep)
  for (r in seq_len(n_rep)) {
    ck <- path_sel_ckpt(code, r)
    if (file.exists(ck) && !force) { out[[r]] <- readRDS(ck); next }
    out[[r]] <- sel_one_repeat(code, clin, res, subset_list, subset_ids, labels, r, seeds[r], cl_workers)
    saveRDS(out[[r]], ck)
  }
  sel <- do.call(rbind, lapply(out, `[[`, "rows"))
  inn <- do.call(rbind, lapply(out, `[[`, "inner"))
  write.csv(sel, path_sel(code), row.names = FALSE)
  gz <- gzfile(path_sel_inner(code), "w"); write.csv(inn, gz, row.names = FALSE); close(gz)
  cat(sprintf("Saved %s (%d folds) and %s\n", path_sel(code), nrow(sel), path_sel_inner(code)))
  invisible(sel)
}

#all cancers; one failing does not stop the rest
full_run_baseline_selection <- function(codes = CANCERS, force = FALSE) {
  failed <- character(0)
  for (code in codes) {
    ok <- tryCatch({ timed_stage(code, "baseline_selection", select_baselines_cancer(code, force)); TRUE },
                   error = function(e) { cat(sprintf("[FAILED]     %-5s selection: %s\n", code, conditionMessage(e))); FALSE })
    if (!ok) failed <- c(failed, code)
    invisible(gc())
  }
  invisible(list(failed = failed))
}