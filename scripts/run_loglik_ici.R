#Log-likelihood and ICI
#Partial-likelihood ratio per death and ICI (median, lower and upper time points) for every
#cell of the lattice, clinical and fused, computed from the saved <CODE>_preds.rds without
#refitting, on N_WORKERS workers (settings.txt). One cancer at a time; writes
#<CODE>_loglik_ici.csv. Run after the lattice.

#ICI time points on the IBS grid (whole thousandths, so no floating-point noise): the median,
#and the symmetric pair closest to the quartiles without lying outside them
.llici_points <- function(gq) {
  g <- round(gq * 1000)
  if (anyNA(g)) stop("the IBS grid contains missing values")
  if (!500 %in% g) stop("the IBS grid has no median point")
  a <- 500 - g[which(g < 500)]  #distances below the median
  a <- a[which(a <= 250 & (500 + a) %in% g)]  #inside the quartiles, mirrored above the median
  if (length(a) == 0) stop("the IBS grid has no symmetric pair inside the quartiles")
  a <- max(a)
  c(q50 = which(g == 500)[1], qlo = which(g == 500 - a)[1], qhi = which(g == 500 + a)[1])
}
LLICI_Q <- .llici_points(GRID_QUANTILES)

path_llici <- function(code) file.path(DIR_RESULTS, paste0(code, "_loglik_ici.csv"))

#per outer fold: the ICI time points (training event-time quantiles) and the test rows
.llici_fold_info <- function(clin, n_rep) {
  info <- list()
  for (r in seq_len(n_rep)) {
    folds <- make_outer_folds(clin$OS, seed = SEED + r)
    for (f in seq_len(OUTER_FOLDS)) {
      trc <- clin[folds != f, , drop = FALSE]
      tq  <- unname(quantile(trc$OS.time[trc$OS == 1], GRID_QUANTILES))[LLICI_Q]  #as the grid is built
      names(tq) <- names(LLICI_Q)
      info[[paste(r, f)]] <- list(t = tq, test_time = clin$OS.time[folds == f],
                                  test_status = clin$OS[folds == f])
    }
  }
  info
}

#metrics of one saved cell
.llici_cell <- function(x, info) {
  k   <- strsplit(x$key, "|", fixed = TRUE)[[1]]
  r   <- as.integer(sub("^r", "", k[2])); f <- as.integer(sub("^f", "", k[3]))
  sid <- as.integer(sub("^s", "", k[4]))
  fi  <- info[[paste(r, f)]]
  tt  <- x$test_time; st <- x$test_status; n <- length(tt); d <- sum(st)
  cols <- match(fi$t, x$times)  #exact: the grid holds these values
  null <- cox_deviance(tt, st, rep(0, n))
  one <- function(risk, surv) {
    ok  <- length(risk) == n && all(is.finite(risk))
    plr <- if (ok && d > 0) pl_ratio(null, cox_deviance(tt, st, risk), n, d) else NA_real_
    ic  <- vapply(seq_along(cols), function(i)
      if (ok && is.matrix(surv) && nrow(surv) == n && !is.na(cols[i])) ici(tt, st, surv[, cols[i]], fi$t[i])
      else NA_real_, numeric(1))
    c(plr, ic)
  }
  a <- one(x$risk_clin, x$surv_clin); b <- one(x$risk_fused, x$surv_fused)
  c(sid, r, f, n, d, a, b)
}

#one chunk of cells (runs on a worker)
.llici_chunk <- function(cells, info) do.call(rbind, lapply(cells, .llici_cell, info = info))

#one cancer end to end
loglik_ici_cancer <- function(code, workers = N_WORKERS, force = FALSE) {
  if (file.exists(path_llici(code)) && !force) {
    cat(sprintf("[loglik_ici] %-5s skipped (%s exists)\n", code, basename(path_llici(code)))); return(invisible(FALSE))
  }
  ppath <- sub("_results\\.csv$", "_preds.rds", path_results(code))
  if (!file.exists(ppath)) stop("no saved predictions for ", code, " (", ppath, ")")
  dat  <- readRDS(path_modelling(code))
  ok   <- !is.na(dat$clinical$OS) & !is.na(dat$clinical$OS.time) & dat$clinical$OS.time > 0
  clin <- dat$clinical[ok, ]; rm(dat)  #same rows, same order as run_cancer()
  info <- .llici_fold_info(clin, n_repeats_for(code))
  
  t0 <- .now(); invisible(gc(reset = TRUE))
  prd <- readRDS(ppath)
  cat(sprintf("=== log-likelihood + ICI %s: %d cells | predictions %.1f GB in memory | loaded in %.0fs ===\n",
              code, length(prd), as.numeric(object.size(prd)) / 1024^3, .now() - t0))
  cat(sprintf("  ICI at percentiles %s of training death times\n",
              paste(sprintf("%s = %.1f", names(LLICI_Q), 100 * GRID_QUANTILES[LLICI_Q]), collapse = ", ")))
  
  #every cell's test rows must be the rebuilt fold's test rows
  for (x in prd[!duplicated(sub("\\|s[0-9]+$", "", vapply(prd, `[[`, "", "key")))]) {
    k <- strsplit(x$key, "|", fixed = TRUE)[[1]]
    fi <- info[[paste(sub("^r", "", k[2]), sub("^f", "", k[3]))]]
    if (!identical(as.numeric(x$test_time), as.numeric(fi$test_time)) ||
        !identical(as.numeric(x$test_status), as.numeric(fi$test_status)))
      stop(sprintf("%s %s: saved test rows differ from the rebuilt fold", code, x$key))
    if (anyNA(match(fi$t, x$times))) stop(sprintf("%s %s: ICI time points not on the saved grid", code, x$key))
  }
  
  if (workers > 1) {
    cl <- parallel::makeCluster(workers)
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    parallel::clusterEvalQ(cl, { source("paths.R"); source_library(); NULL })
    parallel::clusterExport(cl, c(".llici_cell", ".llici_chunk"), envir = environment())
    chunks <- split(prd, cut(seq_along(prd), workers, labels = FALSE))
    rm(prd); invisible(gc())
    out <- do.call(rbind, parallel::parLapply(cl, chunks, .llici_chunk, info = info))
    rm(chunks)
  } else {
    out <- .llici_chunk(prd, info); rm(prd)
  }
  invisible(gc())
  
  q <- names(LLICI_Q)
  colnames(out) <- c("subset_id", "repeat_id", "fold", "n_test", "n_deaths",
                     "plr_clin", paste0("ici_clin_", q), "plr_fused", paste0("ici_fused_", q))
  out <- as.data.frame(out)
  out$plr_lift <- out$plr_fused / out$plr_clin  #>1 = RNA puts more probability on who died
  for (qq in q) out[[paste0("delta_ici_", qq)]] <- out[[paste0("ici_fused_", qq)]] - out[[paste0("ici_clin_", qq)]]
  out <- cbind(cancer = code, out[order(out$subset_id, out$repeat_id, out$fold), ])
  write.csv(out, path_llici(code), row.names = FALSE)
  cat(sprintf("Saved %s (%d cells, %.1f min, peak memory %.1f GB)\n", path_llici(code), nrow(out),
              (.now() - t0) / 60, .peak_mb() / 1024))
  invisible(out)
}

#all cancers, one at a time; one failing does not stop the rest
full_run_loglik_ici <- function(codes = CANCERS, workers = N_WORKERS, force = FALSE) {
  failed <- character(0)
  for (code in codes) {
    ok <- tryCatch({ timed_stage(code, "loglik_ici", loglik_ici_cancer(code, workers, force)); TRUE },
                   error = function(e) { cat(sprintf("[FAILED]     %-5s loglik_ici: %s\n", code, conditionMessage(e))); FALSE })
    if (!ok) failed <- c(failed, code)
    invisible(gc())
  }
  invisible(list(failed = failed))
}