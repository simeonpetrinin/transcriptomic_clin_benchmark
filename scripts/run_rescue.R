#Rescue run: rerun fused (clinical + RNA) cells with glmnet as in the main run,
#switching to the penalized package only for the glmnet fits that fail.
#Same seeds, folds, preprocessing and relative lambda grid; writes results/<CODE>_rescue.csv per cancer.
#Run from the project root: source("scripts/run_rescue.R"). Needs install.packages("penalized").

source("paths.R"); source_library()
suppressMessages(library(penalized))

CODES        <- c("SARC", "COAD", "CESC", "READ", "UCEC", "KIRP", "BRCA")  #smallest first; KIRC already done
CELLS        <- "failed_short"  #"failed" = cells recorded as failed, "failed_short" = also short-path cells, "all"
START_CHECK  <- "rescue"  #"off", "flag" or "rescue": compare glmnet's first path step with coxph
SHORT_RESCUE <- TRUE      #TRUE = a glmnet path shorter than NLAMBDA steps counts as failed
START_TOL    <- 0.05      #largest allowed clinical coefficient gap for the start check
PEN_MAXITER  <- 1000      #iteration cap per penalized path step
RESUME       <- TRUE      #TRUE = keep folds already in a cancer's rescue file and skip them
FR <- exp(seq(0, log(LAMBDA_MIN_RATIO), length.out = NLAMBDA))  #same relative grid as .fit_cox_path

#penalized along the same relative path (lambda1 relative to the first gene's entry point)
pen_path <- function(X, y, pf, fr) {
  ci <- which(pf == 0); n <- nrow(X)
  cph <- coef(survival::coxph(y ~ X[, ci], ties = "breslow")); cph[!is.finite(cph)] <- 0
  r0 <- glmnet::coxgrad(drop(X[, ci] %*% cph), y, rep(1, n), std.weights = FALSE, cox.ties = "breslow")
  l1max <- max(abs(drop(crossprod(X[, -ci, drop = FALSE], r0))))
  B <- matrix(NA_real_, ncol(X), length(fr)); sb <- rep(0, ncol(X) - length(ci)); sg <- cph
  for (j in seq_along(fr)) {
    f <- penalized(y, penalized = X[, -ci, drop = FALSE], unpenalized = X[, ci, drop = FALSE],
                   lambda1 = l1max * fr[j], model = "cox", standardize = FALSE, trace = FALSE,
                   startbeta = sb, startgamma = sg, maxiter = PEN_MAXITER)
    if (!isTRUE(f@converged)) break  #later steps stay NA, so they are not eligible
    sb <- coefficients(f, "penalized"); sg <- coefficients(f, "unpenalized")
    B[ci, j] <- sg; B[-ci, j] <- sb
  }
  B
}

#glmnet path as in the main run; NULL if it fails (or, with START_CHECK = "rescue", if its start is off)
glmnet_path <- function(X, y, pf, state) {
  g <- tryCatch(.fit_cox_path(X, y, pf, NLAMBDA, LAMBDA_MIN_RATIO), error = function(e) NULL)
  if (is.null(g)) { state$why <- c(state$why, "error"); return(NULL) }
  if (length(g$lambda) < NLAMBDA) {
    state$short_flag <- TRUE
    if (SHORT_RESCUE) { state$why <- c(state$why, "short"); return(NULL) }
  }
  if (START_CHECK == "off") return(g)
  ci <- which(pf == 0)
  ref <- tryCatch(coef(survival::coxph(y ~ X[, ci], ties = COX_TIES)), error = function(e) NULL)
  b1 <- as.matrix(g$beta)[ci, 1]
  use <- !is.null(ref) && all(is.finite(ref))
  bad <- use && any(abs(b1 - ref)[abs(ref) < 10] > START_TOL)  #ignore separated (near-infinite) terms
  if (bad) state$start_flag <- TRUE
  if (bad && START_CHECK == "rescue") { state$why <- c(state$why, "start"); return(NULL) }
  g
}

#one cell: clinical-only, then fused with inner tuning and final fit, penalized where glmnet fails
run_rescue_cell <- function(vars, sid, r, seed_r, f, fo) {
  t0 <- proc.time()[["elapsed"]]; state <- new.env(); state$start_flag <- FALSE; state$short_flag <- FALSE; state$why <- character(0); rescued <- character(0)
  trc <- fo$trc; tec <- fo$tec; y <- survival::Surv(trc$OS.time, trc$OS)
  cl <- build_clinical(vars, trc, tec); fx <- fuse(cl, fo$rn)
  m1 <- fit_clinical(cl$train, trc$OS.time, trc$OS, cl$test, fo$grid)
  c1 <- uno_c(tec$OS.time, tec$OS, m1$risk, fo$tau)
  i1 <- integrated_brier(tec$OS.time, tec$OS, m1$surv, fo$grid, fo$tau, cens_time = trc$OS.time, cens_status = trc$OS)
  res <- tryCatch({
    #inner tuning: same loss as tune_clin_rna
    loss <- matrix(NA_real_, INNER_FOLDS, NLAMBDA); vn <- integer(INNER_FOLDS)
    for (k in seq_len(INNER_FOLDS)) {
      tri <- which(fo$fid != k); vai <- which(fo$fid == k)
      ic <- build_clinical(vars, trc[tri, , drop = FALSE], trc[vai, , drop = FALSE])
      ix <- fuse(ic, fo$inner[[k]]$rna)
      all_x <- rbind(as.matrix(ix$train), as.matrix(ix$test))
      g <- glmnet_path(as.matrix(ix$train), y[tri, ], ix$penalty, state)
      if (!is.null(g)) {
        eta <- as.matrix(predict(g, newx = all_x, type = "link"))
      } else {
        rescued <- c(rescued, paste0("inner", k))
        B <- pen_path(as.matrix(ix$train), y[tri, ], ix$penalty, FR)
        eta <- all_x %*% B[, colSums(is.finite(B)) == nrow(B), drop = FALSE]
      }
      np <- ncol(eta); if (!np || any(!is.finite(eta))) next
      d_all <- glmnet::coxnet.deviance(pred = eta, y = y[c(tri, vai), ], std.weights = FALSE, cox.ties = COX_TIES)
      d_tr  <- glmnet::coxnet.deviance(pred = eta[seq_along(tri), , drop = FALSE], y = y[tri, ],
                                       std.weights = FALSE, cox.ties = COX_TIES)
      loss[k, seq_len(np)] <- (d_all - d_tr) / length(vai); vn[k] <- length(vai)
    }
    elig <- which(colSums(is.finite(loss)) == INNER_FOLDS)
    if (!length(elig)) stop("No relative lambda has a finite score in every inner fold.")
    cvm <- colSums(loss[, elig, drop = FALSE] * vn) / sum(vn); sel <- elig[which.min(cvm)]
    #final fit on the outer training set
    Xtr <- as.matrix(fx$train); Xte <- as.matrix(fx$test)
    g <- glmnet_path(Xtr, y, fx$penalty, state)
    fin <- if (!is.null(g) && sel <= length(g$lambda)) tryCatch({
      lam <- g$lambda[sel]
      sf <- survival::survfit(g, s = lam, x = Xtr, y = y, newx = Xte)
      list(risk = as.numeric(predict(g, newx = Xte, s = lam, type = "link")),
           S = t(summary(sf, times = fo$grid, extend = TRUE)$surv),
           genes = sum(as.numeric(coef(g, s = lam))[fx$penalty > 0] != 0))
    }, error = function(e) NULL) else NULL
    if (is.null(fin)) {
      rescued <- c(rescued, "final")
      b <- pen_path(Xtr, y, fx$penalty, FR[seq_len(sel)])[, sel]
      if (any(!is.finite(b))) stop("penalized did not converge at the selected lambda")
      lt <- drop(Xtr %*% b); lv <- drop(Xte %*% b); et <- sort(unique(trc$OS.time[trc$OS == 1]))
      dH <- sapply(et, function(u) sum(trc$OS[trc$OS.time == u]) / sum(exp(lt[trc$OS.time >= u])))  #Breslow baseline
      H <- sapply(fo$grid, function(u) sum(dH[et <= u]))
      fin <- list(risk = lv, S = exp(-outer(exp(lv), H)), genes = sum(b[fx$penalty > 0] != 0))
    }
    c2 <- uno_c(tec$OS.time, tec$OS, fin$risk, fo$tau)
    i2 <- integrated_brier(tec$OS.time, tec$OS, fin$S, fo$grid, fo$tau, cens_time = trc$OS.time, cens_status = trc$OS)
    list(c2 = c2, i2 = as.numeric(i2), frac = FR[sel], genes = fin$genes, err = NA_character_)
  }, error = function(e) list(c2 = NA_real_, i2 = NA_real_, frac = NA_real_, genes = NA_integer_,
                              err = conditionMessage(e)))
  data.frame(cancer = CODE, subset_id = sid, variables = paste(vars, collapse = "+"), repeat_id = r,
             seed = seed_r, fold = f, rescued = if (length(rescued)) paste(rescued, collapse = ",") else "none",
             start_flag = state$start_flag, short_flag = state$short_flag,
             why = if (length(state$why)) paste(unique(state$why), collapse = ",") else "none", lambda_fraction = res$frac, n_genes = res$genes,
             c_clin = c1, ibs_clin = as.numeric(i1), c_fused = res$c2, ibs_fused = res$i2,
             delta_c = res$c2 - c1, delta_ibs = res$i2 - as.numeric(i1), error = res$err,
             secs = round(proc.time()[["elapsed"]] - t0, 1), stringsAsFactors = FALSE)
}

#workers load the same library once; cancer-specific objects are sent per cancer
cl_w <- if (N_WORKERS > 1) parallel::makeCluster(N_WORKERS) else NULL
if (!is.null(cl_w)) {
  parallel::clusterEvalQ(cl_w, { source("paths.R"); source_library(); suppressMessages(library(penalized)); NULL })
  parallel::clusterExport(cl_w, c("START_CHECK", "START_TOL", "PEN_MAXITER", "FR", "SHORT_RESCUE",
                                  "pen_path", "glmnet_path", "run_rescue_cell"))
}
dir.create(file.path(DIR_RESULTS, "timing"), showWarnings = FALSE, recursive = TRUE)
summ <- list()
t_all <- proc.time()[["elapsed"]]

tryCatch({
  for (CODE in CODES) {
    OUT   <- file.path(DIR_RESULTS, paste0(CODE, "_rescue.csv"))
    OUT_T <- file.path(DIR_RESULTS, "timing", paste0(CODE, "_rescue_fold_timing.csv"))  #per-fold wall time
    dat  <- readRDS(path_modelling(CODE))
    ok   <- !is.na(dat$clinical$OS) & !is.na(dat$clinical$OS.time) & dat$clinical$OS.time > 0
    clin <- dat$clinical[ok, ]; expr <- dat$expr[ok, , drop = FALSE]; rm(dat)
    RARE_MISSING_VARS <<- rare_missing_vars(clin, names(CANCER_VARS[[CODE]]), MIN_LEVEL_N)  #as run_cancer
    
    #which cells to rerun, read from the main results
    rec <- read.csv(path_results(CODE), stringsAsFactors = FALSE)
    rec <- rec[!as.logical(toupper(as.character(rec$is_anchor)) == "TRUE"), ]
    pick <- switch(CELLS,
                   failed       = is.na(rec$c_fused),
                   failed_short = is.na(rec$c_fused) | as.character(rec$fused_path_short) == "TRUE",
                   all          = rep(TRUE, nrow(rec)))
    todo <- rec[pick, c("subset_id", "variables", "repeat_id", "seed", "fold", "c_clin", "c_fused", "delta_c")]
    rm(rec)
    
    #resume: folds already written are skipped; otherwise start this cancer fresh
    done <- character(0)
    if (RESUME && file.exists(OUT)) {
      old <- read.csv(OUT, stringsAsFactors = FALSE); done <- unique(paste(old$repeat_id, old$fold))
    } else {
      for (p in c(OUT, OUT_T)) if (file.exists(p)) file.remove(p)
    }
    keys <- setdiff(unique(paste(todo$repeat_id, todo$fold)), done)
    cat(sprintf("\n===== %s: %d cells to rerun (%s), %d of %d folds left, %d workers =====\n", CODE,
                nrow(todo), CELLS, length(keys), length(unique(paste(todo$repeat_id, todo$fold))), N_WORKERS))
    if (!is.null(cl_w)) parallel::clusterExport(cl_w, c("CODE", "RARE_MISSING_VARS"))
    t_c <- proc.time()[["elapsed"]]
    
    for (key in keys) {
      fold_t0 <- proc.time()[["elapsed"]]
      part <- todo[paste(todo$repeat_id, todo$fold) == key, ]
      r <- part$repeat_id[1]; f <- part$fold[1]; seed_r <- part$seed[1]
      #fold objects built once, exactly as run_one_repeat does
      folds <- make_outer_folds(clin$OS, seed = seed_r); tri <- which(folds != f); tei <- which(folds == f)
      trc <- clin[tri, ]; tec <- clin[tei, ]
      tau <- make_tau(trc$OS.time, trc$OS)
      grid <- sort(unique(quantile(trc$OS.time[trc$OS == 1], GRID_QUANTILES))); grid <- grid[grid <= tau]
      fid <- make_inner_folds(trc$OS, k = INNER_FOLDS, seed = inner_fold_seed(r, f))
      fo <- list(trc = trc, tec = tec, tau = tau, grid = grid, fid = fid,
                 rn = build_rna(expr[tri, , drop = FALSE], expr[tei, , drop = FALSE], TOP_N_GENES),
                 inner = build_inner_rna(expr[tri, , drop = FALSE], fid, TOP_N_GENES))
      job <- function(i) run_rescue_cell(strsplit(part$variables[i], "\\+")[[1]], part$subset_id[i], r, seed_r, f, fo)
      out <- if (is.null(cl_w)) lapply(seq_len(nrow(part)), job) else {
        parallel::clusterExport(cl_w, c("part", "r", "f", "seed_r", "fo"), envir = environment())
        parallel::parLapplyLB(cl_w, seq_len(nrow(part)), job)
      }
      out <- do.call(rbind, out)
      out$rec_c_clin <- part$c_clin; out$rec_c_fused <- part$c_fused; out$rec_delta_c <- part$delta_c
      write.table(out, OUT, sep = ",", row.names = FALSE, col.names = !file.exists(OUT), append = file.exists(OUT))
      tm <- data.frame(repeat_id = r, fold = f, n_cells = nrow(out), n_rescued = sum(out$rescued != "none"),
                       n_failing = sum(is.na(out$c_fused)), workers = N_WORKERS,
                       wall_secs = round(proc.time()[["elapsed"]] - fold_t0, 2), host = Sys.info()[["nodename"]])
      write.table(tm, OUT_T, sep = ",", row.names = FALSE, col.names = !file.exists(OUT_T), append = file.exists(OUT_T))
      cat(sprintf("%s repeat %d fold %d: %d cells, %d rescued, %d still failing, %.1fs (%.0f min this cancer)\n",
                  CODE, r, f, nrow(out), tm$n_rescued, tm$n_failing, tm$wall_secs, (proc.time()[["elapsed"]] - t_c) / 60))
    }
    
    #short summary for this cancer
    res <- read.csv(OUT, stringsAsFactors = FALSE)
    cat(sprintf("\n%s rescue: %d cells | now with a result: %d | still failing: %d | clinical C matches main run: %s\n",
                CODE, nrow(res), sum(!is.na(res$c_fused)), sum(is.na(res$c_fused)),
                isTRUE(all.equal(res$c_clin, res$rec_c_clin, tolerance = 1e-8))))
    cat(sprintf("rescue reasons: %s\n", paste(names(table(res$why)), table(res$why), sep = " = ", collapse = " | ")))
    cat(sprintf("rescued where: inner only %d | final only %d | both %d | glmnet start flagged in %d cells\n",
                sum(grepl("inner", res$rescued) & !grepl("final", res$rescued)),
                sum(res$rescued == "final"), sum(grepl("inner", res$rescued) & grepl("final", res$rescued)),
                sum(res$start_flag)))
    cat(sprintf("median lift in C of rerun cells: %.4f | written to %s and %s\n",
                median(res$delta_c, na.rm = TRUE), OUT, OUT_T))
    summ[[CODE]] <- data.frame(cancer = CODE, cells = nrow(res), with_result = sum(!is.na(res$c_fused)),
                               still_failing = sum(is.na(res$c_fused)), core_min = round(sum(res$secs) / 60, 1),
                               wall_min_this_call = round((proc.time()[["elapsed"]] - t_c) / 60, 1))
  }
}, finally = if (!is.null(cl_w)) parallel::stopCluster(cl_w))

#all cancers in one table
summ <- do.call(rbind, summ)
write.csv(summ, file.path(DIR_RESULTS, "timing", "rescue_summary.csv"), row.names = FALSE)
cat("\n===== rescue summary =====\n"); print(summ, row.names = FALSE)
cat(sprintf("total wall time %.0f min\n", (proc.time()[["elapsed"]] - t_all) / 60))