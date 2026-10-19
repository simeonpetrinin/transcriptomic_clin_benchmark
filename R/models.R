#Models
#The fitters: Kaplan-Meier, clinical Cox, fused Cox (lasso or ridge), offset, RNA score,
#random survival forest, blockForest and CoxBoost. Each returns risk (higher = more
#hazard) and a patients x times survival matrix; a failed fit returns NA fields.

suppressPackageStartupMessages(library(survival))

#every fitter returns train_secs: inner tuning + final fit, not prediction or scoring
.now <- function() proc.time()[["elapsed"]]

LAMBDA_MIN_RATIO <- 0.05  #path floor = LAMBDA_MIN_RATIO * lambda_max
NLAMBDA          <- 100  #lambdas along the path, standard in glmnet
INNER_FOLDS      <- 5  #inner CV folds for lambda selection
RIDGE_LAMBDA_MIN_RATIO <- 0.01

#tie handling, set explicitly for coxph and glmnet so both arms use the same one
COX_TIES <- "efron"  #"breslow" or "efron"

#|beta| at or above this marks a clinical fit as separated (a declared convention)
COEF_SEPARATION <- 20

#Kaplan-Meier reference (the null survival model/anchor)
fit_km <- function(time_tr, status_tr, n_test, pred_times) {
  t0 <- .now()
  tryCatch({
    sf   <- survfit(Surv(time_tr, status_tr) ~ 1)
    train_secs <- .now() - t0
    s_at <- summary(sf, times = pred_times, extend = TRUE)$surv  #1 x times
    S    <- matrix(rep(s_at, each = n_test), nrow = n_test)  #patients x times
    list(risk = rep(0, n_test), surv = S, times = pred_times, train_secs = round(train_secs, 4))
  }, error = function(e)
    list(risk = NA, surv = NA, times = pred_times, train_secs = NA_real_,
         error = conditionMessage(e)))
}

#clinical-only: plain Cox, with a separation self-report
fit_clinical <- function(Xtr, time_tr, status_tr, Xte, pred_times) {
  warn <- character(0); t0 <- .now()
  tryCatch(
    withCallingHandlers({
      dtr  <- data.frame(.t = time_tr, .s = status_tr, Xtr, check.names = TRUE)
      dte  <- data.frame(Xte, check.names = TRUE)
      #model = TRUE keeps the model frame inside the fit.
      fit  <- coxph(Surv(.t, .s) ~ ., data = dtr, ties = COX_TIES, model = TRUE)
      train_secs <- .now() - t0
      risk <- as.numeric(predict(fit, newdata = dte, type = "lp"))
      sf   <- survfit(fit, newdata = dte)
      S    <- t(summary(sf, times = pred_times, extend = TRUE)$surv)  #patients x times
      if (is.null(dim(S))) S <- matrix(S, nrow = nrow(Xte))
      
      cf       <- stats::coef(fit)
      finite   <- is.finite(cf)
      max_abs  <- if (any(finite)) max(abs(cf[finite])) else NA_real_
      n_bad    <- sum(!finite)
      n_iter   <- if (is.null(fit$iter)) NA_integer_ else as.integer(fit$iter[length(fit$iter)])
      warn_txt <- if (length(warn)) paste(unique(warn), collapse = " | ") else NA_character_
      sep_warn <- any(grepl("infinite|did not converge|ran out of iterations",
                            warn, ignore.case = TRUE))
      separated <- isTRUE(sep_warn) || n_bad > 0 ||
        (is.finite(max_abs) && max_abs >= COEF_SEPARATION)
      
      list(risk = risk, surv = S, times = pred_times, train_secs = round(train_secs, 4),
           n_coef = length(cf), max_abs_coef = max_abs, n_bad_coef = n_bad,
           n_iter = n_iter, separated = separated, warning = warn_txt,
           error = NA_character_)
    },
    warning = function(w) {
      warn <<- c(warn, conditionMessage(w))
      invokeRestart("muffleWarning")
    }),
    error = function(e)
      list(risk = NA, surv = NA, times = pred_times, train_secs = NA_real_,
           n_coef = NA_integer_, max_abs_coef = NA_real_, n_bad_coef = NA_integer_,
           n_iter = NA_integer_, separated = NA, warning = NA_character_,
           error = conditionMessage(e)))
}

#stratified inner folds, deterministic from a supplied seed
make_inner_folds <- function(status, k = INNER_FOLDS, seed) {
  n_ev <- sum(status == 1)
  k    <- max(3L, min(k, n_ev))  #glmnet requires nfolds >= 3
  if (n_ev < k)
    stop("Insufficient events to place at least one event in every inner fold.")  #debugging comment
  old  <- if (exists(".Random.seed", .GlobalEnv)) get(".Random.seed", .GlobalEnv) else NULL
  set.seed(seed)
  fold <- integer(length(status))
  for (lev in unique(status)) {
    idx <- which(status == lev)
    fold[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
  }
  #restore the caller's RNG state: this function must not perturb anything else
  if (!is.null(old)) assign(".Random.seed", old, envir = .GlobalEnv)
  fold
}

#IPF-LASSO fallback: glmnet first; the penalized package (Goeman 2010) refits only the glmnet fits that fail
FUSED_FALLBACK     <- TRUE   #FALSE = glmnet only, failed fits stay failed cells
FALLBACK_SHORT     <- TRUE   #a glmnet path shorter than NLAMBDA counts as failed
FALLBACK_START_TOL <- 0.05   #glmnet's first step must match coxph's clinical fit (|coef| < 10 only)
FALLBACK_MAXITER   <- 1000L  #penalized iteration cap per path step

#glmnet path for the IPF-LASSO arm; NULL (reason kept in state$why) when the fit counts as failed
.fused_glmnet_path <- function(X, y, penalty, nlambda, lambda_min_ratio, state) {
  g <- tryCatch(.fit_cox_path(X, y, penalty, nlambda, lambda_min_ratio), error = function(e) NULL)
  if (is.null(g)) { state$why <- c(state$why, "error"); return(NULL) }
  if (FALLBACK_SHORT && length(g$lambda) < nlambda) { state$why <- c(state$why, "short"); return(NULL) }
  ci <- which(penalty == 0)
  if (length(ci)) {  #start check: catches paths that run but start from a wrong clinical fit
    ref <- tryCatch(stats::coef(survival::coxph(y ~ X[, ci, drop = FALSE], ties = COX_TIES)),
                    error = function(e) NULL)
    b1 <- as.matrix(g$beta)[ci, 1]
    if (!is.null(ref) && all(is.finite(ref)) &&
        any(abs(b1 - ref)[abs(ref) < 10] > FALLBACK_START_TOL)) {
      state$why <- c(state$why, "start"); return(NULL)
    }
  }
  g
}

#penalized lasso path, clinical unpenalised, on the same relative grid; steps after a non-converged one stay NA
.fused_penalized_path <- function(X, y, penalty, fractions) {
  if (!requireNamespace("penalized", quietly = TRUE)) stop("the penalized package is required for the fallback")
  X <- as.matrix(X); ci <- which(penalty == 0); pj <- which(penalty > 0); n <- nrow(X)
  cph <- numeric(0); eta0 <- rep(0, n)
  if (length(ci)) {
    cph <- stats::coef(survival::coxph(y ~ X[, ci, drop = FALSE], ties = "breslow"))
    cph[!is.finite(cph)] <- 0
    eta0 <- drop(X[, ci, drop = FALSE] %*% cph)
  }
  r0 <- glmnet::coxgrad(eta0, y, rep(1, n), std.weights = FALSE, cox.ties = "breslow")
  l1max <- max(abs(drop(crossprod(X[, pj, drop = FALSE], r0))))  #first gene's entry point
  B <- matrix(NA_real_, ncol(X), length(fractions)); sb <- rep(0, length(pj)); sg <- cph
  for (j in seq_along(fractions)) {
    args <- list(response = y, penalized = X[, pj, drop = FALSE], lambda1 = l1max * fractions[j],
                 model = "cox", standardize = FALSE, trace = FALSE, startbeta = sb,
                 maxiter = FALLBACK_MAXITER)
    if (length(ci)) { args$unpenalized <- X[, ci, drop = FALSE]; args$startgamma <- sg }
    f <- do.call(penalized::penalized, args)
    if (!isTRUE(f@converged)) break
    sb <- f@penalized; B[pj, j] <- sb
    if (length(ci)) { sg <- f@unpenalized; B[ci, j] <- sg }
  }
  B
}

#Breslow survival for a fixed linear predictor (used when the final fit is from the fallback)
.breslow_surv <- function(lp_tr, time_tr, status_tr, lp_te, times) {
  et <- sort(unique(time_tr[status_tr == 1]))
  dH <- vapply(et, function(u) sum(status_tr[time_tr == u]) / sum(exp(lp_tr[time_tr >= u])), numeric(1))
  H  <- vapply(times, function(u) sum(dH[et <= u]), numeric(1))
  exp(-outer(exp(lp_te), H))
}

#fit one training-only path on the fixed relative-lambda grid
.fit_cox_path <- function(X, y, penalty, nlambda, lambda_min_ratio,
                          offset = NULL, alpha = 1) {
  X <- as.matrix(X)
  if (ncol(X) < 2L || length(penalty) != ncol(X) || nrow(X) != nrow(y))
    stop("Invalid Cox design dimensions or penalty vector.")
  bad <- which(!vapply(seq_len(ncol(X)), function(j) all(is.finite(X[, j])), logical(1)))
  if (length(bad))
    stop("non-finite values in input column(s): ",
         paste(head(colnames(X)[bad], 5), collapse = ", "))
  if (!is.null(offset)) {
    if (length(offset) != nrow(X) || any(!is.finite(offset)))
      stop("Offset must be finite and one value per training row.")
  }
  #separate glmnet calls on purpose: survfit() re-reads the recorded call
  if (!is.null(offset) && alpha != 1) stop("The offset arm is lasso only (alpha = 1).")
  fit <- if (is.null(offset) && alpha != 1)
    glmnet::glmnet(X, y, family = "cox", alpha = alpha,
                   penalty.factor = penalty, standardize = TRUE,
                   nlambda = nlambda, lambda.min.ratio = lambda_min_ratio,
                   cox.ties = COX_TIES)
  else if (is.null(offset))
    glmnet::glmnet(X, y, family = "cox", alpha = 1,
                   penalty.factor = penalty, standardize = TRUE,
                   nlambda = nlambda, lambda.min.ratio = lambda_min_ratio,
                   cox.ties = COX_TIES)
  else
    glmnet::glmnet(X, y, family = "cox", alpha = 1,
                   penalty.factor = penalty, standardize = TRUE,
                   nlambda = nlambda, lambda.min.ratio = lambda_min_ratio,
                   cox.ties = COX_TIES, offset = offset)
  #glmnet's automatic path uses a geometric lambda grid.
  lam <- fit$lambda
  if (!length(lam) || any(!is.finite(lam)) || any(lam <= 0) || length(lam) > nlambda)
    stop("Cox fit did not return a usable lambda path.")
  fractions <- exp(seq(0, log(lambda_min_ratio), length.out = nlambda))
  if (!isTRUE(all.equal(as.numeric(lam / lam[1]),
                        fractions[seq_along(lam)], tolerance = 1e-6)))
    stop("Returned Cox path does not match the relative-lambda grid.")
  fit
}

#inner tuning: preprocessing rebuilt in every inner split
tune_clin_rna <- function(clinical_subset, clinical_train, inner_rna, foldid,
                          nlambda = NLAMBDA, lambda_min_ratio = LAMBDA_MIN_RATIO,
                          alpha = 1, state = NULL) {
  if (is.null(state)) { state <- new.env(); state$fallback <- character(0); state$why <- character(0) }
  use_fb <- isTRUE(FUSED_FALLBACK) && alpha == 1  #IPF-LASSO only, not the ridge arm
  if (!requireNamespace("glmnet", quietly = TRUE)) stop("glmnet is required.")
  if (!"cox.ties" %in% names(formals(glmnet::glmnet)) ||
      !"cox.ties" %in% names(formals(glmnet::coxnet.deviance)))
    stop("This inner tuner requires glmnet with cox.ties support (project lockfile: 5.0).")
  if (length(nlambda) != 1L || !is.finite(nlambda) || nlambda < 3 || nlambda != as.integer(nlambda) ||
      length(lambda_min_ratio) != 1L || !is.finite(lambda_min_ratio) ||
      lambda_min_ratio <= 0 || lambda_min_ratio >= 1)
    stop("Use nlambda >= 3 and 0 < lambda_min_ratio < 1.")
  if (length(alpha) != 1L || !is.finite(alpha) || alpha < 0 || alpha > 1)
    stop("alpha must be between 0 (ridge) and 1 (lasso).")
  n <- nrow(clinical_train)
  if (length(foldid) != n || anyNA(foldid) || any(foldid != as.integer(foldid))) stop("Inner fold IDs do not match clinical rows.")
  ids <- sort(unique(foldid))
  if (length(ids) < 3L || !identical(as.integer(ids), seq_along(ids)) ||
      length(inner_rna) != length(ids)) stop("Invalid inner RNA cache or fold IDs.")
  if (!all(c("OS.time", "OS", clinical_subset) %in% names(clinical_train)))
    stop("Inner clinical data are missing required columns.")
  y <- survival::Surv(clinical_train$OS.time, clinical_train$OS)
  fractions <- exp(seq(0, log(lambda_min_ratio), length.out = nlambda))
  loss <- matrix(NA_real_, length(ids), nlambda)
  valid_n <- integer(length(ids)); inner_jerr <- integer(length(ids))
  #Per-inner-split clinical width.
  inner_clin_cols <- integer(length(ids))
  inner_dropped   <- vector("list", length(ids))
  for (k in ids) {
    part <- inner_rna[[k]]
    tri <- which(foldid != k); vai <- which(foldid == k)
    if (!identical(part$train_idx, tri) || !identical(part$valid_idx, vai))
      stop("Inner RNA cache row order disagrees with foldid in fold ", k, ".")
    if (!is.null(part$rna$error)) stop("Inner fold ", k, " RNA: ", part$rna$error)
    tryCatch({
      cl <- build_clinical(clinical_subset,
                           clinical_train[tri, , drop = FALSE],
                           clinical_train[vai, , drop = FALSE])
      inner_clin_cols[k] <- ncol(cl$train)
      inner_dropped[[k]] <- cl$dropped
      fx <- fuse(cl, part$rna)
      fit <- if (use_fb) .fused_glmnet_path(fx$train, y[tri, , drop = FALSE], fx$penalty,
                                            nlambda, lambda_min_ratio, state) else
                                              .fit_cox_path(fx$train, y[tri, , drop = FALSE], fx$penalty,
                                                            nlambda, lambda_min_ratio, alpha = alpha)
      all_x <- rbind(fx$train, fx$test)
      if (is.null(fit)) {  #fallback for this inner split
        state$fallback <- c(state$fallback, paste0("inner", k))
        B <- .fused_penalized_path(fx$train, y[tri, , drop = FALSE], fx$penalty, fractions)
        np <- sum(colSums(is.finite(B)) == nrow(B))
        if (!np) stop("fallback: penalized did not converge at the first lambda")
        eta <- as.matrix(all_x) %*% B[, seq_len(np), drop = FALSE]
      } else {
        np <- length(fit$lambda)
        eta <- as.matrix(predict(fit, newx = all_x, type = "link"))
      }
      if (!identical(dim(eta), c(n, np)) || any(!is.finite(eta)))
        stop("Non-finite or incorrectly sized inner predictions.")
      #Validation predictors are only APPLIED to the training-fitted transform.
      d_all <- glmnet::coxnet.deviance(pred = eta, y = y[c(tri, vai), , drop = FALSE],
                                       std.weights = FALSE, cox.ties = COX_TIES)
      d_train <- glmnet::coxnet.deviance(pred = eta[seq_along(tri), , drop = FALSE],
                                         y = y[tri, , drop = FALSE],
                                         std.weights = FALSE, cox.ties = COX_TIES)
      loss[k, seq_len(np)] <- (d_all - d_train) / length(vai)
      valid_n[k] <- length(vai)
      inner_jerr[k] <- if (is.null(fit) || is.null(fit$jerr)) 0L else fit$jerr
    }, error = function(e) stop("Inner fold ", k, ": ", conditionMessage(e), call. = FALSE))
  }
  eligible <- which(colSums(is.finite(loss)) == length(ids))
  if (!length(eligible)) stop("No relative lambda has a finite score in every inner fold.")
  cvm <- rep(NA_real_, nlambda)
  cvm[eligible] <- colSums(loss[, eligible, drop = FALSE] * valid_n) / sum(valid_n)
  #which.min takes the first minimum: ties select the stronger penalty.
  selected <- eligible[which.min(cvm[eligible])]
  dropped_any <- sort(unique(unlist(inner_dropped)))
  list(sel_idx = selected, fraction = fractions[selected], fractions = fractions,
       cvm = cvm, fold_loss = loss, n_valid_path = length(eligible),
       at_boundary = selected >= max(eligible) - 1L, inner_jerr = inner_jerr,
       inner_clin_cols = inner_clin_cols,
       inner_vars_dropped = if (length(dropped_any))
         paste(dropped_any, collapse = "+") else NA_character_)
}

#clinical + RNA : inner tune, then final outer-training fit
fit_clin_rna <- function(Xtr, time_tr, status_tr, Xte, penalty, pred_times,
                         foldid, nlambda = NLAMBDA,
                         lambda_min_ratio = LAMBDA_MIN_RATIO,
                         clinical_subset, clinical_train, inner_rna, alpha = 1) {
  t0 <- .now()
  state <- new.env(); state$fallback <- character(0); state$why <- character(0)  #fallback record
  fb_txt <- function(x) if (length(x)) paste(unique(x), collapse = ",") else "none"
  tryCatch({
    if (missing(clinical_subset) || missing(clinical_train) || missing(inner_rna))
      stop("Nested tuning requires clinical_subset, clinical_train and inner_rna.")
    if (!identical(as.numeric(clinical_train$OS.time), as.numeric(time_tr)) ||
        !identical(as.numeric(clinical_train$OS), as.numeric(status_tr)))
      stop("Final training outcomes disagree with inner clinical row order.")
    tuning <- tune_clin_rna(clinical_subset, clinical_train, inner_rna, foldid,
                            nlambda, lambda_min_ratio, alpha = alpha, state = state)
    Xtr_m <- as.matrix(Xtr); Xte_m <- as.matrix(Xte)
    ytr <- survival::Surv(time_tr, status_tr)
    use_fb <- isTRUE(FUSED_FALLBACK) && alpha == 1
    sel_idx <- tuning$sel_idx
    fit <- if (use_fb) .fused_glmnet_path(Xtr_m, ytr, penalty, nlambda, lambda_min_ratio, state) else
      .fit_cox_path(Xtr_m, ytr, penalty, nlambda, lambda_min_ratio, alpha = alpha)
    train_secs <- .now() - t0  #inner tuning + final path
    if (!is.null(fit) && sel_idx > length(fit$lambda)) {
      if (!use_fb) stop("Final outer-training path did not reach the selected relative lambda.")
      state$why <- c(state$why, "beyond"); fit <- NULL
    }
    fin <- NULL
    if (!is.null(fit)) fin <- tryCatch({
      lam <- fit$lambda[sel_idx]
      sf <- survival::survfit(fit, s = lam, x = Xtr_m, y = ytr, newx = Xte_m)
      S <- t(summary(sf, times = pred_times, extend = TRUE)$surv)
      if (is.null(dim(S))) S <- matrix(S, nrow = nrow(Xte_m))
      list(cf = as.numeric(coef(fit, s = lam)), lam = lam, S = S,
           risk = as.numeric(predict(fit, newx = Xte_m, s = lam, type = "link")),
           n_path = length(fit$lambda), jerr = fit$jerr)
    }, error = function(e) {
      if (!use_fb) stop(e)
      state$why <- c(state$why, "predict"); NULL
    })
    if (is.null(fin)) {  #fallback for the final outer-training fit
      state$fallback <- c(state$fallback, "final")
      B <- .fused_penalized_path(Xtr_m, ytr, penalty, tuning$fractions[seq_len(sel_idx)])
      b <- B[, sel_idx]
      if (any(!is.finite(b))) stop("fallback: penalized did not converge at the selected lambda")
      train_secs <- .now() - t0
      risk_fb <- drop(Xte_m %*% b)
      fin <- list(cf = b, lam = NA_real_, risk = risk_fb, n_path = NA_integer_, jerr = NA,
                  S = .breslow_surv(drop(Xtr_m %*% b), time_tr, status_tr, risk_fb, pred_times))
    }
    n_path <- fin$n_path; lam <- fin$lam; jerr <- fin$jerr; risk <- fin$risk; S <- fin$S
    fit_problem <- !is.null(jerr) && !is.na(jerr) && jerr != 0
    path_end <- if (is.na(n_path)) nlambda else n_path  #fallback: the full relative grid
    at_boundary <- sel_idx >= path_end - 1L
    near_boundary <- sel_idx >= path_end - 4L
    rna_i <- which(penalty > 0)
    keep <- rna_i[fin$cf[rna_i] != 0]
    nm <- colnames(Xtr_m)
    genes <- if (is.null(nm)) as.character(keep) else nm[keep]
    #Outer clinical width = the unpenalised columns of the final design.
    outer_clin_cols <- sum(penalty == 0)
    inner_clin_mismatch <- any(tuning$inner_clin_cols != outer_clin_cols)
    list(risk = risk, surv = S, times = pred_times, lambda = lam,
         train_secs = round(train_secs, 4),
         n_path = n_path, sel_idx = sel_idx, jerr = jerr,
         fit_problem = fit_problem, near_boundary = near_boundary, at_boundary = at_boundary,
         n_genes = length(keep), genes = genes, inner_folds = max(foldid),
         lambda_fraction = tuning$fraction, inner_n_path = tuning$n_valid_path,
         inner_at_boundary = tuning$at_boundary,
         inner_solver_error = any(tuning$inner_jerr != 0),
         outer_clin_cols = outer_clin_cols,
         inner_clin_cols_min = min(tuning$inner_clin_cols),
         inner_clin_cols_max = max(tuning$inner_clin_cols),
         inner_clin_mismatch = inner_clin_mismatch,
         inner_vars_dropped = tuning$inner_vars_dropped,
         fallback = fb_txt(state$fallback), fallback_why = fb_txt(state$why),
         error = NA_character_)
  }, error = function(e)
    list(risk = NA, surv = NA, times = pred_times, lambda = NA, train_secs = NA_real_,
         n_path = NA, sel_idx = NA, jerr = NA, fit_problem = NA,
         near_boundary = NA, at_boundary = NA, n_genes = NA, genes = character(0),
         inner_folds = NA, lambda_fraction = NA, inner_n_path = NA,
         inner_at_boundary = NA, inner_solver_error = NA,
         outer_clin_cols = NA_integer_, inner_clin_cols_min = NA_integer_,
         inner_clin_cols_max = NA_integer_, inner_clin_mismatch = NA,
         inner_vars_dropped = NA_character_,
         fallback = fb_txt(state$fallback), fallback_why = fb_txt(state$why),
         error = conditionMessage(e)))
}

#offset arm (Volkmann et al. 2019): clinical Cox as a fixed offset, lasso on RNA

OFFSET_NLAMBDA          <- NLAMBDA  #kept equal to the fused arm
OFFSET_LAMBDA_MIN_RATIO <- LAMBDA_MIN_RATIO  #so no difference comes from tuning

#clinical linear predictor used as the offset
.clin_offset <- function(Xtr, Xte, time_tr, status_tr) {
  if (ncol(Xtr) == 0L)
    return(list(lp_train = rep(0, nrow(Xtr)), lp_test = rep(0, nrow(Xte)),
                beta = numeric(0), fit = NULL, null_clin = TRUE))
  dtr <- data.frame(.t = time_tr, .s = status_tr, Xtr, check.names = TRUE)
  cf  <- coxph(Surv(.t, .s) ~ ., data = dtr, ties = COX_TIES, model = TRUE)
  b   <- stats::coef(cf)
  #Aliased columns give NA coefficients.
  b[!is.finite(b)] <- 0
  list(lp_train = as.numeric(as.matrix(Xtr) %*% b),
       lp_test  = as.numeric(as.matrix(Xte) %*% b),
       beta = b, fit = cf, null_clin = FALSE)
}

#inner tuning for the offset arm
tune_offset_rna <- function(clinical_subset, clinical_train, inner_rna, foldid,
                            nlambda = OFFSET_NLAMBDA,
                            lambda_min_ratio = OFFSET_LAMBDA_MIN_RATIO) {
  if (!requireNamespace("glmnet", quietly = TRUE)) stop("glmnet is required.")
  n <- nrow(clinical_train)
  if (length(foldid) != n || anyNA(foldid) || any(foldid != as.integer(foldid)))
    stop("Inner fold IDs do not match clinical rows.")
  ids <- sort(unique(foldid))
  if (length(ids) < 3L || !identical(as.integer(ids), seq_along(ids)) ||
      length(inner_rna) != length(ids)) stop("Invalid inner RNA cache or fold IDs.")
  if (!all(c("OS.time", "OS", clinical_subset) %in% names(clinical_train)))
    stop("Inner clinical data are missing required columns.")
  y <- survival::Surv(clinical_train$OS.time, clinical_train$OS)
  fractions <- exp(seq(0, log(lambda_min_ratio), length.out = nlambda))
  loss <- matrix(NA_real_, length(ids), nlambda)
  valid_n <- integer(length(ids)); inner_jerr <- integer(length(ids))
  inner_clin_cols <- integer(length(ids))
  inner_dropped   <- vector("list", length(ids))
  for (k in ids) {
    part <- inner_rna[[k]]
    tri <- which(foldid != k); vai <- which(foldid == k)
    if (!identical(part$train_idx, tri) || !identical(part$valid_idx, vai))
      stop("Inner RNA cache row order disagrees with foldid in fold ", k, ".")
    if (!is.null(part$rna$error)) stop("Inner fold ", k, " RNA: ", part$rna$error)
    tryCatch({
      cl <- build_clinical(clinical_subset,
                           clinical_train[tri, , drop = FALSE],
                           clinical_train[vai, , drop = FALSE])
      inner_clin_cols[k] <- ncol(cl$train)
      inner_dropped[[k]] <- cl$dropped
      #offset from THIS split's clinical fit only
      off <- .clin_offset(cl$train, cl$test,
                          clinical_train$OS.time[tri], clinical_train$OS[tri])
      Rtr <- as.matrix(part$rna$train); Rva <- as.matrix(part$rna$test)
      pen <- rep(1, ncol(Rtr))  #no unpenalised block here
      use_off <- if (off$null_clin) NULL else off$lp_train
      fit <- .fit_cox_path(Rtr, y[tri, , drop = FALSE], pen,
                           nlambda, lambda_min_ratio, offset = use_off)
      np <- length(fit$lambda)
      all_x <- rbind(Rtr, Rva)
      all_o <- c(off$lp_train, off$lp_test)
      eta <- if (off$null_clin)
        as.matrix(predict(fit, newx = all_x, type = "link"))
      else
        as.matrix(predict(fit, newx = all_x, newoffset = all_o, type = "link"))
      if (!identical(dim(eta), c(n, np)) || any(!is.finite(eta)))
        stop("Non-finite or incorrectly sized inner predictions.")
      d_all <- glmnet::coxnet.deviance(pred = eta, y = y[c(tri, vai), , drop = FALSE],
                                       std.weights = FALSE, cox.ties = COX_TIES)
      d_train <- glmnet::coxnet.deviance(pred = eta[seq_along(tri), , drop = FALSE],
                                         y = y[tri, , drop = FALSE],
                                         std.weights = FALSE, cox.ties = COX_TIES)
      loss[k, seq_len(np)] <- (d_all - d_train) / length(vai)
      valid_n[k] <- length(vai)
      inner_jerr[k] <- if (is.null(fit$jerr)) 0L else fit$jerr
    }, error = function(e) stop("Inner fold ", k, ": ", conditionMessage(e), call. = FALSE))
  }
  eligible <- which(colSums(is.finite(loss)) == length(ids))
  if (!length(eligible)) stop("No relative lambda has a finite score in every inner fold.")
  cvm <- rep(NA_real_, nlambda)
  cvm[eligible] <- colSums(loss[, eligible, drop = FALSE] * valid_n) / sum(valid_n)
  selected <- eligible[which.min(cvm[eligible])]
  dropped_any <- sort(unique(unlist(inner_dropped)))
  list(sel_idx = selected, fraction = fractions[selected], fractions = fractions,
       cvm = cvm, n_valid_path = length(eligible),
       at_boundary = selected >= max(eligible) - 1L, inner_jerr = inner_jerr,
       inner_clin_cols = inner_clin_cols,
       inner_vars_dropped = if (length(dropped_any))
         paste(dropped_any, collapse = "+") else NA_character_)
}

#offset arm : inner tune, outer fit, frozen combined model
fit_offset_rna <- function(Xtr_clin, Xte_clin, Rtr, Rte, time_tr, status_tr,
                           pred_times, foldid,
                           nlambda = OFFSET_NLAMBDA,
                           lambda_min_ratio = OFFSET_LAMBDA_MIN_RATIO,
                           clinical_subset, clinical_train, inner_rna) {
  t0 <- .now()
  tryCatch({
    if (missing(clinical_subset) || missing(clinical_train) || missing(inner_rna))
      stop("Nested tuning requires clinical_subset, clinical_train and inner_rna.")
    if (!identical(as.numeric(clinical_train$OS.time), as.numeric(time_tr)) ||
        !identical(as.numeric(clinical_train$OS), as.numeric(status_tr)))
      stop("Final training outcomes disagree with inner clinical row order.")
    tuning <- tune_offset_rna(clinical_subset, clinical_train, inner_rna, foldid,
                              nlambda, lambda_min_ratio)
    Rtr <- as.matrix(Rtr); Rte <- as.matrix(Rte)
    ytr <- survival::Surv(time_tr, status_tr)
    off <- .clin_offset(Xtr_clin, Xte_clin, time_tr, status_tr)
    pen <- rep(1, ncol(Rtr))
    use_off <- if (off$null_clin) NULL else off$lp_train
    fit <- .fit_cox_path(Rtr, ytr, pen, nlambda, lambda_min_ratio, offset = use_off)
    n_path  <- length(fit$lambda)
    sel_idx <- tuning$sel_idx
    if (sel_idx > n_path)
      stop("Final outer-training path did not reach the selected relative lambda.")
    lam  <- fit$lambda[sel_idx]
    jerr <- fit$jerr
    gb   <- as.numeric(coef(fit, s = lam))
    keep <- which(gb != 0)
    nmR  <- colnames(Rtr); if (is.null(nmR)) nmR <- paste0("g", seq_len(ncol(Rtr)))
    genes <- nmR[keep]
    
    #assemble the FROZEN combined model (Volkmann build_models)
    cn_clin <- colnames(Xtr_clin)
    if (is.null(cn_clin) && ncol(Xtr_clin) > 0)
      cn_clin <- paste0("c", seq_len(ncol(Xtr_clin)))
    dtr <- data.frame(.t = time_tr, .s = status_tr)
    dte <- data.frame(row.names = seq_len(nrow(Rte)))
    init <- numeric(0)
    if (ncol(Xtr_clin) > 0) {
      dtr <- cbind(dtr, setNames(as.data.frame(Xtr_clin), make.names(cn_clin, unique = TRUE)))
      dte <- cbind(dte, setNames(as.data.frame(Xte_clin), make.names(cn_clin, unique = TRUE)))
      init <- c(init, as.numeric(off$beta))
    }
    if (length(keep)) {
      gn <- make.names(paste0("G_", genes), unique = TRUE)
      dtr <- cbind(dtr, setNames(as.data.frame(Rtr[, keep, drop = FALSE]), gn))
      dte <- cbind(dte, setNames(as.data.frame(Rte[, keep, drop = FALSE]), gn))
      init <- c(init, gb[keep])
    }
    if (length(init) == 0L) {
      #Volkmann's Null case: no clinical, no genes -> Kaplan-Meier reference.
      train_secs <- .now() - t0
      km <- fit_km(time_tr, status_tr, nrow(Rte), pred_times)
      risk <- km$risk; S <- km$surv
    } else {
      md_comb <- coxph(Surv(.t, .s) ~ ., data = dtr, init = init,
                       iter.max = 0, ties = COX_TIES, model = TRUE)
      train_secs <- .now() - t0  #tuning + clinical fit + RNA path + frozen model
      risk <- as.numeric(predict(md_comb, newdata = dte, type = "lp"))
      sf   <- survfit(md_comb, newdata = dte)
      S    <- t(summary(sf, times = pred_times, extend = TRUE)$surv)
      if (is.null(dim(S))) S <- matrix(S, nrow = nrow(Rte))
    }
    list(risk = risk, surv = S, times = pred_times, lambda = lam,
         train_secs = round(train_secs, 4),
         n_path = n_path, sel_idx = sel_idx, jerr = jerr,
         at_boundary = sel_idx >= n_path - 1L,
         n_genes = length(keep), genes = genes,
         lambda_fraction = tuning$fraction, inner_n_path = tuning$n_valid_path,
         inner_at_boundary = tuning$at_boundary,
         inner_solver_error = any(tuning$inner_jerr != 0),
         clin_cols = ncol(Xtr_clin),
         inner_clin_mismatch = any(tuning$inner_clin_cols != ncol(Xtr_clin)),
         inner_vars_dropped = tuning$inner_vars_dropped,
         error = NA_character_)
  }, error = function(e)
    list(risk = NA, surv = NA, times = pred_times, lambda = NA, train_secs = NA_real_,
         n_path = NA, sel_idx = NA, jerr = NA, at_boundary = NA,
         n_genes = NA, genes = character(0), lambda_fraction = NA,
         inner_n_path = NA, inner_at_boundary = NA, inner_solver_error = NA,
         clin_cols = NA_integer_, inner_clin_mismatch = NA,
         inner_vars_dropped = NA_character_,
         error = conditionMessage(e)))
}

#RNA-score arm: pre-validated ridge RNA score added to the clinical Cox (Tibshirani and Efron 2002)
SCORE_ALPHA            <- 0
SCORE_LAMBDA_MIN_RATIO <- if (SCORE_ALPHA == 1) LAMBDA_MIN_RATIO else RIDGE_LAMBDA_MIN_RATIO

#tune the RNA-only path and keep each inner model's predictions for its validation rows
.tune_rna_score <- function(y, inner_rna, foldid, alpha, nlambda, lambda_min_ratio) {
  n   <- length(foldid); ids <- sort(unique(foldid))
  if (length(ids) < 3L || length(inner_rna) != length(ids))
    stop("Invalid inner RNA cache or fold IDs.")
  loss  <- matrix(NA_real_, length(ids), nlambda)
  preval <- matrix(NA_real_, n, nlambda)
  valid_n <- integer(length(ids)); inner_jerr <- integer(length(ids))
  for (k in ids) {
    part <- inner_rna[[k]]
    tri <- which(foldid != k); vai <- which(foldid == k)
    if (!identical(part$train_idx, tri) || !identical(part$valid_idx, vai))
      stop("Inner RNA cache row order disagrees with foldid in fold ", k, ".")
    if (!is.null(part$rna$error)) stop("Inner fold ", k, " RNA: ", part$rna$error)
    Rtr <- as.matrix(part$rna$train); Rva <- as.matrix(part$rna$test)
    fit <- .fit_cox_path(Rtr, y[tri, , drop = FALSE], rep(1, ncol(Rtr)),
                         nlambda, lambda_min_ratio, alpha = alpha)
    np  <- length(fit$lambda)
    eta <- as.matrix(predict(fit, newx = rbind(Rtr, Rva), type = "link"))
    if (!identical(dim(eta), c(n, np)) || any(!is.finite(eta)))
      stop("Non-finite or incorrectly sized inner predictions in fold ", k, ".")
    d_all   <- glmnet::coxnet.deviance(pred = eta, y = y[c(tri, vai), , drop = FALSE],
                                       std.weights = FALSE, cox.ties = COX_TIES)
    d_train <- glmnet::coxnet.deviance(pred = eta[seq_along(tri), , drop = FALSE],
                                       y = y[tri, , drop = FALSE],
                                       std.weights = FALSE, cox.ties = COX_TIES)
    loss[k, seq_len(np)] <- (d_all - d_train) / length(vai)
    preval[vai, seq_len(np)] <- eta[length(tri) + seq_along(vai), , drop = FALSE]
    valid_n[k] <- length(vai)
    inner_jerr[k] <- if (is.null(fit$jerr)) 0L else fit$jerr
  }
  eligible <- which(colSums(is.finite(loss)) == length(ids))
  if (!length(eligible)) stop("No relative lambda has a finite score in every inner fold.")
  cvm <- rep(NA_real_, nlambda)
  cvm[eligible] <- colSums(loss[, eligible, drop = FALSE] * valid_n) / sum(valid_n)
  selected <- eligible[which.min(cvm[eligible])]  #first minimum: ties pick the stronger penalty
  list(sel_idx = selected, preval = preval[, selected], n_valid_path = length(eligible),
       at_boundary = selected >= max(eligible) - 1L, inner_jerr = inner_jerr)
}

fit_rna_score <- function(Xtr_clin, Xte_clin, Rtr, Rte, time_tr, status_tr, pred_times,
                          foldid, inner_rna, alpha = SCORE_ALPHA,
                          nlambda = NLAMBDA, lambda_min_ratio = SCORE_LAMBDA_MIN_RATIO) {
  t0 <- .now()
  tryCatch({
    if (!requireNamespace("glmnet", quietly = TRUE)) stop("glmnet is required.")
    if (missing(inner_rna) || is.null(inner_rna)) stop("The score arm requires the inner RNA cache.")
    ytr <- survival::Surv(time_tr, status_tr)
    tun <- .tune_rna_score(ytr, inner_rna, foldid, alpha, nlambda, lambda_min_ratio)
    #the score model refit on the whole training fold, used for the test patients
    Rtr_m <- as.matrix(Rtr); Rte_m <- as.matrix(Rte)
    fit <- .fit_cox_path(Rtr_m, ytr, rep(1, ncol(Rtr_m)), nlambda, lambda_min_ratio, alpha = alpha)
    if (tun$sel_idx > length(fit$lambda))
      stop("Final RNA path did not reach the selected relative lambda.")
    lam <- fit$lambda[tun$sel_idx]
    s_tr <- tun$preval
    s_te <- as.numeric(predict(fit, newx = Rte_m, s = lam, type = "link"))
    n_genes <- sum(as.numeric(coef(fit, s = lam)) != 0)
    #a constant score carries nothing: fall back to the clinical columns alone
    score_used <- n_genes > 0 && is.finite(stats::sd(s_tr)) && stats::sd(s_tr) > 0
    Xtr2 <- as.data.frame(as.matrix(Xtr_clin)); Xte2 <- as.data.frame(as.matrix(Xte_clin))
    if (score_used) { Xtr2$.rna_score <- s_tr; Xte2$.rna_score <- s_te }
    if (ncol(Xtr2) == 0L) {
      m <- fit_km(time_tr, status_tr, nrow(Xte2), pred_times)
      m$separated <- NA
    } else {
      m <- fit_clinical(Xtr2, time_tr, status_tr, Xte2, pred_times)
    }
    if (!is.null(m$error) && !is.na(m$error)) stop("clinical + score Cox: ", m$error)
    m$train_secs    <- round(.now() - t0, 4)  #score tuning + score refit + Cox fit
    m$score_used    <- score_used
    m$score_n_genes <- n_genes
    m$sel_idx       <- tun$sel_idx
    m$at_boundary   <- tun$at_boundary
    m$error         <- NA_character_
    m
  }, error = function(e)
    list(risk = NA, surv = NA, times = pred_times, train_secs = NA_real_,
         score_used = NA, score_n_genes = NA_integer_, separated = NA,
         error = conditionMessage(e)))
}

#forest arm: blockForest (clinical + RNA) against a clinical-only RSF

FOREST_TREES     <- c(250L, 500L, 1000L, 2000L)  #all always reported; max = blockfor() default
BF_NUM_TREES_PRE <- 1500L  #PACKAGE DEFAULT (blockfor): trees per trial forest
BF_NSETS         <- 300L  #PACKAGE DEFAULT (blockfor): random candidates, used when BF_WEIGHT_GRID is NULL
BF_WEIGHT_GRID   <- seq(0.05, 1, by = 0.05)  #RNA weight candidates with clinical fixed at 1 (tree test, 28 Sep)
RF_SPLITRULE     <- "extratrees"  #PACKAGE DEFAULT (blockfor), also used for the clinical RSF
RF_THREADS       <- 1L  #single-threaded: cells are already parallel in cv.R
BF_FIX_CLINICAL_WEIGHT <- TRUE  #TRUE: clinical weight fixed at 1, only the RNA weight tuned

#map forest survival onto the IBS grid
.forest_on_grid <- function(S_raw, death_times, pred_times) {
  S_raw <- as.matrix(S_raw)
  idx <- findInterval(pred_times, death_times)
  out <- matrix(1, nrow = nrow(S_raw), ncol = length(pred_times))
  has <- idx > 0
  out[, has] <- S_raw[, idx[has], drop = FALSE]
  out
}

#predictions at every tree count from ONE grown forest
.forest_preds <- function(forest, newdata, pred_times, trees = FOREST_TREES) {
  out <- lapply(trees, function(nt) {
    pr <- predict(forest, data = newdata, num.trees = nt, num.threads = RF_THREADS)
    if (!is.null(pr$num.trees) && pr$num.trees != nt)
      stop(sprintf("predict() used %s trees when %d were requested.", pr$num.trees, nt))
    chf <- as.matrix(pr$chf)
    list(risk = rowSums(chf),
         surv = .forest_on_grid(pr$survival, pr$unique.death.times, pred_times))
  })
  names(out) <- as.character(trees)
  out
}

#plain random survival forest (ranger): the clinical comparator (mtry = all columns)
fit_rsf <- function(Xtr, time_tr, status_tr, Xte, pred_times, seed, mtry = NULL) {
  warn <- character(0); t0 <- .now()
  tryCatch(
    withCallingHandlers({
      if (!requireNamespace("ranger", quietly = TRUE)) stop("ranger is required.")
      dtr <- data.frame(as.matrix(Xtr), check.names = TRUE)
      dte <- data.frame(as.matrix(Xte), check.names = TRUE)
      dtr$.t <- time_tr; dtr$.s <- status_tr
      fit <- ranger::ranger(data = dtr,
                            dependent.variable.name = ".t",
                            status.variable.name    = ".s",
                            num.trees   = max(FOREST_TREES),
                            mtry        = mtry,
                            splitrule   = RF_SPLITRULE,
                            num.threads = RF_THREADS,
                            seed        = seed)
      train_secs <- .now() - t0
      list(by_trees = .forest_preds(fit, dte, pred_times), times = pred_times,
           train_secs = round(train_secs, 4), mtry = fit$mtry, w_clin = NA_real_, w_rna = NA_real_,
           warning = if (length(warn)) paste(unique(warn), collapse = " | ") else NA_character_,
           error = NA_character_)
    },
    warning = function(w) {
      warn <<- c(warn, conditionMessage(w))
      invokeRestart("muffleWarning")
    }),
    error = function(e)
      list(by_trees = NULL, times = pred_times, train_secs = NA_real_, mtry = NA_integer_,
           w_clin = NA_real_, w_rna = NA_real_,
           warning = NA_character_, error = conditionMessage(e)))
}

#blockfor() block-weight search with the clinical weight fixed at 1
.blockfor_fixed_clinical <- function(Xtr_m, time_tr, status_tr, blocks) {
  df <- data.frame(time = time_tr, status = status_tr, Xtr_m, check.names = FALSE)
  mtry <- sapply(blocks, function(x) sqrt(length(x)))
  mtry[1] <- length(blocks[[1]])  #clinical: all variables, as always.select.block
  cand <- if (is.null(BF_WEIGHT_GRID)) runif(BF_NSETS) else BF_WEIGHT_GRID
  errs <- numeric(length(cand)); sets <- vector("list", length(cand))
  for (l in seq_along(cand)) {
    w <- c(1, cand[l])
    sets[[l]] <- w
    fo <- blockForest::blockForest(dependent.variable.name = "time",
                                   status.variable.name = "status", data = df,
                                   num.trees = BF_NUM_TREES_PRE, blocks = blocks,
                                   block.weights = w, mtry = as.integer(mtry),
                                   keep.inbag = TRUE, block.method = "BlockForest",
                                   splitrule = RF_SPLITRULE, write.forest = FALSE,
                                   num.threads = RF_THREADS)
    errs[l] <- fo$prediction.error
  }
  best <- which(errs == min(errs))
  best <- if (length(best) > 1) sample(best, 1) else best  #random tie-break, as the package
  w <- sets[[best]]
  forest <- blockForest::blockForest(dependent.variable.name = "time",
                                     status.variable.name = "status", data = df,
                                     num.trees = max(FOREST_TREES), blocks = blocks,
                                     block.weights = w, mtry = mtry, keep.inbag = TRUE,
                                     block.method = "BlockForest", splitrule = RF_SPLITRULE,
                                     num.threads = RF_THREADS)
  list(forest = forest, paramvalues = w)
}

#blockForest on clinical + RNA
fit_blockforest <- function(Xtr, time_tr, status_tr, Xte, block, pred_times, seed) {
  if (!any(block == "clinical"))  #anchor: RNA only
    return(fit_rsf(Xtr, time_tr, status_tr, Xte, pred_times, seed, mtry = NULL))
  warn <- character(0); t0 <- .now()
  tryCatch(
    withCallingHandlers({
      if (!requireNamespace("blockForest", quietly = TRUE))
        stop("blockForest is required.")
      Xtr_m <- as.matrix(Xtr); Xte_m <- as.matrix(Xte)
      if (ncol(Xtr_m) != length(block) || ncol(Xte_m) != length(block))
        stop("Block labels do not match the fused design.")
      safe <- paste0("V", seq_len(ncol(Xtr_m)))
      colnames(Xtr_m) <- safe; colnames(Xte_m) <- safe
      blocks <- list(clinical = which(block == "clinical"),
                     rna      = which(block == "rna"))
      set.seed(seed)  #candidate weight sets (and ranger seeds) come from R's RNG
      if (isTRUE(BF_FIX_CLINICAL_WEIGHT)) {
        bf <- .blockfor_fixed_clinical(Xtr_m, time_tr, status_tr, blocks)
      } else {
        bf <- blockForest::blockfor(Xtr_m, cbind(time_tr, status_tr),
                                    blocks = blocks, block.method = "BlockForest",
                                    num.trees = max(FOREST_TREES),
                                    mtry = NULL,  #package default
                                    nsets = BF_NSETS,
                                    num.trees.pre = BF_NUM_TREES_PRE,
                                    splitrule = RF_SPLITRULE,
                                    always.select.block = 1L,  #= blocks$clinical
                                    num.threads = RF_THREADS)
      }
      train_secs <- .now() - t0  #block-weight tuning + final forest
      pv <- as.numeric(bf$paramvalues)
      list(by_trees = .forest_preds(bf$forest, Xte_m, pred_times), times = pred_times,
           train_secs = round(train_secs, 4), mtry = NA_integer_,
           w_clin = if (length(pv) >= 1) pv[1] else NA_real_,
           w_rna  = if (length(pv) >= 2) pv[2] else NA_real_,
           warning = if (length(warn)) paste(unique(warn), collapse = " | ") else NA_character_,
           error = NA_character_)
    },
    warning = function(w) {
      warn <<- c(warn, conditionMessage(w))
      invokeRestart("muffleWarning")
    }),
    error = function(e)
      list(by_trees = NULL, times = pred_times, train_secs = NA_real_, mtry = NA_integer_,
           w_clin = NA_real_, w_rna = NA_real_,
           warning = NA_character_, error = conditionMessage(e)))
}


#boosting arm: CoxBoost favoring clinical, steps tuned on the fused arm's inner folds
COXBOOST_MAXSTEPNO <- 100L  #package default

fit_coxboost <- function(Xtr, time_tr, status_tr, Xte, block, pred_times, inner_foldid) {
  warn <- character(0); t0 <- .now()
  tryCatch(
    withCallingHandlers({
      if (!requireNamespace("CoxBoost", quietly = TRUE)) stop("CoxBoost is required.")
      Xtr_m <- as.matrix(Xtr)
      Xte_m <- as.matrix(Xte)
      if (ncol(Xtr_m) != length(block)) stop("Block labels do not match the fused design.")
      unpen <- which(block == "clinical")
      if (!length(unpen)) unpen <- NULL  #anchor: RNA only
      M <- COXBOOST_MAXSTEPNO
      folds <- sort(unique(inner_foldid))
      contrib <- vapply(folds, function(k) {
        tr <- inner_foldid != k
        fk <- CoxBoost::CoxBoost(time_tr[tr], status_tr[tr], Xtr_m[tr, , drop = FALSE],
                                 unpen.index = unpen, stepno = M, return.score = FALSE)
        full <- predict(fk, newdata = Xtr_m, newtime = time_tr, newstatus = status_tr,
                        type = "logplik", at.step = 0:M)
        part <- predict(fk, newdata = Xtr_m[tr, , drop = FALSE], newtime = time_tr[tr],
                        newstatus = status_tr[tr], type = "logplik", at.step = 0:M)
        as.numeric(full - part)
      }, numeric(M + 1))
      cvm  <- rowSums(contrib)
      if (all(!is.finite(cvm))) stop("CoxBoost CV criterion not finite at any step.")
      step <- which.max(replace(cvm, !is.finite(cvm), -Inf)) - 1L  #0..M
      fit  <- CoxBoost::CoxBoost(time_tr, status_tr, Xtr_m, unpen.index = unpen,
                                 stepno = step, return.score = FALSE)
      train_secs <- .now() - t0  #step tuning on inner folds + final fit
      risk <- as.numeric(predict(fit, newdata = Xte_m, type = "lp"))
      surv <- predict(fit, newdata = Xte_m, times = pred_times, type = "risk")
      b    <- as.numeric(coef(fit))
      pen_cols <- setdiff(seq_len(ncol(Xtr_m)), unpen)
      list(risk = risk, surv = as.matrix(surv), times = pred_times,
           train_secs = round(train_secs, 4), stepno = step, at_boundary = step >= M,
           n_genes = sum(b[pen_cols] != 0),
           warning = if (length(warn)) paste(unique(warn), collapse = " | ") else NA_character_,
           error = NA_character_)
    },
    warning = function(w) {
      warn <<- c(warn, conditionMessage(w))
      invokeRestart("muffleWarning")
    }),
    error = function(e)
      list(risk = NA_real_, surv = NULL, times = pred_times, train_secs = NA_real_,
           stepno = NA_integer_,
           at_boundary = NA, n_genes = NA_integer_, warning = NA_character_,
           error = conditionMessage(e)))
}