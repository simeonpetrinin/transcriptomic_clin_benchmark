#Metrics
#Uno's IPCW C-index (primary), the integrated Brier score, SurvBoard's D-calibration,
#the ICI (calibration at one time) and the partial-likelihood ratio per death.
#IPCW weights use G(T-) with a floor (G_FLOOR), counted when it binds.

suppressPackageStartupMessages(library(survival))

G_FLOOR <- 0.05  #censoring survival is not allowed below this in the weights

#tau = maximum event (death) time
make_tau <- function(time, status) max(time[status == 1], na.rm = TRUE)

#Uno's IPCW C-index
uno_c <- function(time, status, risk, tau = NULL) {
  if (is.null(tau)) tau <- make_tau(time, status)
  cc <- concordance(Surv(time, status) ~ risk,
                    reverse = TRUE, timewt = "n/G2", ymax = tau)
  as.numeric(cc$concordance)
}

#censoring survival: LEFT LIMIT and value AT t
make_G <- function(cens_time, cens_status, floor_at = G_FLOOR) {
  cg <- survfit(Surv(cens_time, 1 - cens_status) ~ 1)
  tt <- c(0, cg$time)  #step function: G is right-continuous
  ss <- c(1, cg$surv)
  n_floored <- 0L
  
  floor_G <- function(g) {
    hit <- g < floor_at
    if (any(hit)) n_floored <<- n_floored + sum(hit)
    pmax(g, floor_at)
  }
  
  G_left <- function(t) {
    #LEFT limit: the value strictly BEFORE t, i.e. the last step at time < t
    i <- findInterval(t - .Machine$double.eps^0.5, tt)
    floor_G(ss[pmax(i, 1L)])
  }
  
  G_at <- function(t) {
    #RIGHT-continuous KM value AT t, including any censoring drop at exactly t
    i <- findInterval(t, tt)
    floor_G(ss[pmax(i, 1L)])
  }
  
  list(G_left = G_left, G_at = G_at, floored = function() n_floored)
}

#Integrated Brier Score (IPCW, Graf et al.)
integrated_brier <- function(time, status, surv_prob, pred_times, tau = NULL,
                             cens_time, cens_status) {
  if (is.null(tau)) tau <- make_tau(time, status)
  keep       <- pred_times <= tau
  pred_times <- pred_times[keep]
  surv_prob  <- surv_prob[, keep, drop = FALSE]
  if (!length(pred_times)) return(NA_real_)
  
  gg <- make_G(cens_time, cens_status)
  
  bs <- numeric(length(pred_times))
  for (k in seq_along(pred_times)) {
    t   <- pred_times[k]; Sk <- surv_prob[, k]
    Gt  <- gg$G_at(t)  #G(t) for the still-at-risk term
    GTi <- gg$G_left(time)  #G(T_i-) for the died-before-t term
    w1  <- (time <= t & status == 1) / GTi  #died on/before t
    w2  <- (time >  t) / Gt  #still at risk after t
    bs[k] <- mean(w1 * Sk^2 + w2 * (1 - Sk)^2)
  }
  out <- if (length(pred_times) < 2) mean(bs) else
    sum(diff(pred_times) * (head(bs, -1) + tail(bs, -1)) / 2) /
    (max(pred_times) - min(pred_times))  #trapezoid, normalised by span
  attr(out, "n_floored") <- gg$floored()
  out
}

#D-calibration, implemented as in SurvBoard (Wissel et al. 2025
D_CAL_BINS <- 5

d_calibration <- function(time, status, surv_prob, pred_times, n_bins = D_CAL_BINS) {
  surv_prob <- as.matrix(surv_prob)
  if (nrow(surv_prob) != length(time) || ncol(surv_prob) != length(pred_times))
    stop("surv_prob must be patients x pred_times.")
  if (anyNA(time) || anyNA(status) || any(!is.finite(surv_prob))) return(NA_real_)
  
  #1. predicted survival at each patient's own time (step read, clipped)
  j <- pmin(pmax(findInterval(time, pred_times), 1L), length(pred_times))
  p <- pmin(pmax(surv_prob[cbind(seq_along(time), j)], 0), 1)
  
  #2. bin: an event counts 1 in its bin
  bin    <- pmin(floor(p * n_bins), n_bins - 1) + 1  #bin b covers [(b-1)/B, b/B)
  event  <- status == 1
  counts <- tabulate(bin[event], nbins = n_bins)
  for (k in which(!event)) {
    b <- bin[k]
    if (p[k] <= 0) { counts[1] <- counts[1] + 1; next }  #SurvBoard divides by 0 here
    counts[b] <- counts[b] + 1 - ((b - 1) / n_bins) / p[k]
    if (b > 1) counts[1:(b - 1)] <- counts[1:(b - 1)] + 1 / (n_bins * p[k])
  }
  
  #3. chi-square against a flat line
  expected <- length(p) / n_bins
  stat <- sum((counts - expected)^2 / expected)
  attr(stat, "p_value") <- pchisq(stat, df = n_bins - 1, lower.tail = FALSE)
  attr(stat, "bins")    <- counts / length(p)
  stat
}

#knots as rms/Hmisc defaults for 3 knots: 10/50/90th percentiles; outer knots 5th from each end if n < 100
.ici_knots <- function(x) {
  k <- unname(quantile(x, c(0.1, 0.5, 0.9)))
  if (length(x) < 100) { xs <- sort(x); k[1] <- xs[5]; k[3] <- xs[length(x) - 4] }
  k
}

#ICI at time t0 (Austin, Harrell & van Klaveren 2020): restricted cubic spline, 3 knots
ici <- function(time, status, surv_t0, t0) {
  tryCatch({
    s <- pmin(pmax(surv_t0, 1e-6), 1 - 1e-6)
    x <- log(-log(s)); p <- 1 - s
    if (sum(status) < 2 || length(unique(x)) < 2) return(NA_real_)
    kn <- if (length(x) >= 10 && length(unique(x)) > 5) .ici_knots(x) else NA
    d  <- data.frame(t = time, e = status, x = x)
    f  <- if (length(unique(kn)) == 3 && !anyNA(kn))
      coxph(Surv(t, e) ~ splines::ns(x, knots = kn[2], Boundary.knots = kn[c(1, 3)]), data = d)
    else coxph(Surv(t, e) ~ x, data = d)  #too few distinct predictions for a spline
    obs <- 1 - summary(survfit(f, newdata = d, se.fit = FALSE), times = t0, extend = TRUE)$surv
    v <- mean(abs(p - as.numeric(obs)))
    if (is.finite(v)) v else NA_real_
  }, error = function(e) NA_real_)
}

#Cox deviance per patient of risk scores on the same patients (glmnet, same ties as the fits)
cox_deviance <- function(time, status, risk) {
  if (length(risk) != length(time) || any(!is.finite(risk))) return(NA_real_)
  glmnet::coxnet.deviance(pred = matrix(risk, ncol = 1), y = Surv(time, status),
                          std.weights = FALSE, cox.ties = COX_TIES) / length(time)
}

#probability ratio per death vs no predictors (Bovelstad et al. 2009 deviance difference, rescaled)
pl_ratio <- function(pl_null, pl_model, n, deaths) exp(n * (pl_null - pl_model) / (2 * deaths))