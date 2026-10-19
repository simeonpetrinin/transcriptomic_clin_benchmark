#Robustness sample check
#After the main lattice: how well do age alone, the full set and n sampled subsets recover each cancer's full-lattice line of fused delta C, against quantity (number of variables) and quality (clinical C)?
#Redraws the sample many times and saves spreads of the error, never the lines themselves.
#The sampling error is set against the full lattice's own between-repeat error (delete-one-repeat jackknife), as Monte Carlo error is judged in simulation studies (Koehler et al. 2009; Morris et al. 2019).
#Each robustness learner is then scaled from fused by its scatter and repeat noise on the robustness subsets (an assumption, not a measurement).
#Also checks the median, P10 and P90 of delta C across subsets. The sample is spread evenly over subset sizes,
#so these are weighted: each sampled subset counts for the lattice subsets of its size it stands in for (N_k / n_k).

RSC_N      <- c(10L, 20L, 40L, 80L, 160L)  #sample sizes to compare
RSC_DRAWS  <- 1000L  #redraws per cancer and sample size
RSC_SEED0  <- 7700000L  #seeds of the redraws (the real sample keeps ROBUST_SEED)

#Theil-Sen slope: median of the slopes between all pairs with different x (Sen 1968)
.rsc_theil_sen <- function(x, y) {
  s <- unlist(lapply(seq_len(length(x) - 1), function(i) {  #row by row, to spare memory
    j <- (i + 1):length(x); dx <- x[j] - x[i]
    ((y[j] - y[i]) / dx)[dx != 0]
  }))
  median(s)
}

#slope of y on x by least squares and by Theil-Sen, NA if too few points
.rsc_fit <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  if (sum(ok) < 3 || length(unique(x[ok])) < 2) return(c(ls = NA_real_, ts = NA_real_))
  c(ls = unname(coef(lm(y[ok] ~ x[ok]))[2]), ts = .rsc_theil_sen(x[ok], y[ok]))
}

#the three lines for a set of subsets: quantity, quality (split repeats), quality (Oldham)
.rsc_slopes <- function(vars_in, sm, pr) {
  s <- sm[sm$variables %in% vars_in, ]
  p <- pr[pr$variables %in% vars_in, ]
  half <- function(odd, col) {  #per-subset mean of col over odd or even repeats
    q <- p[(p$repeat_id %% 2 == 1) == odd, ]
    tapply(q[[col]], q$variables, mean, na.rm = TRUE)[s$variables]
  }
  split1 <- .rsc_fit(half(TRUE, "c_ref"), half(FALSE, "delta_c"))  #x odd repeats, y even
  split2 <- .rsc_fit(half(FALSE, "c_ref"), half(TRUE, "delta_c"))  #and the swap
  f <- list(quantity       = .rsc_fit(s$n_vars, s$delta_c),
            quality_split  = (split1 + split2) / 2,
            quality_oldham = .rsc_fit((s$c_ref + s$c_model) / 2, s$delta_c))
  unlist(f)  #named quantity.ls, quantity.ts, quality_split.ls, ...
}

#weighted quantile (inverse of the weighted distribution); equal weights give an ordinary quantile
.rsc_wq <- function(x, w, p) {
  ok <- !is.na(x) & !is.na(w); x <- x[ok]; w <- w[ok]
  if (!length(x)) return(rep(NA_real_, length(p)))
  o <- order(x); x <- x[o]; cw <- cumsum(w[o]) / sum(w)
  vapply(p, function(pp) x[which(cw >= pp - 1e-12)[1]], numeric(1))
}

#median, P10 and P90 of delta C over a set of subsets; w = weights named by subset (NULL = all equal)
.rsc_quants <- function(vars_in, sm, w = NULL) {
  s <- sm[sm$variables %in% vars_in, ]
  ww <- if (is.null(w)) rep(1, nrow(s)) else unname(w[s$variables])
  q <- .rsc_wq(s$delta_c, ww, c(0.5, 0.1, 0.9))
  c(median.wq = q[1], p10.wq = q[2], p90.wq = q[3])
}

#every number checked: the three lines (unweighted) and the three quantiles (weighted)
.rsc_numbers <- function(vars_in, sm, pr, w = NULL) c(.rsc_slopes(vars_in, sm, pr), .rsc_quants(vars_in, sm, w))

#per-subset means over the repeats in pr (repeats are already fold-collapsed, as in 02_summarise.R)
.rsc_means <- function(pr)
  aggregate(cbind(c_ref, c_model, delta_c) ~ variables + n_vars, data = pr,
            FUN = function(v) mean(v, na.rm = TRUE), na.action = na.pass)

#delete-one-repeat jackknife SE of every full-lattice slope: the lattice's own between-repeat error
.rsc_jack_se <- function(pr, labels) {
  reps <- sort(unique(pr$repeat_id)); R <- length(reps)
  th <- t(vapply(reps, function(r) {
    p <- pr[pr$repeat_id != r, ]; s <- .rsc_means(p)
    unname(.rsc_numbers(s$variables, s, p))
  }, numeric(length(labels))))
  colnames(th) <- labels
  sqrt((R - 1) / R * colSums(sweep(th, 2, colMeans(th, na.rm = TRUE))^2, na.rm = TRUE))
}

#each robustness learner against fused on the same subsets: scatter about the line per unit of x, and repeat noise
.rsc_learner_ratios <- function(code) {
  key <- robust_key(code)
  rs <- read.csv(path_summary(key), stringsAsFactors = FALSE); rs <- rs[!rs$is_anchor, ]
  rp <- read.csv(file.path(DIR_RESULTS, paste0(key, "_per_repeat.csv")), stringsAsFactors = FALSE); rp <- rp[!rp$is_anchor, ]
  scat <- function(x, y) if (sum(!is.na(x) & !is.na(y)) < 4) NA_real_ else sd(resid(lm(y ~ x))) / sd(x[!is.na(y)], na.rm = TRUE)
  noise <- function(p) mean(tapply(p$delta_c, p$variables, sd, na.rm = TRUE), na.rm = TRUE)  #mean over subsets of the SD across repeats
  one <- function(arm, keep) {
    s <- rs[rs$arm == arm & rs$variables %in% keep, ]; p <- rp[rp$arm == arm & rp$variables %in% keep, ]
    c(quantity = scat(s$n_vars, s$delta_c), quality = scat((s$c_ref + s$c_model) / 2, s$delta_c),
      spread = sd(s$delta_c, na.rm = TRUE),  #how far subsets differ: drives the quantiles' sampling error
      noise = noise(p))
  }
  ok_vars <- function(arm) rs$variables[rs$arm == arm & !is.na(rs$delta_c)]
  do.call(rbind, lapply(setdiff(unique(rs$arm), "fused"), function(arm) {
    keep <- intersect(ok_vars(arm), ok_vars("fused"))  #matched subsets only
    a <- one(arm, keep); f <- one("fused", keep)
    data.frame(cancer = code, arm = arm, n_subsets = length(keep),
               scatter_ratio_quantity = a[["quantity"]] / f[["quantity"]],
               scatter_ratio_quality = a[["quality"]] / f[["quality"]],
               spread_ratio = a[["spread"]] / f[["spread"]],
               noise_ratio = a[["noise"]] / f[["noise"]])
  }))
}

#one cancer: full-lattice lines against redrawn samples, summarised as spreads only
.rsc_cancer <- function(code, n, jse) {
  sm <- read.csv(path_summary(code), stringsAsFactors = FALSE)
  pr <- read.csv(file.path(DIR_RESULTS, paste0(code, "_per_repeat.csv")), stringsAsFactors = FALSE)
  sm <- sm[sm$arm == "fused" & !sm$is_anchor, ]
  pr <- pr[pr$arm == "fused" & !pr$is_anchor, ]
  full <- .rsc_numbers(sm$variables, sm, pr)
  #how far each x-axis runs across the lattice, to turn a slope error into a delta C error (quantiles are already in delta C)
  span <- c(quantity = length(CANCER_VARS[[code]]) - 1,
            quality_split = diff(range(sm$c_ref, na.rm = TRUE)),
            quality_oldham = diff(range((sm$c_ref + sm$c_model) / 2, na.rm = TRUE)),
            median = 1, p10 = 1, p90 = 1)
  pool <- robust_pool(code)
  N_k <- table(lengths(pool))
  pick <- function(seed) {  #sampled specs without the omics-only anchor, as variable strings, with their weights
    sp <- robust_specs(code, n, seed = seed, pool = pool)[-1]
    v <- vapply(sp, paste, character(1), collapse = "+")
    k <- as.character(lengths(sp)); smp <- grepl("^sample_", names(sp)); n_k <- table(k[smp])
    w <- ifelse(smp, as.numeric(N_k[k] / n_k[k]), 1)  #age alone and the full set stand for themselves
    list(v = v, w = setNames(w, v))
  }
  num <- function(p) .rsc_numbers(p$v, sm, pr, p$w)
  real <- num(pick(ROBUST_SEED + sum(utf8ToInt(code))))
  whole <- length(pool) <= n
  draws <- if (whole) matrix(full, nrow = 1, dimnames = list(NULL, names(full))) else
    t(vapply(RSC_SEED0 + seq_len(RSC_DRAWS), function(seed) num(pick(seed)), numeric(length(full))))
  do.call(rbind, lapply(names(full), function(a) {
    e  <- draws[, a] - full[a]  #each draw's error
    ci <- if (whole) c(0, 0) else unname(quantile(e, c(0.025, 0.975), na.rm = TRUE))
    sdv <- if (whole) 0 else sd(e, na.rm = TRUE)
    rmse <- if (whole) 0 else sqrt(mean(e^2, na.rm = TRUE))  #scatter plus any bias of the weighted estimate
    an <- sub("\\..*$", "", a)
    real_err <- if (n == ROBUST_N_SAMPLED) unname(real[a] - full[a]) else NA_real_  #only n = 10 is the sample actually run
    #the same errors in delta C over the whole x-axis, comparable across the three lines
    data.frame(cancer = code, analysis = an,
               method = c(ls = "least_squares", ts = "theil_sen", wq = "weighted_quantile")[[sub("^.*\\.", "", a)]],
               n_sampled = n, lattice = nrow(sm), whole_lattice = whole,
               full_value = unname(full[a]),
               slope_err_lo95 = ci[1], slope_err_hi95 = ci[2],
               range_err_lo95 = ci[1] * span[[an]], range_err_hi95 = ci[2] * span[[an]],
               real_sample_range_err = real_err * span[[an]],
               sampling_sd = sdv, sampling_rmse = rmse, full_lattice_se = unname(jse[a]),
               ratio = rmse / unname(jse[a]),  #sampling error per unit of the lattice's own error
               se_inflation = sqrt(1 + (rmse / unname(jse[a]))^2))  #how much the sample widens the total error
  }))
}

#every cancer with main-lattice and robustness summaries; writes results/robust_sample_check.csv and _learners.csv
check_robust_sample <- function(codes = CANCERS, n = RSC_N) {
  main_ok <- file.exists(vapply(codes, path_summary, character(1))) &
    file.exists(file.path(DIR_RESULTS, paste0(codes, "_per_repeat.csv")))
  rob_ok <- file.exists(vapply(robust_key(codes), path_summary, character(1))) &
    file.exists(file.path(DIR_RESULTS, paste0(robust_key(codes), "_per_repeat.csv")))
  if (any(!main_ok)) cat("[sample check] no main summary or per-repeat file, skipped:", paste(codes[!main_ok], collapse = ", "), "\n")
  if (any(!rob_ok)) cat("[sample check] no robustness summary, learners not scaled for:", paste(codes[!rob_ok], collapse = ", "), "\n")
  out <- list(); lrn <- list()
  for (code in codes[main_ok]) {
    pr <- read.csv(file.path(DIR_RESULTS, paste0(code, "_per_repeat.csv")), stringsAsFactors = FALSE)
    pr <- pr[pr$arm == "fused" & !pr$is_anchor, ]
    labels <- names(.rsc_numbers(unique(pr$variables), .rsc_means(pr), pr))
    jse <- .rsc_jack_se(pr, labels)
    fu <- do.call(rbind, lapply(n, function(k) .rsc_cancer(code, k, jse)))
    out[[code]] <- fu
    cat(sprintf("[sample check] %-5s done (%d repeats in the jackknife)\n", code, length(unique(pr$repeat_id))))
    if (!rob_ok[match(code, codes)]) next
    lr <- .rsc_learner_ratios(code)
    #scale fused's sampling error by each learner's scatter, and the lattice's own error by its repeat noise
    lrn[[code]] <- do.call(rbind, lapply(seq_len(nrow(lr)), function(i) {
      r <- lr[i, ]
      sc <- ifelse(fu$analysis == "quantity", r$scatter_ratio_quantity,
                   ifelse(fu$method == "weighted_quantile", r$spread_ratio, r$scatter_ratio_quality))
      rt <- fu$ratio * sc / r$noise_ratio
      data.frame(fu[, c("cancer", "analysis", "method", "n_sampled", "whole_lattice")], arm = r$arm,
                 n_matched_subsets = r$n_subsets, scatter_ratio = sc, noise_ratio = r$noise_ratio,
                 ratio = rt, se_inflation = sqrt(1 + rt^2))
    }))
  }
  out <- do.call(rbind, out); lrn <- do.call(rbind, lrn)
  f1 <- file.path(DIR_RESULTS, "robust_sample_check.csv"); f2 <- file.path(DIR_RESULTS, "robust_sample_check_learners.csv")
  write.csv(out, f1, row.names = FALSE)
  if (!is.null(lrn)) write.csv(lrn, f2, row.names = FALSE)
  print(out[out$method %in% c("least_squares", "weighted_quantile"),
            c("cancer", "analysis", "n_sampled", "whole_lattice", "range_err_lo95", "range_err_hi95", "ratio", "se_inflation")],
        row.names = FALSE, digits = 3)
  cat("\nratio: sampling error (root mean square, includes bias) over the full lattice's jackknife SE. se_inflation = sqrt(1 + ratio^2).\n")
  cat("Written to", f1, if (!is.null(lrn)) paste("and", f2) else "", "\n")
  invisible(list(fused = out, learners = lrn))
}
