# DOMAIN SIMILARITY CHECK  ->  do the proposed clinical domains hold together?
# OUTCOME-BLIND: survival is never read, so nothing here can bias the analysis.
#
# MEASURE: bias-corrected Cramer's V (Bergsma 2013) for every pair of variables
#   in a cohort. 0 = unrelated, 1 = one variable determines the other.
#   Numeric variables are binned into quintiles so every pair uses one measure.
#   Missing values (NA, TCGA "[...]" placeholders, not-assessed codes) are
#   excluded pairwise: each pair uses the patients who have both values.
#
# PRE-SPECIFIED RULE (set before looking at any result):
#   Tested only for CONSTRUCT domains (domains claiming to measure one thing).
#   COLLECTION domains (Demographics, Medical & family history, Location) group
#   the same KIND of information and are not expected to correlate.
#   A construct domain is SUPPORTED if
#     (1) main:      median within-domain V  >  median V between its variables
#                    and all variables outside the domain, and
#     (2) secondary: median within-domain V  >= MIN_V (0.3, conventional
#                    "moderate" association; Cohen 1988).
#   A variable is FLAGGED if its mean V with another construct domain exceeds
#   its mean V with its own domain.
#
# INPUT : domain_assignment.xlsx, sheet "Long (for code)"
#         (Cancer | Domain | Variable (label) | Source column)
# OUTPUT: results/domain_check/<CODE>_pairwise_V.csv
#         results/domain_check/<CODE>_heatmap.png
#         results/domain_check/domain_summary.csv
#         results/domain_check/variable_flags.csv
#
# USAGE (project root):  source("scripts/check_domain_similarity.R")
# NOTE: uses all CDR/PanCanAtlas patients of each cancer type, not only the
# RNA-matched analysis cohort; association structure should be very similar.

source("paths.R")
suppressPackageStartupMessages({ library(readxl); library(data.table) })

ASSIGN_PATH <- "domain_assignment.xlsx"          # adjust if stored elsewhere
OUT_DIR     <- file.path(DIR_RESULTS, "domain_check")
MIN_V       <- 0.3
N_BINS      <- 5
MIN_PAIR_N  <- 30                                # pairs with fewer patients -> NA
COLLECTION_DOMAINS <- c("Demographics", "Medical & family history", "Location")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

#missing-value handling
NOT_ASSESSED <- c("mx", "nx", "tx", "gx", "rx", "n.d.", "nd", "unknown", "indeterminate",
                  "not evaluated", "not available", "discrepancy", "")
clean <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | grepl("^\\[.*\\]$", x) | tolower(x) %in% NOT_ASSESSED] <- NA
  x
}
# numeric if >= 80% of non-missing values parse as numbers
as_analysable <- function(x) {
  x <- clean(x); nn <- x[!is.na(x)]
  num <- suppressWarnings(as.numeric(nn))
  if (length(nn) && mean(!is.na(num)) >= 0.8) {
    v <- suppressWarnings(as.numeric(x))
    br <- unique(quantile(v, probs = seq(0, 1, length.out = N_BINS + 1), na.rm = TRUE))
    if (length(br) < 3) return(factor(v))          # too few distinct values: treat as categories
    return(cut(v, br, include.lowest = TRUE))
  }
  factor(tolower(x))
}

#bias-corrected Cramer's V (Bergsma 2013)
cramers_v_bc <- function(a, b) {
  ok <- !is.na(a) & !is.na(b)
  if (sum(ok) < MIN_PAIR_N) return(NA_real_)
  tab <- table(droplevels(a[ok]), droplevels(b[ok]))
  r <- nrow(tab); k <- ncol(tab); n <- sum(tab)
  if (r < 2 || k < 2) return(NA_real_)
  chi2  <- suppressWarnings(chisq.test(tab, correct = FALSE)$statistic)
  phi2  <- max(0, chi2 / n - (k - 1) * (r - 1) / (n - 1))
  kcorr <- k - (k - 1)^2 / (n - 1); rcorr <- r - (r - 1)^2 / (n - 1)
  den   <- min(kcorr - 1, rcorr - 1)
  if (den <= 0) return(NA_real_)
  unname(sqrt(phi2 / den))
}

#load sources once
assign <- as.data.frame(read_xlsx(ASSIGN_PATH, sheet = "Long (for code)"))
names(assign) <- c("cancer", "domain", "label", "column")
cdr <- as.data.frame(read_xlsx(CDR_PATH, sheet = "TCGA-CDR", guess_max = 11000))
pan_cols <- unique(c("bcr_patient_barcode", "acronym", assign$column))
pan_hdr  <- names(fread(PANCAN_PATH, sep = "\t", nrows = 0))
pan <- fread(PANCAN_PATH, sep = "\t", select = intersect(pan_cols, pan_hdr),
             colClasses = "character", data.table = FALSE, showProgress = FALSE)

get_column <- function(code, col) {
  # CDR first (rule 4: CDR preferred over a one-to-one PAN copy), else PanCanAtlas
  d1 <- cdr[cdr$type == code, , drop = FALSE]
  if (col %in% names(d1) && any(!is.na(clean(d1[[col]]))))
    return(setNames(d1[[col]], d1$bcr_patient_barcode))
  d2 <- pan[pan$acronym == code, , drop = FALSE]
  if (col %in% names(d2)) return(setNames(d2[[col]], d2$bcr_patient_barcode))
  NULL
}

summ <- list(); flags <- list()
for (code in unique(assign$cancer)) {
  a <- assign[assign$cancer == code, ]
  cols <- list()
  for (i in seq_len(nrow(a))) {
    v <- get_column(code, a$column[i])
    if (is.null(v)) { cat(sprintf("  %s: column not found -> %s\n", code, a$column[i])); next }
    cols[[a$label[i]]] <- v
  }
  if (length(cols) < 2) next
  ids <- Reduce(union, lapply(cols, names))
  X <- lapply(cols, function(v) as_analysable(v[ids]))
  labs <- names(X); dom <- setNames(a$domain[match(labs, a$label)], labs)

  # pairwise V
  p <- length(labs); V <- matrix(NA_real_, p, p, dimnames = list(labs, labs))
  for (i in seq_len(p - 1)) for (j in (i + 1):p)
    V[i, j] <- V[j, i] <- cramers_v_bc(X[[i]], X[[j]])
  write.csv(round(V, 3), file.path(OUT_DIR, paste0(code, "_pairwise_V.csv")))

  # domain-level test (construct domains only)
  for (d in unique(dom)) {
    inn <- labs[dom == d]; out <- labs[dom != d]
    type <- if (d %in% COLLECTION_DOMAINS) "collection" else "construct"
    w <- if (length(inn) > 1) V[inn, inn][upper.tri(V[inn, inn])] else numeric(0)
    b <- if (length(out)) as.vector(V[inn, out, drop = FALSE]) else numeric(0)
    mw <- if (length(w)) median(w, na.rm = TRUE) else NA_real_
    mb <- if (length(b)) median(b, na.rm = TRUE) else NA_real_
    supported <- if (type == "collection" || length(inn) < 2) NA else
      isTRUE(mw > mb) && isTRUE(mw >= MIN_V)
    summ[[length(summ) + 1]] <- data.frame(cancer = code, domain = d, type = type,
      n_vars = length(inn), median_within = round(mw, 3), median_between = round(mb, 3),
      rule1_within_gt_between = if (type == "construct" && length(inn) > 1) isTRUE(mw > mb) else NA,
      rule2_within_ge_minV    = if (type == "construct" && length(inn) > 1) isTRUE(mw >= MIN_V) else NA,
      supported = supported, stringsAsFactors = FALSE)
  }

  # variable-level flags: own construct domain vs best other construct domain
  cons <- unique(dom[!dom %in% COLLECTION_DOMAINS])
  for (v in labs[!dom %in% COLLECTION_DOMAINS]) {
    own <- setdiff(labs[dom == dom[v]], v)
    if (!length(own)) next
    m_own <- mean(V[v, own], na.rm = TRUE)
    others <- sapply(setdiff(cons, dom[v]), function(d)
      mean(V[v, labs[dom == d]], na.rm = TRUE))
    if (length(others) && any(!is.na(others)) && max(others, na.rm = TRUE) > m_own)
      flags[[length(flags) + 1]] <- data.frame(cancer = code, variable = v,
        own_domain = dom[v], mean_V_own = round(m_own, 3),
        closer_domain = names(which.max(others)), mean_V_closer = round(max(others, na.rm = TRUE), 3))
  }

  # heatmap ordered by domain
  ord <- order(dom); Vo <- V[ord, ord]; diag(Vo) <- 1
  png(file.path(OUT_DIR, paste0(code, "_heatmap.png")), width = 1400, height = 1300, res = 150)
  par(mar = c(10, 10, 3, 1))
  image(1:p, 1:p, t(Vo[p:1, ]), col = hcl.colors(20, "YlOrRd", rev = TRUE),
        zlim = c(0, 1), axes = FALSE, xlab = "", ylab = "", main = paste(code, "- Cramer's V"))
  axis(1, 1:p, colnames(Vo), las = 2, cex.axis = 0.6)
  axis(2, 1:p, rev(rownames(Vo)), las = 2, cex.axis = 0.6)
  cuts <- cumsum(table(factor(dom[ord], levels = unique(dom[ord])))) + 0.5
  abline(v = cuts, h = p + 1 - cuts, lwd = 2)
  dev.off()
  cat(sprintf("%s: %d variables, %d domains done\n", code, p, length(unique(dom))))
}

S <- do.call(rbind, summ); write.csv(S, file.path(OUT_DIR, "domain_summary.csv"), row.names = FALSE)
if (length(flags)) write.csv(do.call(rbind, flags), file.path(OUT_DIR, "variable_flags.csv"), row.names = FALSE)
cat("\nConstruct domains tested:", sum(!is.na(S$supported)),
    "| supported:", sum(S$supported, na.rm = TRUE), "\n")
print(S[!is.na(S$supported) & !S$supported, ], row.names = FALSE)
