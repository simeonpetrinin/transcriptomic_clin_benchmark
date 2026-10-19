#Discover clinical variables
#Screens every raw clinical column in TCGA-CDR and PanCanAtlas on the modelling cohort (tumour RNA, valid OS)
#to find candidate variables; it does not choose the final variable sets (R/enumerate.R holds those).
#Completeness is measured on the raw column, before the loader's recoding, so it can differ from
#audit_clinical.R (e.g. ESCA grade: raw passes, 78.3% after GX/GB -> missing). For variables already in
#VAR_SPEC, post_transform_completeness gives the modelled figure.
#Tiers: PASS_80 (>=80% complete and informative), REVIEW_50_80 (50-80%, reviewed by hand).
#Outcome, post-baseline, technical and identifier flags are warnings, not exclusions.
#Writes results/clinical_discovery_all.csv and results/clinical_discovery_candidates.csv.

source("paths.R")

# Optional: brings in VAR_SPEC so post-transform completeness can be reported.
# The script runs without it -- that column is simply left NA.
HAVE_LOADER <- tryCatch({ source_library(); exists("VAR_SPEC") },
                        error = function(e) FALSE)

suppressPackageStartupMessages({
  library(readxl)
  library(data.table)
})

COMPL_PASS     <- 0.80
COMPL_REVIEW   <- 0.50
MAX_TOP_SHARE  <- 0.95
NUMERIC_FRAC   <- 0.95
N_EXAMPLES     <- 6L
N_TOP_VALUES   <- 5L


# Missing-value normalisation

SENTINELS <- c(
  "", "na", "n/a", "#n/a",
  "not available", "[not available]",
  "not applicable", "[not applicable]",
  "unknown", "[unknown]",
  "not evaluated", "[not evaluated]",
  "pending", "[pending]",
  "discrepancy", "[discrepancy]",
  "--"
)

# NOTE: some PanCanAtlas columns contain bytes that are not valid UTF-8 (free-text
# fields typed at the tissue source sites). trimws() and tolower() use perl
# regexes, which refuse invalid input with "input string N is invalid UTF-8" --
# the run died at PCPG. iconv() drops the offending bytes rather than failing;
# the affected values are free text we only ever count, never parse, so nothing
# meaningful is lost. Applied first, before anything touches the string.
.to_utf8 <- function(x) {
  z <- as.character(x)
  bad <- is.na(iconv(z, "UTF-8", "UTF-8"))
  if (any(bad, na.rm = TRUE))
    z[which(bad)] <- iconv(z[which(bad)], "UTF-8", "UTF-8", sub = "")
  z
}

clean_raw <- function(x) {
  z  <- trimws(.to_utf8(x))
  zl <- tolower(z)
  z[is.na(z) | zl %in% SENTINELS] <- NA_character_
  z
}


# Semantic FLAGS ONLY -- these do not remove columns from the complete output
# NOTE: these patterns previously ended in a dangling "|", which in a regex means
# "or the empty string" and therefore matches EVERY column. flag_outcome also had
# an unclosed "(" and would have thrown "invalid regular expression" on the first
# call. Patterns are now built from character vectors and joined with paste, so a
# missing or trailing element cannot silently turn into "match everything".

RX_IDENTIFIER <- c("barcode", "uuid", "patient_id", "sample_id", "case_id")

RX_OUTCOME <- c("^os$", "os[._]?time", "survival", "vital_status",
                "days_to_death", "follow.?up", "tumou?r_status",
                "recurrence", "progression", "disease_free", "new_tumou?r")

RX_POSTBASELINE <- c("residual", "margin", "treatment", "therapy", "radiation",
                     "response", "surgery", "surgical", "postoperative",
                     "post_operative")

RX_TECHNICAL <- c("method", "protocol", "platform", "calculation_method",
                  "assay", "aliquot", "shipment", "batch", "slide")

.any_match <- function(z, pats) grepl(paste(pats, collapse = "|"), z)

column_flags <- function(nm) {
  z <- tolower(.to_utf8(nm))
  data.frame(
    flag_identifier   = grepl(paste0("(^|_)(", paste(RX_IDENTIFIER, collapse = "|"),
                                     ")(_|$)"), z),
    flag_outcome      = .any_match(z, RX_OUTCOME),
    flag_postbaseline = .any_match(z, RX_POSTBASELINE),
    flag_technical    = .any_match(z, RX_TECHNICAL),
    stringsAsFactors  = FALSE
  )
}

# fail loudly at load time if a pattern is malformed, rather than mid-run
local({
  chk <- try(column_flags("bcr_patient_barcode"), silent = TRUE)
  if (inherits(chk, "try-error"))
    stop("column_flags() regex is invalid: ", conditionMessage(attr(chk, "condition")))
  if (!chk$flag_identifier || chk$flag_postbaseline)
    stop("column_flags() sanity check failed -- a pattern is matching too much ",
         "or too little. Check for a trailing '|' in one of the RX_ vectors.")
})


# Post-transform completeness, for variables the loader already knows
# Maps a raw column back to its VAR_SPEC entry (if any) and reports completeness
# AFTER the loader's collapsing rule. This is the number the model actually sees.
post_transform <- function(raw_col, source, x, code, n_model) {
  if (!HAVE_LOADER) return(c(var = NA_character_, compl = NA_real_))
  src <- if (source == "PANCAN") "PAN" else "CDR"
  hit <- NULL
  for (v in names(VAR_SPEC))
    if (VAR_SPEC[[v]]$src == src && identical(VAR_SPEC[[v]]$col, raw_col)) { hit <- v; break }
  if (is.null(hit)) return(c(var = NA_character_, compl = NA_real_))
  z <- tryCatch(VAR_SPEC[[hit]]$fn(.to_utf8(x), code),
                error = function(e) rep(NA, length(x)))
  c(var = hit, compl = if (n_model) sum(!is.na(z)) / n_model else NA_real_)
}


# Profile ONE raw variable on ONE cancer cohort

profile_variable <- function(x, cancer, source, raw_column, n_model) {
  
  raw_class <- class(x)[1]
  z   <- clean_raw(x)
  obs <- z[!is.na(z)]
  n_nonmissing <- length(obs)
  
  completeness <- if (n_model) n_nonmissing / n_model else NA_real_
  
  if (!n_nonmissing) {
    n_unique <- 0L; top_share <- NA_real_; numeric_fraction <- NA_real_
    examples <- NA_character_; top_values <- NA_character_
  } else {
    tb <- sort(table(obs), decreasing = TRUE)
    n_unique  <- length(tb)
    top_share <- max(tb) / n_nonmissing
    parsed    <- suppressWarnings(as.numeric(obs))
    numeric_fraction <- mean(!is.na(parsed))
    examples  <- paste(head(unique(obs), N_EXAMPLES), collapse = " | ")
    take <- seq_len(min(length(tb), N_TOP_VALUES))
    top_values <- paste(
      sprintf("%s (%d; %.1f%%)", names(tb)[take], as.integer(tb[take]),
              100 * as.numeric(tb[take]) / n_nonmissing),
      collapse = " | ")
  }
  
  inferred_type <- if (!is.na(numeric_fraction) && numeric_fraction >= NUMERIC_FRAC)
    "numeric-like" else "categorical/text"
  
  informative <- !is.na(top_share) && n_unique >= 2 && top_share <= MAX_TOP_SHARE
  
  # NOTE: tiering uses the UNROUNDED completeness. Rounding first would let
  # 0.7995 round to 0.800 and pass an 0.80 threshold it actually fails.
  tier <- if (informative && completeness >= COMPL_PASS) "PASS_80"
  else if (informative && completeness >= COMPL_REVIEW) "REVIEW_50_80"
  else "OTHER"
  
  pt <- post_transform(raw_column, source, x, cancer, n_model)
  
  cbind(
    data.frame(
      cancer = cancer, source = source, raw_column = raw_column,
      n_model = n_model, n_nonmissing = n_nonmissing,
      completeness = round(completeness, 4),
      var_spec_name = unname(pt["var"]),
      post_transform_completeness = round(as.numeric(pt["compl"]), 4),
      n_unique = n_unique, top_share = round(top_share, 4),
      raw_class = raw_class,
      numeric_fraction = round(numeric_fraction, 4),
      inferred_type = inferred_type,
      example_values = examples, top_values = top_values,
      tier = tier, stringsAsFactors = FALSE
    ),
    column_flags(raw_column)
  )
}


# Load source data ONCE

cat("Reading TCGA-CDR...\n")
cdr <- as.data.frame(read_xlsx(CDR_PATH, sheet = "TCGA-CDR", guess_max = 11000))

cat("Reading complete PanCan clinical table...\n")
pan <- fread(PANCAN_PATH, sep = "\t", header = TRUE,
             data.table = FALSE, showProgress = TRUE)


# Identify patients with tumour RNA

cat("Reading tumour RNA barcodes...\n")
con <- gzfile(RNA_PATH); header <- readLines(con, n = 1); close(con)
samples     <- strsplit(header, "\t")[[1]][-1]
sample_type <- suppressWarnings(as.integer(substr(samples, 14, 15)))
tumour_samples <- samples[!is.na(sample_type) & sample_type < 10]
rna_patients   <- unique(substr(tumour_samples, 1, 12))


# Sanity checks

required_cdr <- c("bcr_patient_barcode", "type", "OS", "OS.time")
required_pan <- c("bcr_patient_barcode", "acronym")

missing_cdr <- setdiff(required_cdr, names(cdr))
missing_pan <- setdiff(required_pan, names(pan))
if (length(missing_cdr))
  stop("Missing required CDR columns: ", paste(missing_cdr, collapse = ", "))
if (length(missing_pan))
  stop("Missing required PANCAN columns: ", paste(missing_pan, collapse = ", "))
if (anyDuplicated(cdr$bcr_patient_barcode))
  stop("TCGA-CDR contains duplicated patient barcodes.")
if (anyDuplicated(pan$bcr_patient_barcode))
  stop("PANCAN contains duplicated patient barcodes.")


# Exhaustive per-cancer discovery

codes <- sort(unique(cdr$type))
rows  <- list()

for (code in codes) {
  
  cat(sprintf("\n%-5s ", code))
  cdr_c <- cdr[cdr$type == code, , drop = FALSE]
  
  os      <- suppressWarnings(as.integer(as.character(cdr_c$OS)))
  os_time <- suppressWarnings(as.numeric(cdr_c$OS.time))
  
  keep <- cdr_c$bcr_patient_barcode %in% rna_patients &
    !is.na(os) & !is.na(os_time) & os_time > 0
  
  cdr_c    <- cdr_c[keep, , drop = FALSE]
  model_ids <- cdr_c$bcr_patient_barcode
  n_model   <- length(model_ids)
  cat(sprintf("n = %d", n_model))
  if (!n_model) next
  
  #CDR: every raw column
  for (nm in names(cdr_c))
    rows[[length(rows) + 1L]] <- profile_variable(
      x = cdr_c[[nm]], cancer = code, source = "CDR",
      raw_column = nm, n_model = n_model)
  
  #PANCAN: every raw column, aligned to the SAME patients
  # match() deliberately leaves NA for a model patient absent from PANCAN, so
  # such a patient counts AGAINST completeness rather than silently dropping out
  # of the denominator. That is the cost of the two-source join, made visible.
  pan_c <- pan[pan$acronym == code, , drop = FALSE]
  idx   <- match(model_ids, pan_c$bcr_patient_barcode)
  
  for (nm in names(pan_c))
    rows[[length(rows) + 1L]] <- profile_variable(
      x = pan_c[[nm]][idx], cancer = code, source = "PANCAN",
      raw_column = nm, n_model = n_model)
}


# Save complete discovery table

out <- do.call(rbind, rows)
tier_order <- match(out$tier, c("PASS_80", "REVIEW_50_80", "OTHER"))
out <- out[order(out$cancer, tier_order, -out$completeness,
                 out$source, out$raw_column), ]

write.csv(out, file.path(DIR_RESULTS, "clinical_discovery_all.csv"),
          row.names = FALSE)


# Human-review shortlist
#
# Identifiers and obvious outcome fields are omitted from THIS convenience table
# only. They remain in clinical_discovery_all.csv.
#
# Post-baseline and technical-looking fields remain visible and FLAGGED because
# they require scientific judgment rather than blind regex exclusion.

cand <- out[out$tier %in% c("PASS_80", "REVIEW_50_80") &
              !out$flag_identifier & !out$flag_outcome, ]

write.csv(cand, file.path(DIR_RESULTS, "clinical_discovery_candidates.csv"),
          row.names = FALSE)


# Report

cat("\n\n================ DISCOVERY COMPLETE ================\n")
cat("Raw columns screened:", length(unique(paste(out$source, out$raw_column))), "\n")
cat("Cancer-column combinations:", nrow(out), "\n")
cat("PASS_80:", sum(out$tier == "PASS_80"), "\n")
cat("REVIEW_50_80:", sum(out$tier == "REVIEW_50_80"), "\n")

if (HAVE_LOADER) {
  drift <- out[!is.na(out$post_transform_completeness) &
                 abs(out$completeness - out$post_transform_completeness) > 0.02, ]
  cat(sprintf("\nColumns where the loader's transform moves completeness by >2pp: %d\n",
              nrow(drift)))
  if (nrow(drift)) {
    d <- drift[order(-abs(drift$completeness - drift$post_transform_completeness)), ]
    cat("  worst offenders (raw -> post-transform):\n")
    for (i in seq_len(min(8, nrow(d))))
      cat(sprintf("    %-5s %-34s %.3f -> %.3f\n", d$cancer[i], d$raw_column[i],
                  d$completeness[i], d$post_transform_completeness[i]))
    cat("  These are why this script and audit_clinical.R disagree. Neither is\n")
    cat("  wrong: discovery screens raw, the audit screens what the model sees.\n")
  }
} else {
  cat("\n(loader not available -- post_transform_completeness left NA)\n")
}

cat("\nWritten:\n",
    "  results/clinical_discovery_all.csv\n",
    "  results/clinical_discovery_candidates.csv\n")
cat("====================================================\n")