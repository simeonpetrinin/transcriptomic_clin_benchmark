#Audit all variables
#Completeness and largest-category share of every raw clinical column in TCGA-CDR and PanCanAtlas, per cancer.
#Raw columns, before the loader's recoding (GX, TX, NX, Tis -> missing), so this is a wide screen; audit_clinical.R
#gives the modelled figures. Run from the project root: source("scripts/audit_all_variables.R").
#Writes results/clinical_gate_scan.csv (one row per cancer x column).

source("paths.R")
source("R/load_data.R")     # SENT, mark_na, lc, VAR_SPEC
source("R/enumerate.R")     # CANCER_VARS

COMPL_MIN <- 0.80                 # completeness gate
TOP_SHARE_MAX <- 0.95             # informativeness gate: largest level <= 95%
NEAR_MISS <- 0.70                 # flag columns between NEAR_MISS and COMPL_MIN
MIN_LEVEL_SHARE <- 0.05           # for the alternative "two levels >= 5%" rule

codes <- tryCatch(readLines(file.path(DIR_RESULTS, "final_cancer_list.txt")),
                  error = function(e) names(CANCER_VARS))

# col_types = "text": otherwise readxl guesses a column's type from its first
# rows, and a column that is empty early on (e.g. residual_tumor) gets guessed
# as logical, turning later values such as "R0" into NA with a warning.
cdr <- as.data.frame(readxl::read_xlsx(CDR_PATH, sheet = "TCGA-CDR",
                                       col_types = "text"))
pan <- as.data.frame(data.table::fread(PANCAN_PATH, colClasses = "character",
                                       data.table = FALSE))

# which raw columns already feed the lattice, per cancer
lattice_cols <- function(code) {
  vars <- names(CANCER_VARS[[code]])
  if (is.null(vars)) return(character(0))
  unlist(lapply(vars, function(v) {
    sp <- VAR_SPEC[[v]]
    if (is.null(sp)) NA_character_ else paste0(sp$src, ":", sp$col)
  }))
}

# summarise one column on one cohort
describe <- function(x, n_cohort) {
  v <- lc(x)                                   # loader's sentinel rules
  ok <- !is.na(v)
  compl <- if (n_cohort > 0) sum(ok) / n_cohort else NA_real_
  num <- suppressWarnings(as.numeric(v[ok]))
  is_num <- length(num) > 0 && mean(!is.na(num)) > 0.95 &&
    length(unique(num[!is.na(num)])) > 15
  if (is_num) {
    data.frame(type = "numeric", n_nonmissing = sum(ok), completeness = compl,
               n_levels = NA_integer_, n_levels_ge5pct = NA_integer_,
               top_level = NA_character_, top_share = NA_real_,
               second_share = NA_real_, level_summary = NA_character_,
               numeric_median = suppressWarnings(median(num, na.rm = TRUE)),
               numeric_min = suppressWarnings(min(num, na.rm = TRUE)),
               numeric_max = suppressWarnings(max(num, na.rm = TRUE)))
  } else {
    tb <- sort(table(v[ok]), decreasing = TRUE)
    sh <- if (length(tb)) as.numeric(tb) / sum(tb) else numeric(0)
    top3 <- if (length(tb)) paste(sprintf("%s (%.0f%%)", names(tb)[1:min(3, length(tb))],
                                          100 * sh[1:min(3, length(tb))]),
                                  collapse = "; ") else NA_character_
    data.frame(type = "categorical", n_nonmissing = sum(ok), completeness = compl,
               n_levels = length(tb), n_levels_ge5pct = sum(sh >= MIN_LEVEL_SHARE),
               top_level = if (length(tb)) names(tb)[1] else NA_character_,
               top_share = if (length(sh)) sh[1] else NA_real_,
               second_share = if (length(sh) > 1) sh[2] else NA_real_,
               level_summary = top3,
               numeric_median = NA_real_, numeric_min = NA_real_, numeric_max = NA_real_)
  }
}

rows <- list()
for (code in codes) {
  f <- path_modelling(code)
  if (!file.exists(f)) { message(code, ": no modelling file, skipped"); next }
  pts <- rownames(readRDS(f)$expr)             # the analysed cohort
  n <- length(pts)
  in_lat <- lattice_cols(code)
  
  for (src in c("CDR", "PAN")) {
    df  <- if (src == "CDR") cdr else pan
    key <- df$bcr_patient_barcode
    sub <- df[match(pts, key), , drop = FALSE]
    for (cl in setdiff(names(df), "bcr_patient_barcode")) {
      d <- describe(sub[[cl]], n)
      rows[[length(rows) + 1]] <- cbind(
        data.frame(cancer = code, source = src, column = cl, n_cohort = n),
        d,
        data.frame(
          passes_completeness  = !is.na(d$completeness) & d$completeness >= COMPL_MIN,
          passes_informative   = d$type == "numeric" |
            (!is.na(d$n_levels) & d$n_levels >= 2 &
               !is.na(d$top_share) & d$top_share <= TOP_SHARE_MAX),
          passes_alt_5pct      = d$type == "numeric" |
            (!is.na(d$n_levels_ge5pct) & d$n_levels_ge5pct >= 2),
          near_miss            = !is.na(d$completeness) &
            d$completeness >= NEAR_MISS &
            d$completeness <  COMPL_MIN,
          in_lattice           = paste0(src, ":", cl) %in% in_lat))
    }
  }
  message(code, ": scanned")
}

scan <- do.call(rbind, rows)
scan$passes_both <- scan$passes_completeness & scan$passes_informative
scan$fail_reason <- ifelse(scan$passes_both, "",
                           ifelse(!scan$passes_completeness & !scan$passes_informative, "completeness+informative",
                                  ifelse(!scan$passes_completeness, "completeness", "informative")))

dir.create(DIR_RESULTS, showWarnings = FALSE, recursive = TRUE)
write.csv(scan, file.path(DIR_RESULTS, "clinical_gate_scan.csv"), row.names = FALSE)

# per-cancer counts
scan$pass_unused <- scan$passes_both & !scan$in_lattice
tab <- aggregate(cbind(passes_completeness, passes_both, pass_unused,
                       near_miss, in_lattice) ~ cancer, data = scan, FUN = sum)
names(tab) <- c("cancer", "pass_80pct", "pass_both_gates", "pass_but_unused",
                "near_miss_70_80", "used_in_lattice")
print(tab[order(-tab$pass_both_gates), ], row.names = FALSE)

# which columns clear both gates in many cohorts but are never used anywhere
cross <- aggregate(cbind(passes_both, in_lattice) ~ source + column,
                   data = scan, FUN = sum)
names(cross) <- c("source", "column", "n_cancers_pass", "n_cancers_used")
cross <- cross[cross$n_cancers_pass > 0 & cross$n_cancers_used == 0, ]
cat("\nColumns passing both gates in several cohorts but never in any lattice:\n")
print(head(cross[order(-cross$n_cancers_pass), ], 40), row.names = FALSE)

#write one workbook, one sheet per cancer
# Sheet order: summary, cross_cancer, then one sheet per cohort.
sheet_cols <- c("source", "column", "type", "completeness", "n_nonmissing",
                "n_levels", "n_levels_ge5pct", "top_level", "top_share",
                "second_share", "level_summary", "numeric_median",
                "numeric_min", "numeric_max", "passes_completeness",
                "passes_informative", "passes_alt_5pct", "passes_both",
                "fail_reason", "near_miss", "in_lattice")

sheets <- c(list(summary = tab, cross_cancer = cross),
            lapply(split(scan, scan$cancer), function(d) {
              d <- d[order(-d$passes_both, -d$in_lattice, -d$completeness), sheet_cols]
              d$completeness <- round(d$completeness, 3)
              d$top_share    <- round(d$top_share, 3)
              d$second_share <- round(d$second_share, 3)
              d
            }))

xlsx_out <- file.path(DIR_RESULTS, "clinical_gate_scan.xlsx")
if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb  <- openxlsx::createWorkbook()
  hdr <- openxlsx::createStyle(textDecoration = "bold")
  for (nm in names(sheets)) {
    openxlsx::addWorksheet(wb, nm)
    openxlsx::writeData(wb, nm, sheets[[nm]], headerStyle = hdr)
    openxlsx::freezePane(wb, nm, firstRow = TRUE)
    openxlsx::setColWidths(wb, nm, cols = seq_along(sheets[[nm]]), widths = "auto")
  }
  openxlsx::saveWorkbook(wb, xlsx_out, overwrite = TRUE)
} else if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(sheets, xlsx_out)   # no freeze panes or column widths
} else {
  stop("Install openxlsx (preferred) or writexl to write the workbook. ",
       "This project uses renv, so use renv::install(\"openxlsx\") and ",
       "renv::snapshot() rather than install.packages().")
}
cat("\nwritten:", xlsx_out, "\n")
cat("written:", file.path(DIR_RESULTS, "clinical_gate_scan.csv"), "\n\n")

# columns that clear both gates but are NOT in the lattice, for one cancer
show_cancer <- function(code, only_unused = TRUE) {
  s <- scan[scan$cancer == code & scan$passes_both, ]
  if (only_unused) s <- s[!s$in_lattice, ]
  s <- s[order(-s$completeness), c("source", "column", "completeness", "type",
                                   "n_levels", "n_levels_ge5pct", "top_share",
                                   "level_summary", "in_lattice")]
  s$completeness <- round(s$completeness, 3)
  s$top_share <- round(s$top_share, 3)
  print(head(s, 40), row.names = FALSE)
}

# every column for one cancer, whether it passes or not, with the failure reason
show_all <- function(code) {
  s <- scan[scan$cancer == code, ]
  s <- s[order(-s$passes_both, -s$completeness),
         c("source", "column", "completeness", "type", "n_levels", "top_share",
           "passes_both", "fail_reason", "near_miss", "in_lattice")]
  s$completeness <- round(s$completeness, 3); s$top_share <- round(s$top_share, 3)
  s
}

cat("\nCESC: columns clearing both gates that are not in the lattice\n")
show_cancer("CESC")