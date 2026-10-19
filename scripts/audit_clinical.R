#Clinical audit
#Checks every clinical variable the loader can build against the eligibility filters
#(50% complete, 2+ levels, largest level at most 95%) on the exact cohort the model
#fits, using load_cancer() itself, and compares the verdict with CANCER_VARS.

if (!exists("source_library")) source("paths.R")
if (!exists("load_cancer")) source_library()

COMPL_MIN     <- 0.50  #observed for at least half of the patients
MIN_LEVELS    <- 2  #a constant variable carries no information
MAX_TOP_SHARE <- 0.95  #the 95% dominance rule

#every variable the loader can build, and the CDR-only ones (for LAML)
CANDIDATES <- names(VAR_SPEC)
CDR_VARS   <- CANDIDATES[vapply(CANDIDATES, function(v) VAR_SPEC[[v]]$src == "CDR", logical(1))]
STANDARD   <- c("age", "sex", "race", "stage", "T", "N", "M", "histology", "grade",
                "smoking", "papillary", "tumor_type")

#the cancers select_cancers.R just selected
cancers <- readLines(file.path(DIR_RESULTS, "final_cancer_list.txt"))

#the modelling cohort, built by the loader: CDR + PanCanAtlas, fallbacks, recodes
audit_cohort <- function(code) {
  load <- function(vars) suppressWarnings(suppressMessages(
    load_cancer(code, verbose = FALSE, vars = vars, read_expr = FALSE)))$clinical
  cl <- tryCatch(load(CANDIDATES), error = function(e) {
    if (grepl("no PanCanAtlas rows", conditionMessage(e))) load(CDR_VARS) else stop(e)
  })
  cl  #load_cancer already applied the survival filter and records n_rna
}

#1. completeness, levels and largest level for every variable
cat("Auditing on the modelling cohort (tumour RNA + valid survival)\n\n")
rows <- list()
for (code in cancers) {
  cl <- tryCatch(audit_cohort(code), error = function(e) {
    cat(sprintf("  %-5s FAILED: %s\n", code, conditionMessage(e))); NULL })
  if (is.null(cl) || !nrow(cl)) next
  n_rna <- attr(cl, "n_rna")
  cat(sprintf("  %-5s n = %4d  (%d with RNA, %d lost to the survival filter)\n",
              code, nrow(cl), n_rna, n_rna - nrow(cl)))
  for (v in intersect(CANDIDATES, names(cl))) {
    x  <- cl[[v]]
    cp <- round(mean(!is.na(x)), 3)
    nl <- if (is.numeric(x)) NA_integer_ else length(unique(x[!is.na(x)]))
    tb <- if (is.numeric(x)) NULL else table(x[!is.na(x)])
    ts <- if (!length(tb)) NA_real_ else max(tb) / sum(tb)  #unrounded for the 95% test
    rows[[length(rows) + 1]] <- data.frame(
      cancer = code, n_model = nrow(cl), n_rna_only = n_rna,
      variable = v, source = VAR_SPEC[[v]]$src,
      type = if (is.numeric(x)) "numeric" else "categorical",
      completeness = cp, n_levels = nl, top_share = round(ts, 4),
      pass = cp >= COMPL_MIN &&
        (is.numeric(x) || (!is.na(nl) && nl >= MIN_LEVELS && !is.na(ts) && ts <= MAX_TOP_SHARE)),
      stringsAsFactors = FALSE)
  }
}
long <- do.call(rbind, rows)
write.csv(long, file.path(DIR_RESULTS, "clinical_completeness_long.csv"), row.names = FALSE)
mat <- reshape(long[, c("cancer", "variable", "completeness")],
               idvar = "cancer", timevar = "variable", direction = "wide")
colnames(mat) <- sub("^completeness\\.", "", colnames(mat))
write.csv(mat, file.path(DIR_RESULTS, "clinical_completeness_matrix.csv"), row.names = FALSE)

#2. the audit verdict next to what CANCER_VARS and the extended set use
ext_of <- function(code) if (exists("EXT_VARS")) unlist(EXT_VARS[[code]], use.names = FALSE) else character(0)
cmp <- list()
for (code in unique(long$cancer)) {
  used   <- names(CANCER_VARS[[code]])
  ext    <- ext_of(code)
  sub    <- long[long$cancer == code, ]
  passes <- sub$variable[sub$pass]
  for (v in union(used, passes)) {
    r <- sub[sub$variable == v, ]
    cmp[[length(cmp) + 1]] <- data.frame(
      cancer = code, variable = v, standard = v %in% STANDARD,
      in_main_set = v %in% used, in_extended_set = v %in% ext, passes_audit = v %in% passes,
      completeness = if (nrow(r)) r$completeness else NA,
      n_levels = if (nrow(r)) r$n_levels else NA,
      top_share = if (nrow(r)) r$top_share else NA,
      verdict = if (v %in% used && v %in% passes) "ok"
      else if (v %in% used) "IN MAIN SET BUT FAILS"
      else if (v %in% ext) "extended set"
      else "PASSES BUT IN NEITHER SET",
      stringsAsFactors = FALSE)
  }
}
cmp <- do.call(rbind, cmp)
write.csv(cmp, file.path(DIR_RESULTS, "audit_vs_cancervars.csv"), row.names = FALSE)

#3. report
cat("\n================ AUDIT vs CANCER_VARS ================\n")
print(table(cmp$verdict))
show <- function(d, title) {
  if (!nrow(d)) return(invisible())
  cat(sprintf("\n%s (%d):\n", title, nrow(d)))
  print(d[order(d$variable, d$cancer), c("cancer", "variable", "completeness", "top_share")],
        row.names = FALSE)
}
show(cmp[cmp$verdict == "IN MAIN SET BUT FAILS", ],
     "In a main set but failing the filters on the modelling cohort")
show(cmp[cmp$verdict != "ok" & cmp$standard & !cmp$in_main_set, ],
     "Standard variables that pass the filters but are not in the main set")
show(cmp[cmp$verdict == "PASSES BUT IN NEITHER SET" & !cmp$standard, ],
     "Other variables that pass the filters but are in neither set")
cat("\nWritten: clinical_completeness_long.csv, _matrix.csv, audit_vs_cancervars.csv\n")
cat("======================================================\n")