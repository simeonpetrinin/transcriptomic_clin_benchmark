#Counts the clinical fields each external cohort records, straight from the raw source files.
#A field counts if it is at least MIN_COMPLETE filled in and has 2+ distinct values; IDs, admin
#and survival/follow-up fields are left out, and each kept field is tagged clinical/treatment/molecular.
#Run from the project root: source("scripts/count_external_clinical.R"). Reads only; writes results/external/clinical_field_counts*.csv.

RAW          <- file.path("data", "raw", "external")
OUT_DIR      <- file.path("results", "external")
MIN_COMPLETE <- 0.50

NA_RX    <- "^\\s*(\\[.*\\]|na|n/a|nan|unknown|not available|not reported|not applicable|-|--|\\.|none|null|nx|tx|mx|undetermined|declined)?\\s*$"
ADMIN_RX <- "(^|[_. ])(id|uuid|barcode|accession|gsm|lcs|sample|specimen|batch|platform|instrument|dataset|organism|individual|tissue|source|cell type|strain|genotype|cre|indexes|cohort|center|centerid|dbgap|manuscript|analysis|rnaseq|exomeseq|totaldrug|oncotree|cancer_type|somatic|gene_panel|extraction|dependancy|disease.?state|organism part)($|[_. ])|^specimen|^analysis"
OUTC_RX  <- "vital|death|days_to|survival|^os$|^os[_. ]|^rfs|^pfs|follow|event|recurr|relapse|progression|cause|months|mths|delay|last_contact|censor|duration|outcome|acquisition|currentstage|currentregimen|mostrecent|cumulativetreatment|responseto|timeofsample|^inferred"
TREAT_RX <- "therap|chemo|radiat|radio|treat|surg|adjuvant|regimen|transplant|induction|drug|margin|resection|residual|endocrine|hormone|intravesical"
MOLEC_RX <- "molecularsubtype|mutation|idh|mgmt|1p19q|codel|hpv|mmr|cimp|cin\\.|kras|braf|tp53|pam50|cluster|intclust|claudin|threegene|prediction|flt3|npm1|karyotype|cytogenetic|fusion|variant|tmb|allelic|cebpa|runx1|asxl1|signature|her2_snp6|eln|ki67|vhl"

#GEO series matrix -> one row per sample, one column per characteristic ("key: value")
read_geo <- function(path) {
  con <- gzfile(path, "r"); on.exit(close(con))
  ch <- character(0)
  repeat {
    l <- readLines(con, n = 1, warn = FALSE)
    if (!length(l) || startsWith(l, "!series_matrix_table_begin")) break
    if (startsWith(l, "!Sample_characteristics_ch1")) ch <- c(ch, l)
  }
  rows <- lapply(ch, function(l) gsub('^"|"$', "", strsplit(l, "\t", fixed = TRUE)[[1]][-1]))
  n <- length(rows[[1]]); out <- list()
  for (r in rows) for (i in seq_len(n)) if (grepl(":", r[i], fixed = TRUE)) {
    k <- tolower(trimws(sub(":.*$", "", r[i]))); v <- trimws(sub("^[^:]*:", "", r[i]))
    if (is.null(out[[k]])) out[[k]] <- rep("", n)
    out[[k]][i] <- v
  }
  as.data.frame(out, check.names = FALSE, stringsAsFactors = FALSE)
}

#cBioPortal clinical table (skips the 4 "#" header lines)
read_cbio <- function(path) read.delim(path, comment.char = "#", colClasses = "character", check.names = FALSE)

#tag and count the fields of one cohort
count_fields <- function(cohort, d) {
  rows <- lapply(names(d), function(f) {
    x <- trimws(as.character(d[[f]])); ok <- !is.na(x) & !grepl(NA_RX, x, ignore.case = TRUE)
    key <- tolower(f)
    kind <- if (grepl(ADMIN_RX, key)) "excluded: id/admin" else if (grepl(OUTC_RX, key)) "excluded: outcome/follow-up" else
      if (grepl(TREAT_RX, key)) "treatment" else if (grepl(MOLEC_RX, key)) "molecular" else "clinical"
    tb <- table(x[ok])
    data.frame(cohort = cohort, field = f, kind = kind, complete = round(mean(ok), 3),
               n_values = length(tb), top_share = if (length(tb)) round(max(tb) / sum(ok), 3) else NA_real_,
               counted = !startsWith(kind, "excluded") && mean(ok) >= MIN_COMPLETE && length(tb) >= 2,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#each cohort: raw clinical table restricted to its tumour patients
cohorts <- list(
  GSE68465 = function() { d <- read_geo(file.path(RAW, "GSE68465", "GSE68465_series_matrix.txt.gz")); d[grepl("adeno", d$disease_state, ignore.case = TRUE), ] },
  GSE39582 = function() { d <- read_geo(file.path(RAW, "GSE39582", "GSE39582_series_matrix.txt.gz")); d[d$dataset != "Non Tumoral", ] },
  GSE65858 = function() read_geo(file.path(RAW, "GSE65858", "GSE65858_series_matrix.txt.gz")),
  GSE13507 = function() { d <- read_geo(file.path(RAW, "GSE13507", "GSE13507_series_matrix.txt.gz")); d[d[["biological source"]] == "primary bladder cancer", ] },
  GSE14520 = function() { d <- read.delim(gzfile(file.path(RAW, "GSE14520", "GSE14520_Extra_Supplement.txt.gz")), colClasses = "character", check.names = FALSE)
  d <- d[d[["Tissue Type"]] == "Tumor", ]; d[!duplicated(d$ID), ] },
  GSE84437 = function() { d <- read_geo(file.path(RAW, "GSE84437", "GSE84437_series_matrix.txt.gz")); d[d$age != "", ] },
  METABRIC = function() { p <- read_cbio(file.path(RAW, "METABRIC", "brca_metabric", "data_clinical_patient.txt"))
  s <- read_cbio(file.path(RAW, "METABRIC", "brca_metabric", "data_clinical_sample.txt"))
  m <- merge(p, s, by = "PATIENT_ID"); m[m$PATIENT_ID %in% readRDS(file.path("data", "processed", "external", "METABRIC_modelling.rds"))$clinical$patient_id, ] },
  CGGA_LGG = function() { z <- file.path(RAW, "CGGA", "CGGA.mRNAseq_693_clinical.20200506.txt.zip")
  d <- read.delim(unz(z, unzip(z, list = TRUE)$Name[1]), colClasses = "character", check.names = FALSE)
  d[d$PRS_type == "Primary" & d$Grade %in% c("WHO II", "WHO III"), ] },
  CGGA_GBM = function() { z <- file.path(RAW, "CGGA", "CGGA.mRNAseq_693_clinical.20200506.txt.zip")
  d <- read.delim(unz(z, unzip(z, list = TRUE)$Name[1]), colClasses = "character", check.names = FALSE)
  d[d$PRS_type == "Primary" & d$Grade == "WHO IV", ] },
  BEATAML  = function() { d <- as.data.frame(readxl::read_excel(file.path(RAW, "BEATAML", "beataml_wv1to4_clinical.xlsx"), col_types = "text"))
  d <- d[d$rnaSeq == "y" & d$diseaseStageAtSpecimenCollection %in% "Initial Diagnosis", ]; d[!duplicated(d$dbgap_subject_id), ] }
)

all_fields <- list(); summ <- list()
for (id in names(cohorts)) {
  d <- tryCatch(cohorts[[id]](), error = function(e) { cat(sprintf("%-9s could not be read: %s\n", id, conditionMessage(e))); NULL })
  if (is.null(d)) next
  f <- count_fields(id, d); all_fields[[id]] <- f
  k <- f[f$counted, ]
  summ[[id]] <- data.frame(cohort = id, source = if (startsWith(id, "GSE")) "GEO" else "other", patients = nrow(d),
                           clinical = sum(k$kind == "clinical"), treatment = sum(k$kind == "treatment"),
                           molecular = sum(k$kind == "molecular"), clinical_fields = paste(k$field[k$kind == "clinical"], collapse = "; "))
}
summ <- do.call(rbind, summ); fields <- do.call(rbind, all_fields)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
write.csv(summ, file.path(OUT_DIR, "clinical_field_counts.csv"), row.names = FALSE)
write.csv(fields, file.path(OUT_DIR, "clinical_field_counts_by_field.csv"), row.names = FALSE)
print(summ[, c("cohort", "source", "patients", "clinical", "treatment", "molecular")], row.names = FALSE)
geo <- summ[summ$source == "GEO", ]
cat(sprintf("\nGEO cohorts with at most 10 clinical fields (>= %.0f%% complete): %d of %d\n",
            100 * MIN_COMPLETE, sum(geo$clinical <= 10), nrow(geo)))
cat("Check the tags in clinical_field_counts_by_field.csv: they come from field names, so a few may be wrong.\n")