#TCGA feasibility audit
#Per cancer: patients with both tumour RNA and CDR clinical data, and the events and
#event ratio for OS and PFI.

if (!exists("DIR_RAW")) source("paths.R")  #standalone-safe

#CDR_PATH and RNA_PATH come from paths.R
dir.create(DIR_RESULTS, recursive = TRUE, showWarnings = FALSE)
out_csv <- file.path(DIR_RESULTS, "tcga_feasibility_audit.csv")

#missing-value handling (CDR sentinels)
is_missing <- function(x) {
  x <- trimws(tolower(as.character(x)))
  is.na(x) | x %in% c(
    "", "na", "n/a", "#n/a", "not available", "[not available]",
    "[not applicable]", "not applicable", "[unknown]", "unknown",
    "[discrepancy]", "[not evaluated]", "[pending]", "--"
  )
}

#Calculate completeness percentage
completeness <- function(x) {
  if (length(x) == 0L) return(NA_real_)
  round(mean(!is_missing(x)), 3)
}

#Pathologic value where recorded, clinical values serve as fallback
coalesce_stage <- function(pathologic, clinical) {
  out <- as.character(pathologic)
  out[is_missing(out)] <- as.character(clinical)[is_missing(out)]
  out
}

#Events and event ratio for one endpoint.
count_events <- function(ev, t) {
  usable <- ev %in% c(0, 1) & !is.na(t) & t > 0
  n_usable <- sum(usable)
  events <- sum(ev[usable] == 1)
  ratio <- if (n_usable > 0) {
    round(events / n_usable, 3)
  } else {
    NA_real_
  }
  list(events = events, ratio = ratio)
}

#1. read the RNA barcodes from the Xena header
cat("Reading RNA barcodes from Xena header ...\n")
connection <- gzfile(RNA_PATH)
header <- readLines(connection, n = 1)  #first line only = all sample IDs
close(connection)

samps <- strsplit(header, "\t")[[1]][-1]  #drop the leading "sample" label
stype  <- suppressWarnings(as.integer(substr(samps, 14, 15)))
tumour <- samps[!is.na(stype) & stype < 10]  #Keeps all tumour codes 01-09 (not only primary 01), so SKCM
rna_patients <- unique(toupper(trimws(substr(tumour, 1, 12))))  #dedup to 12-char patient IDs, cleaned like the CDR IDs
cat("  RNA samples:", length(samps),  #RNA Samples
    " | tumour samples:", length(tumour),  #Tumor (<10) Samples
    " | unique tumour patients:", length(rna_patients), "\n\n")  #Unique tumor patients

#2. read the CDR clinical/survival table
cat("Reading TCGA-CDR ...\n")

cdr <- readxl::read_xlsx(
  CDR_PATH,
  sheet = "TCGA-CDR",
  guess_max = 11000,
  .name_repair = "unique_quiet"
)

cdr$barcode <- toupper(trimws(cdr$bcr_patient_barcode))
cdr$OS_ev   <- suppressWarnings(as.numeric(as.character(cdr$OS)))
cdr$OS_t    <- suppressWarnings(as.numeric(as.character(cdr$OS.time)))
cdr$PFI_ev  <- suppressWarnings(as.numeric(as.character(cdr$PFI)))
cdr$PFI_t   <- suppressWarnings(as.numeric(as.character(cdr$PFI.time)))

cancer_codes <- sort(unique(cdr$type))
cat("  CDR covers", length(cancer_codes), "cancer types,", nrow(cdr), "patients.\n\n")

#3. per-cancer audit on the analysis cohort
audit <- data.frame()
for (code in cancer_codes) {
  sub <- cdr[which(cdr$type == code), , drop = FALSE]
  ana <- sub[sub$barcode %in% rna_patients, ]  #CDR patients WITH tumour RNA
  
  n_ana <- nrow(ana)
  os    <- count_events(ana$OS_ev,  ana$OS_t)
  pfi   <- count_events(ana$PFI_ev, ana$PFI_t)
  
  audit <- rbind(audit, data.frame(
    cancer          = code,
    n_cdr           = nrow(sub),
    n_analysis      = n_ana,
    events_OS       = os$events,
    event_ratio     = os$ratio,
    events_PFI      = pfi$events,
    event_ratio_PFI = pfi$ratio,
    stringsAsFactors = FALSE
  ))
}

write.csv(audit, out_csv, row.names = FALSE)
cat("=== AUDIT COMPLETE ===\n")
cat("Saved to:", normalizePath(out_csv), "\n\n")
print(audit)