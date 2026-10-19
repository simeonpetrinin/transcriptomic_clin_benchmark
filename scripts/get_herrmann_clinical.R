#Herrmann clinical columns
#Downloads the clinical block of Herrmann et al. (2021) exactly as they coded it, from their OpenML
#datasets (IDs from github.com/HerrMo/multi-omics_benchmark_study, data/datset_ids.txt).
#Their clinical-only Cox uses every column whose name contains "_clinical" (get_blocks2() in their
#ancillary_code_bench.R), already numeric and complete. Only those columns and the patient code are kept.
#Writes data/processed/published/herrmann_clinical_<CODE>.csv (patient_id + hm_* columns).
#Usage, from the project root: source("paths.R"); source("scripts/get_herrmann_clinical.R"); get_herrmann_all()

for (p in c("jsonlite", "data.table"))
  if (!requireNamespace(p, quietly = TRUE)) stop(sprintf("Package '%s' is missing.", p))

DIR_PUBLISHED <- file.path(DIR_PROCESSED, "published")

#OpenML dataset IDs per cancer (each dataset was split into 2-3 parts to fit OpenML's size limit)
HERRMANN_IDS <- list(BLCA = c(42279, 42280), BRCA = c(42323, 42324, 42325), COAD = c(42281, 42282),
                     ESCA = c(42283, 42284), HNSC = c(42285, 42286), KIRC = c(42287, 42288),
                     KIRP = c(42289, 42290), LAML = c(42291, 42292), LGG = c(42293, 42294),
                     LIHC = c(42295, 42296), LUAD = c(42297, 42298), LUSC = c(42299, 42300),
                     OV = c(42301, 42302), PAAD = c(42303, 42304), SARC = c(42309, 42310),
                     SKCM = c(42311, 42312), STAD = c(42313, 42314), UCEC = c(42321, 42322))

#clinical columns per cancer in their Table 2 ("clin."), used to check the download
HERRMANN_N_CLIN <- c(BLCA = 5, BRCA = 8, COAD = 7, ESCA = 6, HNSC = 11, KIRC = 9, KIRP = 6, LAML = 7, LGG = 10,
                     LIHC = 11, LUAD = 9, LUSC = 9, OV = 6, PAAD = 10, SARC = 11, SKCM = 9, STAD = 7, UCEC = 11)

#the download link of one OpenML dataset (its ARFF file)
.oml_url <- function(id) {
  d <- jsonlite::fromJSON(sprintf("https://www.openml.org/api/v1/json/data/%d", id))
  u <- d$data_set_description$url
  if (is.null(u) || !nzchar(u)) stop("OpenML gave no download link for dataset ", id)
  u
}

#ARFF header: attribute names in column order, and the line where the data start
.arff_header <- function(f) {
  con <- file(f, "r"); on.exit(close(con))
  att <- character(0); n <- 0L
  repeat {
    l <- readLines(con, n = 1, warn = FALSE); n <- n + 1L
    if (!length(l)) stop("no @data line in ", f)
    if (grepl("^\\s*@data", l, ignore.case = TRUE)) break
    if (grepl("^\\s*@attribute", l, ignore.case = TRUE)) {
      rest <- sub("^\\s*@attribute\\s+", "", l, ignore.case = TRUE)
      nm <- if (grepl("^['\"]", rest)) sub("^['\"]([^'\"]*)['\"].*", "\\1", rest) else sub("\\s.*$", "", rest)  #quoted or bare name
      att <- c(att, nm)
    }
  }
  list(names = att, data_line = n)
}

#one cancer: find the part that holds the "_clinical" columns, keep them with the patient code
get_herrmann <- function(code, force = FALSE) {
  dir.create(DIR_PUBLISHED, showWarnings = FALSE, recursive = TRUE)
  out <- file.path(DIR_PUBLISHED, sprintf("herrmann_clinical_%s.csv", code))
  if (file.exists(out) && !force) { cat(sprintf("[herrmann] %-5s skipped (exists)\n", code)); return(invisible(out)) }
  for (id in HERRMANN_IDS[[code]]) {
    f <- tempfile(fileext = ".arff"); on.exit(unlink(f), add = TRUE)
    cat(sprintf("[herrmann] %-5s downloading OpenML %d ...\n", code, id))
    download.file(.oml_url(id), f, mode = "wb", quiet = TRUE)
    h <- .arff_header(f)
    clin <- grep("_clinical", h$names)
    if (!length(clin)) { unlink(f); next }  #clinical block is in another part
    id_col <- 1L  #their code drops column 1 as the patient code (bench_experiment.R, task[, -1])
    x <- data.table::fread(f, skip = h$data_line, header = FALSE, sep = ",", quote = "'",
                           select = c(id_col, clin), na.strings = "?", data.table = FALSE)
    names(x) <- c("patient_raw", h$names[clin])
    for (v in names(x)) if (is.character(x[[v]])) x[[v]] <- gsub("^[\"']|[\"']$", "", x[[v]])  #drop ARFF quotes
    pid <- toupper(gsub("\\.", "-", as.character(x$patient_raw)))
    x$patient_id <- substr(pid, 1, 12)  #TCGA patient barcode
    if (!all(grepl("^TCGA-", x$patient_id))) stop(code, ": column 1 does not look like TCGA barcodes: ", paste(head(x$patient_raw, 3), collapse = ", "))
    if (anyDuplicated(x$patient_id)) stop(code, ": duplicated patients in Herrmann's data")
    keep <- h$names[clin]
    if (length(keep) != HERRMANN_N_CLIN[[code]])
      warning(code, ": ", length(keep), " clinical columns found, their Table 2 says ", HERRMANN_N_CLIN[[code]], call. = FALSE)
    y <- data.frame(patient_id = x$patient_id, x[, keep, drop = FALSE], check.names = FALSE)
    names(y)[-1] <- paste0("hm_", make.names(sub("_clinical$", "", keep)))
    write.csv(y, out, row.names = FALSE)
    cat(sprintf("[herrmann] %-5s %d patients, %d clinical columns: %s\n", code, nrow(y), length(keep),
                paste(keep, collapse = ", ")))
    return(invisible(out))
  }
  stop(code, ": no '_clinical' columns found in any part")
}

#all 18 cancers; one failing does not stop the rest
get_herrmann_all <- function(codes = names(HERRMANN_IDS), force = FALSE) {
  failed <- character(0)
  for (code in codes) {
    ok <- tryCatch({ get_herrmann(code, force); TRUE },
                   error = function(e) { cat(sprintf("[FAILED]     %-5s herrmann: %s\n", code, conditionMessage(e))); FALSE })
    if (!ok) failed <- c(failed, code)
  }
  invisible(list(failed = failed))
}