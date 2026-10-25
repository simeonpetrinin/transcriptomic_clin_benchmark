# EXTERNAL DATA  ->  reproducible loaders for independent validation cohorts
# Output contract for every loader:
#   list(clinical = data.frame(patient_id, OS, OS.time, ...),
#        expr     = numeric matrix [patients x HGNC gene symbols],
#        metadata = list(...))

if (!exists("DIR_RAW")) DIR_RAW <- file.path("data", "raw")
if (!exists("DIR_PROCESSED")) DIR_PROCESSED <- file.path("data", "processed")
if (!exists("DIR_RESULTS")) DIR_RESULTS <- "results"

DIR_EXTERNAL_RAW       <- file.path(DIR_RAW, "external")
DIR_EXTERNAL_PROCESSED <- file.path(DIR_PROCESSED, "external")
DIR_EXTERNAL_RESULTS   <- file.path(DIR_RESULTS, "external")
DIR_EXTERNAL_CHECKPOINTS <- file.path(
  if (exists("DIR_CHECKPOINTS")) DIR_CHECKPOINTS else file.path("data", "checkpoints"),
  "external"
)

for (d in c(DIR_EXTERNAL_RAW, DIR_EXTERNAL_PROCESSED, DIR_EXTERNAL_RESULTS,
            DIR_EXTERNAL_CHECKPOINTS))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

path_external_modelling <- function(cohort_id)
  file.path(DIR_EXTERNAL_PROCESSED, paste0(cohort_id, "_modelling.rds"))
path_external_results <- function(cohort_id)
  file.path(DIR_EXTERNAL_RESULTS, paste0(cohort_id, "_results.csv"))
path_external_predictions <- function(cohort_id)
  file.path(DIR_EXTERNAL_RESULTS, paste0(cohort_id, "_preds.rds"))
path_external_variable_audit <- function(cohort_id)
  file.path(DIR_EXTERNAL_RESULTS, paste0(cohort_id, "_variable_audit.csv"))
path_external_genes <- function(cohort_id)
  file.path(DIR_EXTERNAL_RESULTS, paste0(cohort_id, "_genes.csv.gz"))
path_external_manifest <- function(cohort_id)
  file.path(DIR_EXTERNAL_RESULTS, paste0(cohort_id, "_manifest.txt"))
path_external_subsets <- function(cohort_id)
  file.path(DIR_EXTERNAL_PROCESSED, paste0(cohort_id, "_subsets.csv"))
path_external_repeat_checkpoint <- function(cohort_id, r)
  file.path(DIR_EXTERNAL_CHECKPOINTS, sprintf("%s_repeat_%02d.rds", cohort_id, r))
path_external_settings <- function(cohort_id)
  file.path(DIR_EXTERNAL_CHECKPOINTS, paste0(cohort_id, "_settings.txt"))
path_external_design_guard <- function(cohort_id)
  file.path(DIR_EXTERNAL_CHECKPOINTS, paste0(cohort_id, "_parallel_design.txt"))
path_external_level_fixes <- function(cohort_id)
  file.path(DIR_EXTERNAL_RESULTS, paste0(cohort_id, "_level_fixes.csv"))

DAYS_PER_MONTH <- 365.25 / 12

# Only cohorts that can support a clinically meaningful shared-variable lattice
# are declared. Variables are re-audited after loading; a variable below its
# declared completeness gate or constant variable is excluded loudly in
# *_variable_audit.csv. Completeness gate 50%, as for TCGA.
# Selection: same inclusion criteria as TCGA (primary tumours untreated at
# sampling, genome-wide expression, overall survival, core clinical
# variables); no minimum size.
EXTERNAL_COHORTS <- list(
  GSE68465 = list(
    cancer = "LUAD", kind = "GEO", accession = "GSE68465", platform = "GPL96",
    log2 = TRUE, expected_n = 442L,
    variables = c("age", "sex", "race", "T", "N", "smoking", "grade"),
    loader = "load_gse68465_external",
    citation = "Shedden et al., Nat Med 2008; GEO GSE68465"
  ),
  GSE39582 = list(
    cancer = "COAD", kind = "GEO", accession = "GSE39582", platform = "GPL570",
    log2 = FALSE, expected_n = 556L,
    variables = c("age", "sex", "stage", "T", "N", "M"),
    loader = "load_gse39582_external",
    citation = "Marisa et al., PLoS Med 2013; GEO GSE39582"
  ),
  GSE65858 = list(
    cancer = "HNSC", kind = "GEO", accession = "GSE65858", platform = "GPL10558",
    log2 = FALSE, expected_n = 253L,  #primary tumours only; 14 secondary and 3 relapse dropped
    variables = c("age", "sex", "stage", "T", "N", "smoking", "hn_site", "packyears", "alcohol", "hpv"),
    loader = "load_gse65858_external",
    citation = "Wichmann et al., Int J Cancer 2015; GEO GSE65858"
  ),
  GSE13507 = list(
    cancer = "BLCA", kind = "GEO", accession = "GSE13507", platform = "GPL6102",
    log2 = FALSE, expected_n = 165L,
    variables = c("age", "sex", "grade", "T", "N", "M", "invasiveness"),
    variable_max_top_share = c(M = 1),  #M kept as an exception to the 95% rule (8/165 M1)
    loader = "load_gse13507_external",
    citation = "Kim et al., Mol Cancer 2010; GEO GSE13507"
  ),
  GSE84437 = list(
    cancer = "STAD", kind = "GEO", accession = "GSE84437", platform = "GPL6947",
    log2 = TRUE, expected_n = NA_integer_,  #GenomeStudio quantile-normalised values are unlogged
    variables = c("age", "sex", "T", "N"),
    loader = "load_gse84437_external",
    citation = "Yoon et al., Clin Cancer Res 2020; GEO GSE84437"
  ),
  # CGGA: two diseases from one resource, split like TCGA (LGG = WHO II-III,
  # GBM = WHO IV). Primary tumours only (recurrent/secondary excluded); children kept.
  # IDH, 1p/19q and MGMT are recorded here but absent from TCGA's clinical
  # files: the lattice therefore shows RNA's lift with and without them.
  # Radio/chemotherapy status is post-diagnosis treatment and is NOT used.
  CGGA_LGG = list(
    cancer = "LGG", kind = "MANUAL", accession = "CGGA mRNAseq_693",
    expected_n = NA_integer_,
    variables = c("age", "sex", "histology", "IDH", "codel_1p19q", "MGMT"),  #grade dropped: histology labels encode it
    loader = "load_cgga_lgg_external",
    citation = "Zhao et al., Genomics Proteomics Bioinformatics 2021 (CGGA)"
  ),
  CGGA_GBM = list(
    cancer = "GBM", kind = "MANUAL", accession = "CGGA mRNAseq_693",
    expected_n = NA_integer_,
    variables = c("age", "sex", "IDH", "MGMT"),
    loader = "load_cgga_gbm_external",
    citation = "Zhao et al., Genomics Proteomics Bioinformatics 2021 (CGGA)"
  ),
  # Beat AML: de novo AML, specimen taken at initial diagnosis, one specimen per
  # patient (bone marrow preferred). ELN 2017 risk is the AML analogue of
  # cytogenetic risk (absent from TCGA-LAML's clinical files).
  BEATAML = list(
    cancer = "LAML", kind = "MANUAL", accession = "Beat AML (waves 1-4)",
    expected_n = NA_integer_,
    variables = c("age", "sex", "race", "ELN", "NPM1", "FLT3_ITD", "blasts_bm", "blasts_pb", "FAB", "who_dx"),  #ELN 2017 + WHO 2016 criteria
    loader = "load_beataml_external",
    citation = "Bottomly et al., Cancer Cell 2022; Tyner et al., Nature 2018"
  ),
  #further cohorts
  GSE14520 = list(
    cancer = "LIHC", kind = "GEO_SUPPL", accession = "GSE14520", platform = "GPL3921",
    log2 = FALSE, expected_n = NA_integer_,
    clinical_file = "GSE14520_Extra_Supplement.txt.gz",
    variables = c("age", "sex", "stage", "bclc", "afp_high", "tumor_size_gt5", "multinodular", "cirrhosis",
                  "hbv", "alt_high", "clip"),
    loader = "load_gse14520_external",
    citation = "Roessler et al., Cancer Res 2010; GEO GSE14520"
  ),
  #ER = IHC only, HER2 = SNP6 copy number: the status fields adjusted with RNA expression are not used
  METABRIC = list(
    cancer = "BRCA", kind = "CBIOPORTAL", accession = "brca_metabric",
    expected_n = NA_integer_,
    variables = c("age", "histology", "ER", "HER2", "grade", "tumor_size", "positive_nodes", "NPI", "stage"),
    loader = "load_metabric_external",
    citation = "Curtis et al., Nature 2012; Pereira et al., Nat Commun 2016"
  )
)

.geo_bucket <- function(accession, prefix) {
  digits <- sub(paste0("^", prefix), "", accession)
  stem <- if (nchar(digits) <= 3L) "" else substr(digits, 1L, nchar(digits) - 3L)
  paste0(prefix, stem, "nnn")
}

.geo_matrix_url <- function(accession) sprintf(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/%s_series_matrix.txt.gz",
  .geo_bucket(accession, "GSE"), accession, accession
)

.geo_annot_url <- function(platform) sprintf(
  "https://ftp.ncbi.nlm.nih.gov/geo/platforms/%s/%s/annot/%s.annot.gz",
  .geo_bucket(platform, "GPL"), platform, platform
)

.download_once <- function(url, dest, force = FALSE) {
  if (file.exists(dest) && file.info(dest)$size > 0 && !force) return(dest)
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(dest, ".part")
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  download.file(url, tmp, mode = "wb", method = "libcurl", quiet = FALSE)
  if (!file.exists(tmp) || file.info(tmp)$size == 0) stop("Download produced an empty file: ", url)
  if (!file.rename(tmp, dest)) stop("Could not move completed download to ", dest)
  dest
}

download_geo_inputs <- function(cohort_id, force = FALSE) {
  sp <- EXTERNAL_COHORTS[[cohort_id]]
  if (is.null(sp) || sp$kind != "GEO") stop("Not a registered GEO cohort: ", cohort_id)
  d <- file.path(DIR_EXTERNAL_RAW, cohort_id)
  matrix_path <- file.path(d, paste0(sp$accession, "_series_matrix.txt.gz"))
  annot_path  <- file.path(d, paste0(sp$platform, ".annot.gz"))
  .download_once(.geo_matrix_url(sp$accession), matrix_path, force)
  .download_once(.geo_annot_url(sp$platform), annot_path, force)
  list(matrix = matrix_path, annotation = annot_path)
}

.unquote_geo <- function(x) {
  x <- sub('^"', "", x)
  sub('"$', "", x)
}

.normalise_key <- function(x) {
  x <- tolower(trimws(x))
  x <- gsub("[^a-z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}

.read_geo_header <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  lines <- character(0)
  repeat {
    z <- readLines(con, n = 1L, warn = FALSE)
    if (!length(z) || startsWith(z, "!series_matrix_table_begin")) break
    lines <- c(lines, z)
  }
  sample_lines <- lines[startsWith(lines, "!Sample_")]
  rows <- strsplit(sample_lines, "\t", fixed = TRUE)
  labels <- vapply(rows, `[`, character(1), 1L)
  get_row <- function(label) {
    hit <- which(labels == label)
    if (!length(hit)) return(NULL)
    .unquote_geo(rows[[hit[1L]]][-1L])
  }
  ids <- get_row("!Sample_geo_accession")
  if (is.null(ids)) stop("GEO matrix has no !Sample_geo_accession row: ", path)
  p <- data.frame(sample_id = ids, stringsAsFactors = FALSE)
  for (label in c("!Sample_title", "!Sample_source_name_ch1", "!Sample_platform_id")) {
    z <- get_row(label)
    if (!is.null(z)) p[[.normalise_key(sub("^!Sample_", "", label))]] <- z
  }
  for (r in rows[labels == "!Sample_characteristics_ch1"]) {
    values <- .unquote_geo(r[-1L])
    for (i in seq_along(values)) {
      pos <- regexpr(":", values[i], fixed = TRUE)[1L]
      if (is.na(pos) || pos < 2L) next
      key <- .normalise_key(substr(values[i], 1L, pos - 1L))
      val <- trimws(substr(values[i], pos + 1L, nchar(values[i])))
      if (!nzchar(key)) next
      if (!key %in% names(p)) p[[key]] <- NA_character_
      p[[key]][i] <- val
    }
  }
  p
}

.geo_table_skip <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  n <- 0L
  repeat {
    z <- readLines(con, n = 1L, warn = FALSE)
    if (!length(z)) stop("No !series_matrix_table_begin marker in ", path)
    n <- n + 1L
    if (startsWith(z, "!series_matrix_table_begin")) return(n)
  }
}

.find_column <- function(nms, candidates) {
  hit <- match(tolower(candidates), tolower(nms), nomatch = 0L)
  hit <- hit[hit > 0L]
  if (!length(hit)) return(NA_integer_)
  hit[1L]
}

.read_platform_annotation <- function(path) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required")
  # GEO annotation tables end with !platform_table_end.  fread correctly
  # discards that non-tabular footer but emits a warning for some platforms
  # (notably GPL6102).  Muffle only that known GEO-footer warning; preserve all
  # other parser warnings.
  a <- withCallingHandlers(
    data.table::fread(path, sep = "\t", header = TRUE, data.table = FALSE,
                      check.names = FALSE, showProgress = FALSE),
    warning = function(w) {
      msg <- conditionMessage(w)
      if (grepl("Discarded single-line footer", msg, fixed = TRUE) &&
          grepl("!platform_table_end", msg, fixed = TRUE))
        invokeRestart("muffleWarning")
    }
  )
  id_i <- .find_column(names(a), c("ID", "ID_REF", "ProbeID"))
  sy_i <- .find_column(names(a), c("Gene symbol", "Gene Symbol", "GENE_SYMBOL"))
  if (is.na(id_i) || is.na(sy_i))
    stop("Platform annotation must contain probe ID and Gene symbol columns: ", path)
  symbols <- trimws(as.character(a[[sy_i]]))
  symbols <- sub("[[:space:]]*///.*$", "", symbols)
  symbols <- sub("[[:space:]]*//.*$", "", symbols)
  symbols <- toupper(symbols)
  symbols[symbols %in% c("", "---", "NA", "N/A")] <- NA_character_
  data.frame(probe = as.character(a[[id_i]]), gene = symbols, stringsAsFactors = FALSE)
}

collapse_expression_to_genes <- function(expr, symbols) {
  expr <- as.matrix(expr)
  if (length(symbols) != ncol(expr)) stop("One gene symbol is required per expression column")
  symbols <- toupper(trimws(as.character(symbols)))
  keep <- !is.na(symbols) & nzchar(symbols) & symbols != "---"
  expr <- expr[, keep, drop = FALSE]
  symbols <- symbols[keep]
  if (!length(symbols)) stop("No expression columns have usable gene symbols")
  # Mean aggregation is fixed in advance and uses no survival outcome.
  sums <- rowsum(t(expr), group = symbols, reorder = FALSE, na.rm = TRUE)
  # rowsum() rejects logical matrices, so make the observed-value mask numeric
  # before summing the number of non-missing probes mapped to each gene.
  observed <- !is.na(expr)
  storage.mode(observed) <- "integer"
  counts <- rowsum(t(observed), group = symbols, reorder = FALSE)
  z <- sums / pmax(counts, 1)
  z[counts == 0] <- NA_real_
  out <- t(z)
  storage.mode(out) <- "double"
  out
}

.read_geo_expression <- function(matrix_path, annotation_path, log2_transform = FALSE) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required")
  d <- data.table::fread(
    matrix_path, skip = .geo_table_skip(matrix_path), header = TRUE,
    fill = TRUE, data.table = FALSE, check.names = FALSE, showProgress = FALSE
  )
  names(d)[1L] <- "probe"
  d <- d[!is.na(d$probe) & nzchar(d$probe) & d$probe != "!series_matrix_table_end", , drop = FALSE]
  probes <- .unquote_geo(as.character(d$probe))
  sample_ids <- .unquote_geo(names(d)[-1L])
  x <- t(as.matrix(d[, -1L, drop = FALSE]))
  suppressWarnings(storage.mode(x) <- "double")
  rownames(x) <- sample_ids
  colnames(x) <- probes
  if (log2_transform) x <- log2(pmax(x, 0) + 1)
  ann <- .read_platform_annotation(annotation_path)
  gene <- ann$gene[match(colnames(x), ann$probe)]
  collapse_expression_to_genes(x, gene)
}

load_geo_matrix <- function(cohort_id, force_download = FALSE) {
  sp <- EXTERNAL_COHORTS[[cohort_id]]
  paths <- download_geo_inputs(cohort_id, force_download)
  p <- .read_geo_header(paths$matrix)
  x <- .read_geo_expression(paths$matrix, paths$annotation, isTRUE(sp$log2))
  common <- intersect(p$sample_id, rownames(x))
  if (!length(common)) stop("GEO phenotype and expression sample IDs do not overlap")
  p <- p[match(common, p$sample_id), , drop = FALSE]
  x <- x[common, , drop = FALSE]
  list(pheno = p, expr = x, paths = paths)
}

.as_num <- function(x) suppressWarnings(as.numeric(mark_na(x)))

.as_event <- function(x) {
  z <- tolower(trimws(as.character(mark_na(x))))
  out <- rep(NA_integer_, length(z))
  out[z %in% c("1", "true", "yes", "dead", "death", "deceased") |
        grepl("(^|:)deceased$", z)] <- 1L
  out[z %in% c("0", "false", "no", "alive", "living", "survival") |
        grepl("(^|:)living$", z)] <- 0L
  out
}

.as_sex <- function(x) {
  z <- tolower(trimws(as.character(mark_na(x))))
  out <- rep(NA_character_, length(z))
  out[z %in% c("m", "male")] <- "male"
  out[z %in% c("f", "female")] <- "female"
  out
}

.as_race <- function(x) {
  z <- tolower(trimws(as.character(mark_na(x))))
  z[grepl("^white", z)] <- "white"
  z[grepl("black|african", z)] <- "black or african american"
  z[grepl("^asian", z)] <- "asian"
  z
}

.as_binary_smoking <- function(x) {
  z <- tolower(trimws(as.character(mark_na(x))))
  out <- rep(NA_character_, length(z))
  out[grepl("never|non-smok|nonsmok|^no$", z)] <- "never"
  # Do not let the word "smoked" in "never smoked" overwrite "never".
  ever <- grepl("former|past|current|smoked|^yes$", z) & is.na(out)
  out[ever] <- "ever"
  out
}

.extract_t <- function(x) {
  z <- toupper(as.character(x))
  hit <- grepl("T[0-4]", z)
  out <- rep(NA_character_, length(z))
  out[hit] <- sub(".*(T[0-4][A-D]?).*", "\\1", z[hit])
  collapse_T(out)
}

.extract_n <- function(x) {
  z <- toupper(as.character(x))
  hit <- grepl("N[0-3]", z)
  out <- rep(NA_character_, length(z))
  out[hit] <- sub(".*(N[0-3][A-D]?).*", "\\1", z[hit])
  collapse_N(out)
}

#GSE68465 smoking history -> never / former / current, as TCGA (Unknown and "--" -> missing)
.smoking_3 <- function(x) {
  z <- tolower(trimws(as.character(mark_na(x))))
  out <- rep(NA_character_, length(z))
  out[z == "never smoked"] <- "never"; out[z == "smoked in the past"] <- "former"; out[z == "currently smoking"] <- "current"
  out
}

#M from a combined TNM text such as "T2N1M0"
.extract_m <- function(x) {
  z <- toupper(as.character(x))
  ifelse(grepl("M[01X]", z), sub(".*(M[01X]).*", "\\1", z), NA_character_)
}

#grade from words: well / moderate / poor differentiation -> G1 / G2 / G3
.grade_from_words <- function(x) {
  z <- tolower(as.character(mark_na(x)))
  out <- rep(NA_character_, length(z))
  out[grepl("well", z)] <- "G1"; out[grepl("moderate", z)] <- "G2"; out[grepl("poor", z)] <- "G3"
  out
}

.collapse_stage_any <- function(x) {
  z <- toupper(trimws(as.character(mark_na(x))))
  z[z %in% c("1", "1.0")] <- "I"
  z[z %in% c("2", "2.0")] <- "II"
  z[z %in% c("3", "3.0")] <- "III"
  z[z %in% c("4", "4.0")] <- "IV"
  collapse_stage(ifelse(grepl("^STAGE", z), z, paste("Stage", z)))
}

.finalize_external <- function(clinical, expr, cohort_id, metadata = list()) {
  if (anyDuplicated(clinical$patient_id)) stop(cohort_id, ": duplicated patient IDs")
  keep <- !is.na(clinical$OS) & !is.na(clinical$OS.time) & clinical$OS.time > 0
  clinical <- clinical[keep, , drop = FALSE]
  common <- intersect(as.character(clinical$patient_id), rownames(expr))
  clinical <- clinical[match(common, clinical$patient_id), , drop = FALSE]
  expr <- expr[common, , drop = FALSE]
  if (!nrow(clinical)) stop(cohort_id, ": no patients remain with expression and valid OS")
  if (!identical(as.character(clinical$patient_id), rownames(expr)))
    stop(cohort_id, ": clinical/expression alignment failed")
  sp <- EXTERNAL_COHORTS[[cohort_id]]
  if (!is.null(sp$expected_n) && !is.na(sp$expected_n) && nrow(clinical) != sp$expected_n)
    warning(cohort_id, ": expected ", sp$expected_n, " evaluable patients; obtained ", nrow(clinical))
  metadata <- c(list(cohort_id = cohort_id, cancer = sp$cancer,
                     citation = sp$citation, prepared_at = as.character(Sys.time()),
                     expression_unit = "gene-level; probe means; log2 where declared"), metadata)
  structure(list(clinical = clinical, expr = expr, metadata = metadata),
            class = c("external_survival_cohort", "list"))
}

load_gse68465_external <- function(force_download = FALSE) {
  z <- load_geo_matrix("GSE68465", force_download); p <- z$pheno
  keep <- tolower(p$disease_state) == "lung adenocarcinoma"
  p <- p[keep, , drop = FALSE]; x <- z$expr[p$sample_id, , drop = FALSE]
  clinical <- data.frame(
    patient_id = p$sample_id,
    OS = .as_event(p$vital_status),
    OS.time = .as_num(p$months_to_last_contact_or_death) * DAYS_PER_MONTH,
    age = .as_num(p$age), sex = .as_sex(p$sex), race = .as_race(p$race),
    T = .extract_t(p$disease_stage), N = .extract_n(p$disease_stage),
    smoking = .smoking_3(p$smoking_history),  #three levels as TCGA
    grade = .grade_from_words(p$histologic_grade), stringsAsFactors = FALSE
  )
  .finalize_external(clinical, x, "GSE68465", list(platform = "GPL96"))
}

load_gse39582_external <- function(force_download = FALSE) {
  z <- load_geo_matrix("GSE39582", force_download); p <- z$pheno
  keep <- grepl("primary colorectal adenocarcinoma", p$source_name_ch1, ignore.case = TRUE)
  p <- p[keep, , drop = FALSE]; x <- z$expr[p$sample_id, , drop = FALSE]
  clinical <- data.frame(
    patient_id = p$sample_id,
    OS = .as_event(p$os_event), OS.time = .as_num(p$os_delay_months) * DAYS_PER_MONTH,
    age = .as_num(p$age_at_diagnosis_year), sex = .as_sex(p$sex),
    stage = .collapse_stage_any(p$tnm_stage),
    T = collapse_T(p$tnm_t), N = collapse_N(p$tnm_n),
    M = mx_to_M0(collapse_M(p$tnm_m)), stringsAsFactors = FALSE  #MX -> M0, as for TCGA (AJCC 8th)
  )
  .finalize_external(clinical, x, "GSE39582", list(platform = "GPL570"))
}

load_gse65858_external <- function(force_download = FALSE) {
  z <- load_geo_matrix("GSE65858", force_download); p <- z$pheno
  p <- p[which(tolower(trimws(p$tumor_type)) == "primary"), , drop = FALSE]  #secondary and relapse tumours dropped
  site <- p$tumor_site
  # Cavum oris is the Latin label used by this study for oral cavity.
  site[grepl("cavum oris", site, ignore.case = TRUE)] <- "oral cavity"
  clinical <- data.frame(
    patient_id = p$sample_id,
    OS = .as_event(p$os_event), OS.time = .as_num(p$os),
    age = .as_num(p$age), sex = .as_sex(p$gender),
    stage = collapse_stage(paste("Stage", p$uicc_stage)),
    T = collapse_T(paste0("T", p$t_category)),
    N = collapse_N(paste0("N", p$n_category)),
    smoking = .as_binary_smoking(p$smoking),
    hn_site = collapse_hnsc_site(site),
    packyears = .as_num(p$packyears),
    alcohol = tolower(mark_na(p$alcohol)),  #recorded only as ranges: no / 1 to 30 / 31 to 60 / >60
    hpv = tolower(mark_na(p$hpv_dna)), stringsAsFactors = FALSE
  )
  .finalize_external(clinical, z$expr[p$sample_id, , drop = FALSE], "GSE65858", list(platform = "GPL10558"))
}

load_gse13507_external <- function(force_download = FALSE) {
  z <- load_geo_matrix("GSE13507", force_download); p <- z$pheno
  biological_source <- tolower(trimws(as.character(p$biological_source)))
  # Logical data-frame indexing retains an all-NA row wherever `keep` is NA.
  # GSE13507 contains samples without this characteristic, so use which() to
  # retain explicit TRUE values only.
  keep <- !is.na(biological_source) &
    biological_source == "primary bladder cancer"
  p <- p[which(keep), , drop = FALSE]
  if (!nrow(p)) stop("GSE13507: no primary bladder cancer samples were identified")
  
  # Use checked integer matching rather than character matrix subscripting.
  # This gives a useful error if phenotype/expression alignment ever changes
  # and cannot accidentally pass an NA sample identifier to `[.matrix`.
  expr_i <- match(as.character(p$sample_id), rownames(z$expr))
  if (anyNA(expr_i)) {
    missing_ids <- unique(as.character(p$sample_id[is.na(expr_i)]))
    stop("GSE13507: ", length(missing_ids),
         " selected phenotype sample(s) are absent from the expression matrix: ",
         paste(head(missing_ids, 5L), collapse = ", "))
  }
  x <- z$expr[expr_i, , drop = FALSE]
  rownames(x) <- as.character(p$sample_id)
  grade <- toupper(trimws(as.character(mark_na(p$grade))))
  grade[grade == "HIGH"] <- "HIGH GRADE"
  grade[grade == "LOW"]  <- "LOW GRADE"
  clinical <- data.frame(
    patient_id = p$sample_id,
    OS = .as_event(p$overall_survival),
    OS.time = .as_num(p$survival_month) * DAYS_PER_MONTH,
    age = .as_num(p$age), sex = .as_sex(p$sex), grade = grade,
    T = .extract_t(p$stage), N = .extract_n(p$stage),
    M = mx_to_M0(collapse_M(.extract_m(p$stage))),
    invasiveness = tolower(mark_na(p$invasiveness)), stringsAsFactors = FALSE  #nearly duplicates T (one exception)
  )
  .finalize_external(clinical, x, "GSE13507", list(platform = "GPL6102"))
}

.read_cbio_table <- function(path) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required")
  con <- file(path, open = "rt"); on.exit(close(con), add = TRUE)
  first <- readLines(con, n = 20L, warn = FALSE)
  skip <- which(!startsWith(first, "#"))[1L] - 1L
  if (!is.finite(skip)) stop("No cBioPortal table header found: ", path)
  data.table::fread(path, sep = "\t", skip = skip, data.table = FALSE,
                    check.names = FALSE, na.strings = c("", "NA"), showProgress = FALSE)
}

.brca_histology <- function(x) {
  z <- tolower(as.character(mark_na(x)))
  out <- rep("other", length(z)); out[is.na(z)] <- NA_character_
  out[grepl("ductal|nst", z)] <- "ductal"
  out[grepl("lobular", z)] <- "lobular"
  out[grepl("mixed", z)] <- "mixed"
  out
}

# Download the brca_metabric study from cBioPortal/Datahub into data_dir first.
# Required files are named exactly as in the public study archive.
load_metabric_external <- function(
    data_dir = file.path(DIR_EXTERNAL_RAW, "METABRIC", "brca_metabric"), id = "METABRIC"
) {
  patient_file <- file.path(data_dir, "data_clinical_patient.txt")
  sample_file  <- file.path(data_dir, "data_clinical_sample.txt")
  expr_file    <- file.path(data_dir, "data_mrna_illumina_microarray.txt")
  need <- c(patient_file, sample_file, expr_file)
  if (any(!file.exists(need))) stop(
    "METABRIC requires the cBioPortal brca_metabric study files in ", data_dir,
    ". Missing: ", paste(basename(need[!file.exists(need)]), collapse = ", ")
  )
  pat <- .read_cbio_table(patient_file)
  sam <- .read_cbio_table(sample_file)
  meta <- merge(sam, pat, by = "PATIENT_ID", all.x = TRUE, suffixes = c("", ".patient"))
  e <- data.table::fread(expr_file, sep = "\t", data.table = FALSE,
                         check.names = FALSE, showProgress = FALSE)
  symbol_i <- .find_column(names(e), c("Hugo_Symbol", "HUGO_SYMBOL"))
  if (is.na(symbol_i)) stop("METABRIC expression file has no Hugo_Symbol column")
  sample_cols <- which(names(e) %in% meta$SAMPLE_ID)
  if (!length(sample_cols)) stop("METABRIC expression and clinical sample IDs do not overlap")
  x0 <- t(as.matrix(e[, sample_cols, drop = FALSE]))
  suppressWarnings(storage.mode(x0) <- "double")
  rownames(x0) <- names(e)[sample_cols]
  x <- collapse_expression_to_genes(x0, e[[symbol_i]])
  meta <- meta[match(rownames(x), meta$SAMPLE_ID), , drop = FALSE]
  clinical <- data.frame(
    patient_id = meta$SAMPLE_ID,
    
    OS = .as_event(meta$OS_STATUS),
    OS.time = .as_num(meta$OS_MONTHS) * DAYS_PER_MONTH,
    
    age = .as_num(meta$AGE_AT_DIAGNOSIS),
    histology = .brca_histology(meta$HISTOLOGICAL_SUBTYPE),
    
    ER = collapse_ER(sub("^positve$", "positive", tolower(meta$ER_IHC))),  #IHC only ("Positve" is the file's spelling)
    HER2 = .her2_snp6(meta$HER2_SNP6),  #DNA copy number: gain / neutral / loss
    stage = .collapse_stage_any(meta$TUMOR_STAGE),  #stage 0 -> missing, as collapse_stage does for TCGA
    
    grade = ifelse(as.character(mark_na(meta$GRADE)) %in% c("1", "2", "3"), paste0("G", mark_na(meta$GRADE)), NA),
    tumor_size = .as_num(meta$TUMOR_SIZE),
    positive_nodes = .as_num(meta$LYMPH_NODES_EXAMINED_POSITIVE),
    NPI = .as_num(meta$NPI),
    
    stringsAsFactors = FALSE
  )
  .finalize_external(clinical, x, id, list(platform = "Illumina HT-12 v3"))
}

#HER2 by SNP6 copy number; "UNDEF" is a missing code, not a category
.her2_snp6 <- function(x) {
  z <- tolower(trimws(as.character(mark_na(x))))
  z[!z %in% c("gain", "neutral", "loss")] <- NA_character_
  z
}

# LOADERS -- GSE84437, CGGA (LGG/GBM), Beat AML, GSE14520
# .pick_field() stops with the list of available fields whenever an
# expected field is missing, so a name mismatch is a one-line fix.

# first matching field (normalised names); required fields stop loudly
.pick_field <- function(d, candidates, cohort_id, what, required = TRUE) {
  nms <- names(d); key <- .normalise_key(nms)
  for (cnd in candidates) {
    hit <- which(key == .normalise_key(cnd))
    if (length(hit)) return(d[[hit[1L]]])
  }
  for (cnd in candidates[nchar(candidates) >= 3]) {   # then a prefix match (never for
    hit <- which(startsWith(key, .normalise_key(cnd)))  # 1-2 letter names like "t", "n")
    if (length(hit)) return(d[[hit[1L]]])
  }
  if (!required) return(rep(NA, nrow(d)))
  stop(cohort_id, ": no field for ", what, " (tried: ", paste(candidates, collapse = ", "),
       "). Available: ", paste(nms, collapse = ", "))
}

.geo_suppl_url <- function(accession, file) sprintf(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/suppl/%s",
  .geo_bucket(accession, "GSE"), accession, file)

.geo_platform_matrix_url <- function(accession, platform) sprintf(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/%s-%s_series_matrix.txt.gz",
  .geo_bucket(accession, "GSE"), accession, accession, platform)

.log2p1 <- function(x) { x <- as.matrix(x); storage.mode(x) <- "double"; log2(pmax(x, 0) + 1) }

.first_file <- function(dir, pattern, cohort_id, hint) {
  f <- list.files(dir, pattern = pattern, full.names = TRUE, ignore.case = TRUE, recursive = TRUE)
  if (!length(f)) stop(cohort_id, ": no file matching '", pattern, "' in ", dir, ". ", hint)
  f[1L]
}

# read a delimited or zipped table (CGGA ships .zip, Beat AML .txt/.xlsx)
.read_any_table <- function(path) {
  if (grepl("\\.zip$", path, ignore.case = TRUE)) {
    td <- tempfile(); dir.create(td); utils::unzip(path, exdir = td)
    inner <- list.files(td, full.names = TRUE, recursive = TRUE)
    inner <- inner[!grepl("__MACOSX", inner)]
    if (length(inner) != 1L) stop("expected one file inside ", path)
    path <- inner
  }
  if (grepl("\\.xlsx?$", path, ignore.case = TRUE)) {
    if (!requireNamespace("readxl", quietly = TRUE)) stop("readxl is required")
    return(as.data.frame(readxl::read_excel(path, guess_max = 5000)))
  }
  data.table::fread(path, data.table = FALSE, check.names = FALSE, showProgress = FALSE,
                    na.strings = c("", "NA", "NaN", "N/A"))
}

#GSE84437 (GPL6947)
load_gse84437_external <- function(force_download = FALSE) {
  id <- "GSE84437"; z <- load_geo_matrix(id, force_download); p <- z$pheno
  clinical <- data.frame(
    patient_id = p$sample_id,
    OS = .as_event(.pick_field(p, c("death", "os event", "vital status"), id, "OS event")),
    OS.time = .as_num(.pick_field(p, c("duration overall survival", "overall survival"), id, "OS time")) * DAYS_PER_MONTH,
    age = .as_num(.pick_field(p, c("age"), id, "age")),
    sex = .as_sex(.pick_field(p, c("gender", "sex"), id, "sex")),
    T = .extract_t(.pick_field(p, c("ptstage", "pathological t", "t"), id, "T")),  #GEO values are already "T3"
    N = .extract_n(.pick_field(p, c("pnstage", "pathological n", "n"), id, "N")),
    stringsAsFactors = FALSE)
  .finalize_external(clinical, z$expr, id, list(platform = "GPL6947"))
}

#CGGA (manual download from http://www.cgga.org.cn, release 20200506)
# Place in data/raw/external/CGGA/:
#   CGGA.mRNAseq_693.RSEM-genes.20200506.txt(.zip)
#   CGGA.mRNAseq_693_clinical.20200506.txt(.zip)
.load_cgga <- function(id, grades, data_dir = file.path(DIR_EXTERNAL_RAW, "CGGA")) {
  hint <- "Download the mRNAseq_693 expression and clinical files from http://www.cgga.org.cn."
  clin <- .read_any_table(.first_file(data_dir, "mRNAseq_693_clinical", id, hint))
  expr <- .read_any_table(.first_file(data_dir, "mRNAseq_693\\.RSEM-genes", id, hint))
  genes <- toupper(as.character(expr[[1L]]))
  x0 <- t(.log2p1(expr[, -1L, drop = FALSE])); rownames(x0) <- names(expr)[-1L]
  x <- collapse_expression_to_genes(x0, genes)
  pid   <- as.character(.pick_field(clin, c("CGGA_ID"), id, "patient ID"))
  prs   <- tolower(as.character(.pick_field(clin, c("PRS_type"), id, "primary/recurrent")))
  grade <- toupper(gsub("\\s+", "", as.character(.pick_field(clin, c("Grade"), id, "grade"))))
  age   <- .as_num(.pick_field(clin, c("Age"), id, "age"))
  keep  <- prs == "primary" & grade %in% grades  #all ages; a missing age is imputed in-fold like any rare missing value
  keep[is.na(keep)] <- FALSE
  clin <- clin[keep, , drop = FALSE]; pid <- pid[keep]; grade <- grade[keep]
  ev <- .pick_field(clin, c("Censor (alive=0; dead=1)", "Censor"), id, "OS event")
  idh <- tolower(as.character(mark_na(.pick_field(clin, c("IDH_mutation_status"), id, "IDH"))))
  cod <- tolower(as.character(mark_na(.pick_field(clin, c("1p19q_codeletion_status", "1p19q"), id, "1p/19q"))))
  mgm <- tolower(as.character(mark_na(.pick_field(clin, c("MGMTp_methylation_status", "MGMT"), id, "MGMT"))))
  clinical <- data.frame(
    patient_id = pid,
    OS = as.integer(.as_num(ev)),
    OS.time = .as_num(.pick_field(clin, c("OS"), id, "OS time")),          # days
    age = age[keep], sex = .as_sex(.pick_field(clin, c("Gender"), id, "sex")),
    grade = grade, histology = tolower(as.character(mark_na(.pick_field(clin, c("Histology"), id, "histology")))),
    IDH = idh, codel_1p19q = cod, MGMT = mgm, stringsAsFactors = FALSE)
  .finalize_external(clinical, x, id, list(platform = "CGGA mRNAseq_693 (RSEM, log2+1)"))
}
load_cgga_lgg_external <- function(...) .load_cgga("CGGA_LGG", c("WHOII", "WHOIII"), ...)
load_cgga_gbm_external <- function(...) .load_cgga("CGGA_GBM", c("WHOIV"), ...)

#Beat AML (manual download: https://biodev.github.io/BeatAML2/)
# Place in data/raw/external/BEATAML/ the waves 1-4 clinical table
# (beataml_wv1to4_clinical.xlsx) and normalised expression
# (beataml_waves1to4_norm_exp_dbgap.txt).
load_beataml_external <- function(data_dir = file.path(DIR_EXTERNAL_RAW, "BEATAML")) {
  id <- "BEATAML"; hint <- "Download the waves 1-4 clinical and normalised expression files from https://biodev.github.io/BeatAML2/."
  clin <- .read_any_table(.first_file(data_dir, "clinical", id, hint))
  expr <- .read_any_table(.first_file(data_dir, "norm_exp", id, hint))
  sym_i <- .find_column(names(expr), c("display_label", "gene_symbol", "symbol"))
  if (is.na(sym_i)) stop(id, ": expression file has no gene-symbol column. Columns: ",
                         paste(head(names(expr), 10), collapse = ", "))
  samp_cols <- setdiff(names(expr)[vapply(expr, is.numeric, logical(1))], names(expr)[sym_i])
  x0 <- t(as.matrix(expr[, samp_cols, drop = FALSE])); rownames(x0) <- samp_cols
  x <- collapse_expression_to_genes(x0, expr[[sym_i]])        # already normalised (log2 scale)
  rna_id  <- as.character(.pick_field(clin, c("dbgap_rnaseq_sample"), id, "RNA-seq sample ID"))
  subj    <- as.character(.pick_field(clin, c("dbgap_subject_id"), id, "subject ID"))
  denovo  <- tolower(as.character(.pick_field(clin, c("isDenovo"), id, "de novo flag")))
  stage   <- tolower(as.character(.pick_field(clin, c("diseaseStageAtSpecimenCollection"), id, "specimen timing")))
  spec    <- tolower(as.character(.pick_field(clin, c("specimenType"), id, "specimen type")))
  keep <- !is.na(rna_id) & rna_id %in% rownames(x) & denovo == "true" & grepl("initial", stage)
  keep[is.na(keep)] <- FALSE
  d <- clin[keep, , drop = FALSE]; rna_id <- rna_id[keep]; subj <- subj[keep]; spec <- spec[keep]
  # one specimen per patient: bone marrow before peripheral blood, then first listed
  ord <- order(subj, !grepl("marrow", spec)); d <- d[ord, , drop = FALSE]
  rna_id <- rna_id[ord]; subj <- subj[ord]
  first <- !duplicated(subj); d <- d[first, , drop = FALSE]; rna_id <- rna_id[first]
  eln <- tolower(as.character(mark_na(.pick_field(d, c("ELN2017", "ELN2022"), id, "ELN risk"))))
  # combined labels such as "FavorableOrIntermediate" mean the risk could not be
  # resolved (missing data) -> NA, i.e. modelled as missing
  ambiguous <- grepl("(favorable|intermediate|adverse).*(favorable|intermediate|adverse)", eln)
  eln[grepl("favorable", eln)] <- "favorable"; eln[grepl("intermediate", eln)] <- "intermediate"
  eln[grepl("adverse", eln)] <- "adverse"; eln[ambiguous] <- NA
  eln[!eln %in% c("favorable", "intermediate", "adverse")] <- NA
  #FLT3-ITD as ELN 2017 uses it: negative / allelic ratio < 0.5 / >= 0.5 (positive without a ratio -> NA)
  itd <- tolower(as.character(mark_na(.pick_field(d, c("FLT3-ITD"), id, "FLT3-ITD"))))
  ar  <- .as_num(.pick_field(d, c("allelic_ratio"), id, "FLT3-ITD allelic ratio"))
  flt3 <- rep(NA_character_, length(itd))
  flt3[itd %in% "negative"] <- "negative"
  flt3[itd %in% "positive" & !is.na(ar)] <- ifelse(ar[itd %in% "positive" & !is.na(ar)] >= 0.5, "itd_high", "itd_low")
  #FAB: sub-labels to their class (M5a -> M5, M4eo -> M4); ambiguous labels such as "M0/M1" -> NA
  fab <- toupper(gsub("\\s+", "", as.character(mark_na(.pick_field(d, c("fabBlastMorphology"), id, "FAB")))))
  fab <- sub("^(M[0-7])[A-Z]+$", "\\1", fab)
  fab[!grepl("^M[0-7]$", fab)] <- NA_character_
  clinical <- data.frame(
    patient_id = rna_id,
    OS = .as_event(.pick_field(d, c("vitalStatus"), id, "vital status")),
    OS.time = .as_num(.pick_field(d, c("overallSurvival"), id, "OS time")),   # days
    age = .as_num(.pick_field(d, c("ageAtDiagnosis"), id, "age")),
    sex = .as_sex(.pick_field(d, c("consensus_sex", "reportedSex"), id, "sex")),
    race = .as_race(.pick_field(d, c("reportedRace", "race"), id, "race", required = FALSE)),
    ELN = eln,
    NPM1 = tolower(as.character(mark_na(.pick_field(d, c("NPM1"), id, "NPM1")))),
    FLT3_ITD = flt3,
    blasts_bm = .as_num(.pick_field(d, c("%.Blasts.in.BM"), id, "marrow blasts")),
    blasts_pb = .as_num(.pick_field(d, c("%.Blasts.in.PB"), id, "blood blasts")),
    FAB = fab,
    who_dx = trimws(as.character(mark_na(.pick_field(d, c("specificDxAtInclusion"), id, "WHO diagnosis")))),
    stringsAsFactors = FALSE)
  .finalize_external(clinical, x[rna_id, , drop = FALSE], id,
                     list(platform = "Beat AML RNA-seq (normalised, as deposited)"))
}

#GSE14520 (HBV-related HCC; GPL3921 series; clinical supplement)
load_gse14520_external <- function(force_download = FALSE) {
  id <- "GSE14520"; sp <- EXTERNAL_COHORTS[[id]]; d <- file.path(DIR_EXTERNAL_RAW, id)
  mat <- .download_once(.geo_platform_matrix_url(sp$accession, sp$platform),
                        file.path(d, sprintf("%s-%s_series_matrix.txt.gz", sp$accession, sp$platform)), force_download)
  ann <- .download_once(.geo_annot_url(sp$platform), file.path(d, paste0(sp$platform, ".annot.gz")), force_download)
  sup <- .download_once(.geo_suppl_url(sp$accession, sp$clinical_file), file.path(d, sp$clinical_file), force_download)
  x  <- .read_geo_expression(mat, ann, isTRUE(sp$log2))
  cl <- .read_any_table(sup)
  gsm  <- as.character(.pick_field(cl, c("Affy_GSM", "GSM"), id, "GEO sample ID"))
  tiss <- tolower(as.character(.pick_field(cl, c("Tissue Type"), id, "tissue type")))
  keep <- grepl("tumor", tiss) & !grepl("non", tiss) & gsm %in% rownames(x)
  cl <- cl[keep, , drop = FALSE]; gsm <- gsm[keep]
  yn <- function(v) { z <- tolower(as.character(mark_na(v))); z[z %in% c("n", "no", "0")] <- "no"
  z[z %in% c("y", "yes", "1")] <- "yes"; z[!z %in% c("no", "yes")] <- NA; z }
  afp <- tolower(as.character(mark_na(.pick_field(cl, c("AFP (>300ng/ml)", "AFP"), id, "AFP"))))
  afp[grepl("high|>", afp)] <- "yes"; afp[grepl("low|<", afp)] <- "no"
  clinical <- data.frame(
    patient_id = gsm,
    OS = .as_event(.pick_field(cl, c("Survival status", "survival_status"), id, "OS event")),
    OS.time = .as_num(.pick_field(cl, c("Survival months", "survival_months"), id, "OS time")) * DAYS_PER_MONTH,
    age = .as_num(.pick_field(cl, c("Age"), id, "age")),
    sex = .as_sex(.pick_field(cl, c("Gender"), id, "sex")),
    stage = .collapse_stage_any(.pick_field(cl, c("TNM staging", "TNM"), id, "TNM stage")),
    bclc = toupper(trimws(as.character(mark_na(.pick_field(cl, c("BCLC staging", "BCLC"), id, "BCLC"))))),
    afp_high = afp,
    tumor_size_gt5 = tolower(as.character(mark_na(.pick_field(cl, c("Main Tumor Size (>5 cm)", "Main Tumor Size"), id, "tumour size")))),
    multinodular = yn(.pick_field(cl, c("Multinodular"), id, "multinodular")),
    cirrhosis = yn(.pick_field(cl, c("Cirrhosis"), id, "cirrhosis")),
    hbv = tolower(as.character(mark_na(.pick_field(cl, c("HBV viral status"), id, "HBV status")))),  #CC / AVR-CC / N
    alt_high = c(high = "yes", low = "no")[tolower(as.character(mark_na(.pick_field(cl, c("ALT(>/<=50U/L)", "ALT"), id, "ALT"))))],
    clip = .as_num(.pick_field(cl, c("CLIP staging", "CLIP"), id, "CLIP score")),  #numeric 0-5; "." -> missing
    stringsAsFactors = FALSE)
  .finalize_external(clinical, x, id, list(platform = sp$platform))
}

validate_external_object <- function(dat, cohort_id = "external") {
  if (!is.list(dat) || !is.data.frame(dat$clinical) || !is.matrix(dat$expr))
    stop(cohort_id, ": expected list(clinical=data.frame, expr=matrix)")
  req <- c("patient_id", "OS", "OS.time")
  if (!all(req %in% names(dat$clinical))) stop(cohort_id, ": missing clinical fields: ",
                                               paste(setdiff(req, names(dat$clinical)), collapse = ", "))
  if (is.null(rownames(dat$expr)) || is.null(colnames(dat$expr)))
    stop(cohort_id, ": expression needs patient row names and gene-symbol column names")
  if (!identical(as.character(dat$clinical$patient_id), rownames(dat$expr)))
    stop(cohort_id, ": clinical and expression rows are not identically aligned")
  if (anyDuplicated(dat$clinical$patient_id) || anyDuplicated(colnames(dat$expr)))
    stop(cohort_id, ": duplicated patient IDs or gene symbols")
  invisible(TRUE)
}

prepare_external_cohort <- function(cohort_id, force = FALSE, ...) {
  sp <- EXTERNAL_COHORTS[[cohort_id]]
  if (is.null(sp)) stop("Unknown external cohort: ", cohort_id)
  out <- path_external_modelling(cohort_id)
  if (file.exists(out) && !force) {
    cat(sprintf("[external prepare] %-10s skipped (%s exists)\n", cohort_id, basename(out)))
    return(invisible(readRDS(out)))
  }
  loader <- get(sp$loader, mode = "function")
  dat <- loader(...)
  validate_external_object(dat, cohort_id)
  saveRDS(dat, out, compress = "xz")
  cat(sprintf("[external prepare] %-10s saved: %d patients, %d genes, %d deaths\n",
              cohort_id, nrow(dat$clinical), ncol(dat$expr), sum(dat$clinical$OS == 1L)))
  invisible(dat)
}

prepare_external_panel <- function(cohort_ids = names(EXTERNAL_COHORTS), force = FALSE) {
  failed <- character(0)
  for (id in cohort_ids) {
    ok <- tryCatch({ prepare_external_cohort(id, force = force); TRUE },
                   error = function(e) { message("[external prepare FAILED] ", id, ": ", conditionMessage(e)); FALSE })
    if (!ok) failed <- c(failed, id)
  }
  invisible(list(failed = failed, prepared = setdiff(cohort_ids, failed)))
}