#Load and join
#Builds one cancer's modelling data: a cleaned clinical table (CDR first, PanCanAtlas
#for the rest) and its RNA matrix, for patients with RNA and valid survival. Which
#variables a cancer gets is set in enumerate.R; this file says how to fetch and recode them.

for (p in c("readxl", "data.table")) {
  if (!requireNamespace(p, quietly = TRUE)) stop(sprintf("Package '%s' is missing. Restore project dependencies with renv::restore().", p))
}
suppressPackageStartupMessages({ library(readxl); library(data.table) })

MIN_LEVEL_N <- 6  #rare-level rule: levels below this are pooled (see header)

#missing-value codes (shared with discover_clinical.R)
SENT <- c("", "na", "n/a", "#n/a", "null",
          "not available", "[not available]",
          "not applicable", "[not applicable]",
          "unknown", "[unknown]",
          "not evaluated", "[not evaluated]",
          "not reported", "[not reported]",
          "not collected", "[not collected]",
          "not permitted", "[not permitted]",
          "pending", "[pending]",
          "discrepancy", "[discrepancy]", "--")

#repair invalid UTF-8 text so string matching does not fail
repair_utf8_vector <- function(x) {
  z <- as.character(x)
  bad <- !is.na(z) & !validUTF8(z)
  if (!any(bad)) return(z)
  repaired <- iconv(z[bad], from = "UTF-8", to = "UTF-8", sub = "byte")
  if (any(is.na(repaired)) || any(!validUTF8(repaired)))
    stop("UTF-8 repair failed to produce valid UTF-8.")
  z[bad] <- repaired
  z
}

#turn every missing-value code into NA
mark_na <- function(x) {
  x <- trimws(repair_utf8_vector(x))
  low <- tolower(x)
  x[is.na(x) | low %in% SENT] <- NA
  x
}
lc <- function(x) tolower(mark_na(x))

#A. universal recodes. NX kept, MX -> M0, TX and GX -> missing, T0 -> missing (not SKCM)

#T: T1a/T2b -> T1/T2
collapse_T <- function(x, keep_T0 = FALSE) {
  x <- toupper(trimws(as.character(mark_na(x))))
  out <- sub("^(T[0-4]).*", "\\1", x)  #T1a/T2b -> T1/T2
  out[!grepl("^T[0-4]$", out)] <- NA  #TX, Tis and the rest -> NA
  if (!keep_T0) out[out == "T0"] <- NA
  out
}
#N: N1a etc -> N1
collapse_N <- function(x) {
  x <- toupper(trimws(as.character(mark_na(x))))
  out <- sub("^(N[0-3]).*", "\\1", x)
  out[grepl("^NX", x)] <- "NX"  #NX kept as its own level
  out[!grepl("^(N[0-3]|NX)$", out)] <- NA
  out[out %in% c("N2", "N3")] <- "N2+"
  out
}
#M: M1a/M1b/M1c -> M1
collapse_M <- function(x) {
  x <- toupper(trimws(as.character(mark_na(x))))
  out <- sub("^C?(M[01X]).*", "\\1", x)
  out[!out %in% c("M0", "M1", "MX")] <- NA
  out
}
mx_to_M0 <- function(v) { v[!is.na(v) & v == "MX"] <- "M0"; v }
#stage: "Stage IIA" etc -> I/II/III/IV
collapse_stage <- function(x) {
  x <- toupper(trimws(as.character(mark_na(x))))
  out <- sub("^STAGE\\s+(IV|III|II|I).*", "\\1", x)
  out[!out %in% c("I", "II", "III", "IV")] <- NA
  out
}
#grade: GX (cannot be assessed) and GB (borderline, not on the G1-G4 scale) -> NA
collapse_grade <- function(x) {
  out <- toupper(trimws(as.character(mark_na(x))))
  out[out %in% c("GX", "GB")] <- NA
  out
}
#smoking: never / former / current.
collapse_smoking <- function(x) {
  low <- lc(x); out <- rep(NA_character_, length(low))
  out[grepl("non-smoker|nonsmoker|never", low)] <- "never"
  out[grepl("reformed", low)]                   <- "former"
  out[grepl("^current smoker", low)]            <- "current"
  out
}
#ER: positive / negative ("indeterminate" = failed assay -> NA)
collapse_ER <- function(x) {
  low <- lc(x); out <- rep(NA_character_, length(low))
  out[low == "positive"] <- "positive"
  out[low == "negative"] <- "negative"  #"indeterminate" -> NA
  out
}
#menopause: pre / peri / post ("indeterminate" -> NA)
collapse_menopause <- function(x) {
  low <- lc(x); out <- rep(NA_character_, length(low))
  out[grepl("^pre",  low)] <- "pre"
  out[grepl("^peri", low)] <- "peri"
  out[grepl("^post", low)] <- "post"  #"indeterminate" -> NA
  out
}

#cancer-specific variables (main set: AJCC 8th tier 1 items)

#PR: same rule as ER.
collapse_PR <- function(x) {
  low <- lc(x); out <- rep(NA_character_, length(low))
  out[low == "positive"] <- "positive"
  out[low == "negative"] <- "negative"  #"indeterminate" -> NA
  out
}

#HER2: Negative 556 / Equivocal 179 / Positive 160 / Indeterminate 12
HER2_EQUIVOCAL_NA <- TRUE
collapse_HER2 <- function(x) {
  low <- lc(x); out <- rep(NA_character_, length(low))
  out[low == "negative"] <- "negative"
  out[low == "positive"] <- "positive"
  if (!HER2_EQUIVOCAL_NA) out[low == "equivocal"] <- "equivocal"
  out  #"indeterminate" -> NA
}

#lymphatic / venous invasion: no / yes
collapse_invasion <- function(x) {
  low <- lc(x); out <- rep(NA_character_, length(low))
  out[grepl("^no",  low)] <- "no"
  out[grepl("^yes", low)] <- "yes"
  out
}

#ESCA tumour location: proximal / mid / distal (Distal 128 / Mid 49 / Proximal 6)
collapse_esoph_loc <- function(x) {
  low <- lc(x); out <- rep(NA_character_, length(low))
  out[grepl("proximal|upper", low)] <- "proximal"
  out[grepl("mid",            low)] <- "mid"
  out[grepl("distal|lower",   low)] <- "distal"
  out
}

#HNSC site: 12 raw levels -> the four anatomical regions of head-and-neck staging.
collapse_hnsc_site <- function(x) {
  low <- lc(x); out <- rep(NA_character_, length(low))
  out[grepl("oral|floor of mouth|buccal|alveolar|palate|lip", low)] <- "oral cavity"
  out[grepl("tonsil|base of tongue|oropharynx",               low)] <- "oropharynx"
  out[grepl("hypopharynx",                                    low)] <- "hypopharynx"
  out[grepl("larynx", low) & !grepl("hypopharynx", low)]            <- "larynx"
  out
}

#OV residual disease (extended set only): R0 none, R1 up to 10 mm, R2+ over 10 mm
collapse_residual <- function(x) {
  low <- lc(x); out <- rep(NA_character_, length(low))
  out[grepl("no macroscopic", low)]      <- "R0"
  out[grepl("^1-10", low)]               <- "R1"
  out[grepl("^11-20|^>20", low)]         <- "R2+"
  out
}

#Child-Pugh: A / B+ (C, 1 patient, merged into B); extended set only
collapse_child_pugh <- function(x) {
  low <- lc(x); out <- rep(NA_character_, length(low))
  out[low == "a"] <- "A"
  out[low %in% c("b", "c")] <- "B+"
  out
}

#SARC necrosis: TCGA's four bands -> FNCLCC necrosis score 0/1/2 (AJCC form 41)
collapse_necrosis <- function(x) {
  low <- lc(x); out <- rep(NA_character_, length(low))
  out[grepl("^0%", low)]                    <- "0"
  out[grepl("^<10%|^moderate", low)]        <- "1"
  out[grepl("^extensive", low)]             <- "2"
  out
}

#SARC body site -> the AJCC 8th soft-tissue sarcoma chapter it is staged in
collapse_sarc_site <- function(x) {
  low <- lc(x); out <- rep(NA_character_, length(low))
  one <- !is.na(low) & !grepl("\\|", low)
  out[one & grepl("^(lower|upper) extremity|^superficial trunk|^chest - chest wall", low)] <- "trunk/extremity"
  out[one & grepl("^retroperitoneum/upper abdominal - retroperitoneum$", low)]            <- "retroperitoneum"
  out[one & grepl("^retroperitoneum/upper abdominal - (kidney|colon|small intestines|gastric)$|^chest - lung/pleura$", low)] <- "visceral"
  out[one & grepl("^gynecological - uterus$", low)]                                      <- "uterus"
  out[one & grepl("^head and neck", low)]                                                <- "head and neck"
  out[!is.na(low) & is.na(out)] <- "other/unassigned"  #recorded, but no single chapter (incl. several sites)
  out
}

#OV one or both ovaries: AJCC 8th T1a vs T1b (form 55); left and right are one ovary
collapse_bilateral <- function(x) {
  low <- lc(x); out <- rep(NA_character_, length(low))
  out[low == "bilateral"]            <- "bilateral"
  out[low %in% c("left", "right")]   <- "unilateral"
  out
}

#B. per-cancer level fixes (Harrell rule: under 6 patients), counted on the modelling cohort
LEVEL_FIXES <- list(
  
  BRCA = list(
    histology = list(
      merge = list("Other, specify" = c("Infiltrating Carcinoma NOS"),
                   "Infiltrating Ductal Carcinoma" = c("Medullary Carcinoma")),
      why = "NOS (n=1) means not otherwise specified, which is TCGA's own 'Other, specify'. Medullary (n=5) into ductal: WHO 2019 calls it no-special-type carcinoma with a medullary pattern (nomenclature). Metaplastic (n=9) kept."),
    race = list(NA_ = c("AMERICAN INDIAN OR ALASKA NATIVE"),
                why = "n=1, nominal with no 'Other' level -> missing.")),
  
  SARC = list(histology = list(
    merge = list(
      "Synovial Sarcoma" = c("Sarcoma; synovial; poorly differentiated",
                             "Synovial Sarcoma - Biphasic",
                             "Synovial Sarcoma - Monophasic"),
      "Undifferentiated Pleomorphic Sarcoma (UPS)" = c(
        "Pleomorphic 'MFH' / Undifferentiated pleomorphic sarcoma",
        "Giant cell 'MFH' / Undifferentiated pleomorphic sarcoma with giant cells")),
    NA_ = c("Desmoid Tumor"),
    why = "synovial variants are one entity (n=10 together) and MFH is the old WHO name for UPS: both are nomenclature merges. Desmoid (n=2) is a different disease with no 'Other' level -> missing. MPNST (n=9) kept."),
    site = list(merge = list("other/unassigned" = c("head and neck")),
                why = "head and neck (n=3) joins 'other/unassigned', the site level for sites without one AJCC chapter. Visceral (n=17) kept.")),
  
  CESC = list(
    histology = list(
      merge = list("Endocervical Type of Adenocarcinoma" =
                     c("Endocervical Adenocarcinoma of the Usual Type")),
      NA_ = c("Adenosquamous", "Endometrioid Adenocarcinoma of Endocervix"),
      why = "'Usual Type' and 'Endocervical Type' are the same entity under two names. Adenosquamous (n=5) and endometrioid (n=3) are distinct subtypes with no 'Other' level -> missing."),
    grade = list(merge = list("G3+" = c("G3", "G4")),
                 why = "G4 (n=1) into G3: ordinal."),
    race = list(NA_ = c("NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER"),
                why = "n=2 -> missing. American Indian (n=7) kept.")),
  
  PAAD = list(
    histology = list(
      merge = list("Pancreas-Adenocarcinoma-Other Subtype" =
                     c("Pancreas-Colloid (mucinous non-cystic) Carcinoma",
                       "Pancreas-Undifferentiated Carcinoma")),
      why = "colloid (n=4) and undifferentiated (n=1) into TCGA's own 'Other Subtype' level."),
    stage = list(merge = list("III+" = c("III", "IV")),
                 why = "III (n=3) and IV (n=4). III's smaller neighbour is IV, giving III+ (n=7)."),
    T     = list(merge = list("T3+" = c("T3", "T4")),
                 why = "T4 (n=3) into T3: ordinal. T1 (n=7) kept."),
    grade = list(merge = list("G3+" = c("G3", "G4")),
                 why = "G4 (n=2) into G3: ordinal."),
    N     = list(NA_ = c("NX"),
                 why = "NX (n=4) is not on the N0-N2 scale, so it has no neighbour -> missing.")),
  
  UCEC = list(
    grade = list(merge = list("G3" = c("G3", "HIGH GRADE")),
                 why = "'High Grade' (n=11) is G3 under another name: nomenclature."),
    race  = list(NA_ = c("AMERICAN INDIAN OR ALASKA NATIVE"),
                 why = "n=4 -> missing. Pacific Islander (n=9) kept.")),
  
  #HNSC: grade G4 (n=7) kept.
  HNSC = list(
    histology = list(NA_ = c("Head & Neck Squamous Cell Carcinoma, Spindle Cell Variant"),
                     why = "n=1, a distinct variant with no 'Other' level -> missing (histology is not in the HNSC main set)."),
    N    = list(NA_ = c("NX"),
                why = "NX (n=4) has no neighbour -> missing."),
    race = list(NA_ = c("AMERICAN INDIAN OR ALASKA NATIVE"),
                why = "n=2 -> missing.")),
  
  LUAD = list(
    location = list(NA_ = c("BRONCHIAL", "OTHER (PLEASE SPECIFY)"),
                    why = "not lobes (Other n=4, Bronchial n=0) -> missing (location is not in the LUAD main set)."),
    race = list(NA_ = c("AMERICAN INDIAN OR ALASKA NATIVE"),
                why = paste("n=1 -> missing. As the alphabetically first level it was",
                            "the reference category, so every other race coefficient",
                            "was estimated against one patient and separated in ~80%",
                            "of cells in a test run. Asian (n=8) kept."))),
  
  #LIHC: M (95% rule) and N (only N0 and NX left) are not in the main set
  LIHC = list(
    stage = list(merge = list("III+" = c("III", "IV")),
                 why = "IV (n=4) into III: ordinal."),
    race  = list(NA_ = c("AMERICAN INDIAN OR ALASKA NATIVE"),
                 why = "n=1 -> missing.")),
  
  COAD = list(race = list(NA_ = c("AMERICAN INDIAN OR ALASKA NATIVE"),
                          why = "n=1 -> missing. Asian (n=11) kept.")),
  
  BLCA = list(
    stage = list(merge = list("I-II" = c("I", "II")),
                 why = "I (n=2) into II. Merging the bottom of an ordinal scale gives 'I-II' ('I or II'), not 'I+'."),
    T     = list(merge = list("T1-T2" = c("T1", "T2")),
                 why = "T1 (n=3) into T2, labelled like stage.")),
  
  KIRP = list(
    T    = list(merge = list("T3+" = c("T3", "T4")),
                why = "T4 (n=2) into T3: ordinal."),
    N    = list(merge = list("N1+" = c("N1", "N2+")),
                why = "N2+ (n=4) into N1: ordinal."),
    race = list(NA_ = c("AMERICAN INDIAN OR ALASKA NATIVE"),
                why = "n=2 -> missing. Asian (n=6) kept.")),
  
  ESCA = list(race = list(NA_ = c("BLACK OR AFRICAN AMERICAN"),
                          why = "n=5 -> missing. T4 (n=6) and proximal location (n=6) kept.")),
  
  OV = list(
    stage = list(merge = list("I-II" = c("I", "II")),
                 why = "I (n=1) into II. Ovarian cancer is diagnosed late; stage III dominates."),
    grade = list(merge = list("G2" = c("G1", "G2"), "G3" = c("G3", "G4")),
                 why = "G1 (n=1) and G4 (n=1) into their ordinal neighbours, giving a G2/G3 binary. (GX and GB are set to missing by the universal grade rule.)"),
    race  = list(NA_ = c("AMERICAN INDIAN OR ALASKA NATIVE",
                         "NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER"),
                 why = "n=2 and n=1 -> missing.")),
  
  STAD = list(race = list(NA_ = c("NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER"),
                          why = "n=1 -> missing.")),
  
  LGG  = list(race = list(NA_ = c("AMERICAN INDIAN OR ALASKA NATIVE"),
                          why = "n=1 -> missing. Asian (n=8) kept.")),
  
  GBM  = list(race = list(NA_ = c("ASIAN"),
                          why = "n=5 -> missing.")),
  
  LAML = list(race = list(NA_ = c("ASIAN"),
                          why = "n=2 -> missing.")),
  
  READ = list(N = list(NA_ = c("NX"),
                       why = "NX (n=2) has no neighbour -> missing. T1 (n=9) kept."))
)

#cancers that keep T0 (SKCM: 23 patients with unknown primary)
KEEP_T0 <- c("SKCM")

#apply one cancer's LEVEL_FIXES and print every change
apply_level_fixes <- function(clinical, code, verbose = TRUE) {
  fx <- LEVEL_FIXES[[code]]
  if (is.null(fx)) return(clinical)
  for (v in names(fx)) {
    if (!v %in% names(clinical)) next
    spec <- fx[[v]]
    if (!is.null(spec$merge)) {
      for (new in names(spec$merge)) {
        hit <- clinical[[v]] %in% spec$merge[[new]]
        if (any(hit)) {
          if (verbose) cat(sprintf("  [merge] %s$%s: %d patient(s) -> '%s'\n",
                                   code, v, sum(hit), new))
          clinical[[v]][hit] <- new
        }
      }
    }
    if (!is.null(spec$NA_)) {
      hit <- clinical[[v]] %in% spec$NA_
      if (any(hit)) {
        if (verbose) cat(sprintf("  [NA]    %s$%s: %d patient(s) {%s} -> NA\n",
                                 code, v, sum(hit), paste(spec$NA_, collapse = ", ")))
        clinical[[v]][hit] <- NA
      }
    }
    if (verbose && !is.null(spec$why))
      cat(sprintf("          reason: %s\n", spec$why))
  }
  clinical
}

#one RNA sample per patient: the lowest sample-type code (primary first)
choose_one_per_patient <- function(samps) {
  pt <- substr(samps, 1, 12); code <- as.integer(substr(samps, 14, 15))
  ord <- order(pt, code); samps[ord][!duplicated(pt[ord])]
}

#CDR priority: any field in both files is read from the CDR

#spec helpers: missing codes -> NA
.num   <- function(x) suppressWarnings(as.numeric(mark_na(x)))
.cat   <- function(x) { y <- toupper(mark_na(x)); y }
.range <- function(v, lo, hi) { v[!is.na(v) & (v < lo | v > hi)] <- NA; v }
spec_cat <- function(col) list(src = "PAN", col = col, fn = function(x, code) .cat(x))
spec_num <- function(col, lo = -Inf, hi = Inf)
  list(src = "PAN", col = col, fn = function(x, code) .range(.num(x), lo, hi))
spec_log <- function(col, lo = 0, hi = Inf)
  list(src = "PAN", col = col, fn = function(x, code) log(.range(.num(x), lo, hi) + 1))
#performance status scored after adjuvant therapy -> NA
PS_TIMING <- list(col = "performance_status_scale_timing", pattern = "post")

#VAR_SPEC: source file, column and cleaning function per variable
VAR_SPEC <- list(
  age        = list(src="CDR", col="age_at_initial_pathologic_diagnosis",
                    fn=function(x, code) suppressWarnings(as.numeric(x))),
  sex        = list(src="CDR", col="gender",     fn=function(x, code) mark_na(x)),
  race       = list(src="CDR", col="race",       fn=function(x, code) mark_na(x)),
  grade      = list(src="CDR", col="histological_grade", fn=function(x, code) collapse_grade(x)),
  histology  = list(src="CDR", col="histological_type",  fn=function(x, code) mark_na(x)),
  stage      = list(src="CDR", col="__stage__",  fn=function(x, code) collapse_stage(x)),
  #T and N: pathologic where recorded, clinical as fallback.
  T          = list(src="PAN", col="pathologic_T", fallback="clinical_T",
                    fn=function(x, code) collapse_T(x, keep_T0 = code %in% KEEP_T0)),
  N          = list(src="PAN", col="pathologic_N", fallback="clinical_N",
                    fallback_over="NX",
                    fn=function(x, code) collapse_N(x)),
  #M: pathologic where recorded, clinical as fallback; leftover MX -> M0
  M          = list(src="PAN", col="pathologic_M", fallback="clinical_M",
                    fallback_over="MX",
                    fn=function(x, code) collapse_M(x),
                    post=function(v) mx_to_M0(v)),
  smoking    = list(src="PAN", col="tobacco_smoking_history",
                    fn=function(x, code) collapse_smoking(x)),
  #BLCA papillary and KIRP type 1/2: the histology subtype field (rule i)
  papillary  = spec_cat("diagnosis_subtype"),
  tumor_type = spec_cat("tumor_type"),
  #BRCA: ER, PR and HER2 are required for AJCC 8th stage grouping (form 48, 5.2-5.4)
  ER         = list(src="PAN", col="breast_carcinoma_estrogen_receptor_status",
                    fn=function(x, code) collapse_ER(x)),
  PR         = list(src="PAN", col="breast_carcinoma_progesterone_receptor_status",
                    fn=function(x, code) collapse_PR(x)),
  HER2       = list(src="PAN", col="lab_proc_her2_neu_immunohistochemistry_receptor_status",
                    fn=function(x, code) collapse_HER2(x)),
  #ESCA: tumour location is required for squamous-cell stage grouping (form 16.1, 5.2)
  esoph_loc  = list(src="PAN", col="esophageal_tumor_cental_location",
                    fn=function(x, code) collapse_esoph_loc(x)),
  #HNSC: the subsite picks the AJCC chapter (oral cavity, larynx, oropharynx, hypopharynx)
  hn_site    = list(src="PAN", col="anatomic_neoplasm_subdivision",
                    fn=function(x, code) collapse_hnsc_site(x)),
  #HNSC: extranodal extension is part of the N definition
  ece        = spec_cat("presence_of_pathological_nodal_extracapsular_spread"),
  #number of positive nodes: part of the N definition
  nodes_pos  = spec_num("number_of_lymphnodes_positive_by_he", 0, 150),
  #UCEC: pelvic (N1, IIIC1) and para-aortic (N2, IIIC2) positive nodes
  pelv_pos   = spec_num("total_pelv_lnp", 0, 100),
  aor_pos    = spec_num("total_aor-lnp", 0, 100),
  #COAD: tumour deposits define N1c
  tumour_deposits = spec_cat("non_nodal_tumor_deposits"),
  #LIHC: vascular invasion is part of the T1/T2 definition
  vascular_inv = spec_cat("vascular_tumor_cell_invasion_type"),
  #LIHC: performance status is a BCLC 2022 stage criterion (Reig et al. 2022).
  ecog       = c(spec_num("eastern_cancer_oncology_group", 0, 4), list(mask = PS_TIMING)),
  #PAAD: T1-T3 are defined by the greatest dimension (cm)
  max_size   = spec_num("maximum_tumor_dimension", 0, 50),
  #SARC: T is defined by the greatest dimension
  tumour_length = spec_num("pathologic_tumor_length", 0, 100),
  necrosis   = list(src="PAN", col="tumor_total_necrosis_percent",
                    fn=function(x, code) collapse_necrosis(x)),
  site       = list(src="PAN", col="tumor_tissue_site",
                    fn=function(x, code) collapse_sarc_site(x)),
  #SKCM: Breslow thickness and ulceration define T
  breslow    = list(src="PAN", col="breslow_depth_value",
                    fn=function(x, code) suppressWarnings(as.numeric(mark_na(x)))),
  ulceration = list(src="PAN", col="melanoma_ulceration_indicator",
                    fn=function(x, code) mark_na(x)),
  #UCEC myometrial invasion (%), the T1a/T1b criterion, kept continuous
  myo_inv    = list(src="PAN", col="pct_tumor_invasion",
                    fn=function(x, code) {
                      v <- suppressWarnings(as.numeric(mark_na(x)))
                      v[!is.na(v) & (v < 0 | v > 100)] <- NA
                      v
                    }),
  #OV: one or both ovaries (T1a vs T1b)
  bilateral  = list(src="PAN", col="anatomic_neoplasm_subdivision",
                    fn=function(x, code) collapse_bilateral(x))
)

#C. extended-set variables (EXT_VARS in extended.R).
collapse_resid_tumor <- function(x) {  #R0/R1/R2 kept
  y <- toupper(mark_na(x)); y[!y %in% c("R0", "R1", "R2")] <- NA; y
}
collapse_barretts <- function(x) {  #YES-USA / YES-UK definitions -> YES
  y <- toupper(mark_na(x)); y[grepl("^YES", y)] <- "YES"; y
}
clean_albumin <- function(x) {
  v <- .num(x); g <- !is.na(v) & v >= 10 & v <= 60; v[g] <- v[g] / 10; .range(v, 1, 6)
}
clean_platelets <- function(x) {
  v <- .num(x); g <- !is.na(v) & v >= 1000; v[g] <- v[g] / 1000; .range(v, 5, 1500)
}
clean_creatinine <- function(x) {
  v <- .num(x); g <- !is.na(v) & v > 20; v[g] <- v[g] / 88.4; .range(v, 0.1, 20)
}

VAR_SPEC_EXT <- list(
  #former main-set variables that are not AJCC tier 1
  menopause  = list(src="CDR", col="menopause_status",
                    fn=function(x, code) collapse_menopause(x)),
  child_pugh = list(src="PAN", col="child_pugh_classification_grade",
                    fn=function(x, code) collapse_child_pugh(x)),
  clark      = list(src="PAN", col="melanoma_clark_level_value",
                    fn=function(x, code) mark_na(x)),
  lymph_inv  = list(src="PAN", col="lymphatic_invasion",
                    fn=function(x, code) collapse_invasion(x)),
  venous_inv = list(src="PAN", col="venous_invasion",
                    fn=function(x, code) collapse_invasion(x)),
  residual   = list(src="PAN", col="tumor_residual_disease",
                    fn=function(x, code) collapse_residual(x)),
  #margins and residual tumour
  resid_tumor     = list(src = "CDR", col = "residual_tumor",  #CDR priority
                         fn = function(x, code) collapse_resid_tumor(x)),
  margin          = list(src = "CDR", col = "margin_status",  #CDR priority
                         fn = function(x, code) .cat(x)),
  #pathology / extent
  lvi             = spec_cat("lymphovascular_invasion_present"),
  perineural      = spec_cat("perineural_invasion_present"),
  synchronous     = spec_cat("synchronous_colon_cancer_present"),
  hist_non_mibc   = spec_cat("hist_of_non_mibc"),
  keratinizing    = spec_cat("keratinizing_squamous_cell_carcinoma_present_indicator"),
  unknown_primary = spec_cat("primary_neoplasm_melanoma_dx"),
  tumour_width    = spec_num("pathologic_tumor_width", 0, 100),
  tumour_depth_cm = spec_num("pathologic_tumor_depth", 0, 100),
  multifocal      = spec_cat("tumor_multifocal"),
  peritoneal_wash = spec_cat("peritoneal_wash"),
  mmr             = spec_cat("loss_expression_of_mismatch_repair_proteins_by_ihc"),
  #location
  location        = spec_cat("anatomic_neoplasm_subdivision"),
  tumor_location  = spec_cat("tumor_location"),
  supratentorial  = spec_cat("supratentorial_localization"),
  #blood values
  cea             = spec_log("preoperative_pretreatment_cea_level", 0, 1e5),
  afp             = spec_log("laboratory_procedure_alpha_fetoprotein_outcome_value", 0, 1e7),
  albumin         = list(src = "PAN", col = "laboratory_procedure_albumin_result_specified_value",
                         fn = function(x, code) clean_albumin(x)),
  plt_count       = list(src = "PAN", col = "lab_procedure_platelet_result_specified_value",
                         fn = function(x, code) clean_platelets(x)),
  creatinine      = list(src = "PAN", col = "hematology_serum_creatinine_laboratory_result_value_in_mg_dl",
                         fn = function(x, code) clean_creatinine(x)),
  prothrombin     = spec_num("laboratory_procedure_prothrombin_time_result_value", 5, 100),
  hb_q            = spec_cat("hemoglobin_result"),
  plt_q           = spec_cat("platelet_qualitative_result"),
  wbc_q           = spec_cat("white_cell_count_result"),
  #performance status and organ function
  kps             = c(spec_num("karnofsky_performance_score", 0, 100), list(mask = PS_TIMING)),
  fibrosis        = spec_cat("liver_fibrosis_ishak_score_category"),
  liver_inflam    = spec_cat("adjacent_hepatic_tissue_inflammation_extent_type"),
  #history, comorbidity, family
  aetiology       = spec_cat("history_hepato_carcinoma_risk_factor"),
  diabetes        = spec_cat("history_of_diabetes"),
  pancreatitis    = spec_cat("history_of_chronic_pancreatitis"),
  reflux          = spec_cat("reflux_history"),
  barretts        = list(src = "PAN", col = "barretts_esophagus",
                         fn = function(x, code) collapse_barretts(x)),
  prior_esoph     = spec_cat("history_of_esophageal_cancer"),
  fam_hist        = spec_cat("family_history_of_cancer"),
  fam_hist_stomach= spec_cat("family_history_of_stomach_cancer"),
  relatives_cancer= spec_num("number_of_first_degree_relatives_with_cancer_diagnosis", 0, 20),
  asthma          = spec_cat("asthma_history"),
  hay_fever       = spec_cat("hay_fever_history"),
  allergy_mould   = spec_cat("mold_or_dust_allergy_history"),
  allergy_food    = spec_cat("food_allergy_history"),
  allergy_animal  = spec_cat("animal_insect_allergy_history"),
  ethnicity       = spec_cat("ethnicity"),
  #symptoms and presentation (LGG)
  presenting_symptom = spec_cat("first_presenting_symptom"),
  symptom_duration   = spec_cat("first_presenting_symptom_longest_duration"),
  seizures        = spec_cat("seizure_history"),
  headache        = spec_cat("headache_history"),
  mental_status   = spec_cat("mental_status_changes"),
  motor           = spec_cat("motor_movement_changes"),
  sensory         = spec_cat("sensory_changes"),
  visual          = spec_cat("visual_changes"),
  preop_antiseizure = spec_cat("preoperative_antiseizure_meds"),
  preop_steroids  = spec_cat("preoperative_corticosteroids"),
  #lifestyle and body size
  pack_years      = spec_num("number_pack_years_smoked", 0, 300),
  alcohol         = spec_cat("alcohol_history_documented"),
  alcohol_amount  = spec_num("amount_of_alcohol_consumption_per_day", 0, 50),
  alcohol_freq    = spec_num("frequency_of_alcohol_consumption", 0, 7),
  alcohol_exposure= spec_cat("alcoholic_exposure_category"),
  log_weight      = spec_log("weight", 20, 300),  #kg
  log_height      = spec_log("height", 100, 230)  #cm
)
#stop if an extended variable reuses a main-set name, then merge the two lists
if (any(names(VAR_SPEC_EXT) %in% names(VAR_SPEC)))
  stop("VAR_SPEC_EXT redefines: ",
       paste(intersect(names(VAR_SPEC_EXT), names(VAR_SPEC)), collapse = ", "))
VAR_SPEC <- c(VAR_SPEC, VAR_SPEC_EXT)

#stop if patient IDs are missing or duplicated
validate_patient_ids <- function(x, source, code) {
  x <- mark_na(x)
  if (anyNA(x))
    stop(source, " contains missing patient IDs for '", code, "'")
  if (anyDuplicated(x)) {
    dup <- unique(x[duplicated(x)])
    stop(source, " contains duplicated patient IDs for '", code, "': ",
         paste(head(dup, 5), collapse = ", "),
         if (length(dup) > 5) ", ..." else "")
  }
  x
}

#load one cancer (vars defaults to the main set; read_expr = FALSE skips the RNA matrix)
load_cancer <- function(code, verbose = TRUE, vars = names(CANCER_VARS[[code]]),
                        read_expr = TRUE) {
  if (is.null(vars)) stop("no CANCER_VARS entry for '", code, "'")
  
  #1. CDR: survival, stage (pathologic, else clinical) and the CDR variables
  cdr <- as.data.frame(read_xlsx(CDR_PATH, sheet = "TCGA-CDR", guess_max = 11000,
                                 .name_repair = "unique_quiet"))
  cdr <- cdr[cdr$type == code, ]
  if (nrow(cdr) == 0) stop("no CDR rows for '", code, "'")
  cdr$bcr_patient_barcode <- validate_patient_ids(
    cdr$bcr_patient_barcode, "CDR", code
  )
  stage_path <- collapse_stage(cdr$ajcc_pathologic_tumor_stage)  #pathologic stage
  stage_clin <- collapse_stage(cdr$clinical_stage)  #else clinical
  stage_path[is.na(stage_path)] <- stage_clin[is.na(stage_path)]
  cdr$`__stage__` <- stage_path
  
  clinical <- data.frame(patient_id = cdr$bcr_patient_barcode,
                         OS      = suppressWarnings(as.integer(cdr$OS)),
                         OS.time = suppressWarnings(as.numeric(cdr$OS.time)),
                         stringsAsFactors = FALSE)
  for (v in vars) {
    sp <- VAR_SPEC[[v]]
    if (is.null(sp)) stop("no VAR_SPEC entry for '", v, "'")
    if (sp$src == "CDR") clinical[[v]] <- sp$fn(cdr[[sp$col]], code)
  }
  
  #2. PanCanAtlas: read the needed columns, clean, apply fallback, post and mask
  pan_vars <- vars[vapply(vars, function(v) VAR_SPEC[[v]]$src == "PAN", logical(1))]
  if (length(pan_vars)) {
    need <- c("bcr_patient_barcode", "acronym",
              unique(unlist(lapply(pan_vars, function(v)
                c(VAR_SPEC[[v]]$col, VAR_SPEC[[v]]$fallback, VAR_SPEC[[v]]$mask$col)))))
    hdr  <- names(fread(PANCAN_PATH, sep = "\t", header = TRUE, nrows = 0))
    absent <- setdiff(need, hdr)
    if (length(absent))
      warning(code, ": PanCanAtlas column(s) not found, read as all-missing: ",
              paste(absent, collapse = ", "), call. = FALSE)
    pan <- fread(PANCAN_PATH, sep = "\t", header = TRUE, select = intersect(need, hdr),
                 na.strings = SENT, data.table = FALSE, showProgress = FALSE)
    pan <- pan[pan$acronym == code, ]
    if (nrow(pan) == 0) stop("no PanCanAtlas rows for '", code, "'")
    pan$bcr_patient_barcode <- validate_patient_ids(
      pan$bcr_patient_barcode, "PanCanAtlas", code
    )
    cp <- data.frame(patient_id = pan$bcr_patient_barcode, stringsAsFactors = FALSE)
    for (col in absent) pan[[col]] <- NA_character_
    for (v in pan_vars) {
      sp  <- VAR_SPEC[[v]]
      val <- sp$fn(pan[[sp$col]], code)
      if (!is.null(sp$fallback)) {  #e.g. M: pathologic, else clinical
        fb   <- sp$fn(pan[[sp$fallback]], code)
        need <- is.na(val)
        #fallback_over: also replace NX / MX when the other route assessed the patient
        if (!is.null(sp$fallback_over))
          need <- need | (!is.na(val) & val %in% sp$fallback_over & !is.na(fb) &
                            !fb %in% sp$fallback_over)
        n_fb <- sum(need & !is.na(fb))
        if (n_fb) message(sprintf("  [%s] %s: %d patient(s) took the %s fallback",
                                  code, v, n_fb, sp$fallback))
        val[need] <- fb[need]
      }
      if (!is.null(sp$post)) val <- sp$post(val)  #e.g. M: leftover MX -> M0
      if (!is.null(sp$mask)) {  #e.g. PS scored after adjuvant therapy -> NA
        mk <- tolower(mark_na(pan[[sp$mask$col]]))
        val[!is.na(mk) & grepl(sp$mask$pattern, mk)] <- NA
      }
      cp[[v]] <- val
    }
    clinical <- merge(clinical, cp, by = "patient_id", all.x = TRUE)
  }
  if (anyDuplicated(clinical$patient_id))
    stop("clinical merge produced duplicated patient IDs for '", code, "'")
  
  #3. expression: tumour samples (codes 01-09), one per patient
  con <- gzfile(RNA_PATH); hdr <- readLines(con, n = 1); close(con)
  all_samps <- strsplit(hdr, "\t")[[1]][-1]
  stype <- suppressWarnings(as.integer(substr(all_samps, 14, 15)))
  keep  <- !is.na(stype) & stype < 10 &
    substr(all_samps, 1, 12) %in% clinical$patient_id
  if (!any(keep)) stop("no tumour RNA for '", code, "'")
  chosen  <- choose_one_per_patient(all_samps[keep])
  if (!read_expr) {
    clinical <- clinical[clinical$patient_id %in% substr(chosen, 1, 12), ]
    n_rna    <- nrow(clinical)
    clinical <- apply_level_fixes(clinical[valid_survival(clinical), ], code, verbose)
    attr(clinical, "n_rna") <- n_rna  #patients with RNA before the survival filter
    return(list(clinical = clinical, expr = NULL))
  }
  expr_dt <- fread(RNA_PATH, select = c("sample", chosen),
                   data.table = FALSE, showProgress = FALSE)
  expr <- t(as.matrix(expr_dt[, -1, drop = FALSE]))
  rownames(expr) <- substr(chosen, 1, 12); colnames(expr) <- expr_dt$sample
  if (anyDuplicated(rownames(expr)))
    stop("RNA selection produced duplicated patient IDs for '", code, "'")
  
  #4. align: keep patients with both clinical data and RNA, in the same order
  common   <- intersect(clinical$patient_id, rownames(expr))
  clinical <- clinical[match(common, clinical$patient_id), ]
  expr     <- expr[common, , drop = FALSE]
  if (anyDuplicated(clinical$patient_id))
    stop("aligned clinical data contain duplicated patient IDs for '", code, "'")
  if (!identical(as.character(clinical$patient_id), rownames(expr)))
    stop("clinical/RNA patient alignment failed for '", code, "'")
  
  #5. survival filter: the modelling cohort
  ok       <- valid_survival(clinical)
  n_rna    <- nrow(clinical)
  clinical <- clinical[ok, ]; expr <- expr[ok, , drop = FALSE]
  
  #6. per-cancer level fixes, counted on the modelling cohort (printed)
  clinical <- apply_level_fixes(clinical, code, verbose)
  attr(clinical, "n_rna") <- n_rna  #patients with RNA before the survival filter
  list(clinical = clinical, expr = expr)
}

#valid survival: OS known and OS.time above 0
valid_survival <- function(clinical)
  !is.na(clinical$OS) & !is.na(clinical$OS.time) & clinical$OS.time > 0

#build one cancer, print the report, save.
build_modelling_data <- function(code, verbose = TRUE,
                                 vars = names(CANCER_VARS[[code]]), key = code) {
  dat  <- load_cancer(code, verbose, vars)
  if (verbose) {
    cl <- dat$clinical
    cat(sprintf("\n=== %s ===\n", code))
    cat(sprintf("patients %d (%d with RNA, %d dropped by the survival filter) | genes %d | events %d (%.0f%%)\n",
                nrow(cl), attr(cl, "n_rna"), attr(cl, "n_rna") - nrow(cl),
                ncol(dat$expr), sum(cl$OS == 1), 100 * mean(cl$OS == 1)))
    if (identical(key, code))
      cat(sprintf("variables (%d) -> %s subsets + anchor\n\n", length(vars),
                  format(2^length(vars) - 1, big.mark = ",")))
    else
      cat(sprintf("variables (%d) for explicit specifications (%s)\n\n", length(vars), key))
    todo <- character(0)
    for (v in vars) {
      x <- cl[[v]]; cmp <- mean(!is.na(x))
      if (is.numeric(x)) {
        cat(sprintf("  %-11s %3.0f%%  numeric [%.0f - %.0f]\n", v, 100*cmp,
                    min(x, na.rm=TRUE), max(x, na.rm=TRUE)))
      } else {
        tb <- table(x, useNA = "no")
        cat(sprintf("  %-11s %3.0f%%  %s\n", v, 100*cmp,
                    paste(sprintf("%s(%d)", names(tb), as.integer(tb)), collapse=" ")))
        thin <- names(tb)[tb < MIN_LEVEL_N]
        if (length(thin)) {
          cat(sprintf("               ** UNDER %d: %s **\n", MIN_LEVEL_N,
                      paste(sprintf("%s(%d)", thin, as.integer(tb[thin])), collapse=", ")))
          todo <- c(todo, paste0(code, ":", v))
        }
      }
    }
    if (length(todo))
      cat(sprintf("\n  ACTION: %s break the rare-level rule -> add to LEVEL_FIXES\n          in R/load_data.R\n",
                  paste(todo, collapse = ", ")))
  }
  saveRDS(dat, path_modelling(key))
  cat("\nSaved ", path_modelling(key), "\n", sep = "")
  invisible(dat)
}

#build every cancer in CANCER_VARS, errors won't stop the run
build_all <- function(codes = names(CANCER_VARS), skip_existing = FALSE) {
  failed <- character(0); built <- character(0)
  for (code in codes) {
    if (skip_existing && file.exists(path_modelling(code))) {
      cat(sprintf("\n%s: exists, skipped\n", code)); next
    }
    cat("\n", strrep("=", 68), "\n", sep = "")
    ok <- tryCatch({ build_modelling_data(code); TRUE },
                   error = function(e) {
                     cat(sprintf("\n!! %s FAILED: %s\n", code, conditionMessage(e))); FALSE })
    if (ok) built <- c(built, code) else failed <- c(failed, code)
  }
  cat("\n", strrep("=", 68), "\n", sep = "")
  cat(sprintf("built %d | failed %d\n", length(built), length(failed)))
  if (length(failed)) cat("failed: ", paste(failed, collapse = ", "), "\n", sep = "")
  invisible(list(built = built, failed = failed))
}