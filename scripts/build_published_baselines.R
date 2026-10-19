#Published baselines
#For each cancer's modelling cohort, builds every published clinical baseline in two versions:
#  original = each paper's own preprocessing, from its own code (columns sb_* / jd_* / hm_*)
#  ours     = the same variables, each fetched and coded by this project's preprocessing (VAR_SPEC in load_data.R,
#             level fixes, missing as its own category); only built when every variable of the set has a
#             counterpart here, otherwise the reason is recorded and only the original version exists
#No variable is dropped in either version. A set whose ours version uses only main-lattice variables is read off
#the lattice; every other set is fitted by run_published_baselines.R
#Original recipes:
#  SurvBoard (Wissel et al. 2025; survboard/R/prep/prepare_tcga_data.R; same recipe in Wissel et al. 2023):
#    CDR age, gender, race, raw pathologic stage, raw clinical stage, histological type; missing = its own category
#  Jardillier et al. 2022 (R_scripts/load_data/load_clinCox_data.R): age, gender (if 2 levels), and grade,
#    pathologic T, N, M each only if >= 90% complete in the cancer, T/N/M cut to two characters ("T2a" -> "T2")
#    Stated deviation: they drop incomplete patients; here missing is its own category so every set is
#    fitted on the same patients and folds as the lattice
#  Herrmann et al. 2021: the "_clinical" columns of their OpenML data, numeric as they coded them
#    (scripts/get_herrmann_clinical.R); patients missing from their data get one shared absent flag
#Writes data/processed/published/<CODE>_published_clinical.csv and data/processed/published/published_specs.csv
#Usage, from the project root: source("master.R"); source("scripts/build_published_baselines.R"); build_published_all()

DIR_PUBLISHED   <- file.path(DIR_PROCESSED, "published")
JD_MIN_COMPLETE <- 0.90  #Jardillier's perc_min_clin = 90 (main.R)

#fields used by Herrmann that VAR_SPEC does not cover yet, coded the project's way (spec_cat / spec_num helpers)
PUB_SPEC <- list(
  pub_surgery_breast  = spec_cat("breast_carcinoma_surgical_procedure_name"),
  pub_neck_dissection = spec_cat("lymphnode_neck_dissection"),
  pub_laterality      = spec_cat("laterality"),
  pub_surgery_type    = spec_cat("surgery_performed_type"),
  pub_radiation       = spec_cat("radiation_therapy"),
  pub_surgical_approach = spec_cat("surgical_approach"),
  pub_aor_removed     = spec_num("total_aor_lnr", 0, 150),
  pub_pelv_removed    = spec_num("total_pelv_lnr", 0, 150))

#our variable for each original column (regex on the column name, first match wins); no match = no counterpart
#Herrmann's names come from his OpenML data, where one-hot dummies keep their variable's name as a prefix
OURS_RULES <- c(
  "^(sb|jd|hm)_age$" = "age", "^(sb|jd)_gender$|^hm_gender_" = "sex", "^sb_race$" = "race",
  "^sb_stage_|^hm_stage_event_" = "stage", "^jd_T$" = "T", "^jd_N$" = "N", "^jd_M$" = "M",
  "^(sb|jd)_histology$|^hm_(primary_pathology_)?histological_type_" = "histology",
  "^jd_grade$|^hm_neoplasm_histologic_grade_" = "grade",
  "^hm_diagnosis_subtype_" = "papillary", "^hm_number_of_lymphnodes_positive_by_he$" = "nodes_pos",
  "^hm_breast_carcinoma_estrogen_receptor_status_" = "ER", "^hm_breast_carcinoma_progesterone_receptor_status_" = "PR",
  "^hm_tobacco_smoking_history$" = "smoking", "^hm_vascular_tumor_cell_type_" = "vascular_inv",
  "^hm_maximum_tumor_dimension$" = "max_size", "^hm_primary_pathology_tumor_tissue_sites_" = "site",
  "^hm_primary_pathology_tumor_total_necrosis_percent_" = "necrosis", "^hm_breslow_depth_value$" = "breslow",
  "^hm_melanoma_ulceration_indicator_" = "ulceration", "^hm_pct_tumor_invasion$" = "myo_inv",
  #extended-set variables already defined in VAR_SPEC_EXT
  "^hm_venous_invasion_" = "venous_inv", "^hm_lymphatic_invasion_" = "lymph_inv",
  "^hm_primary_pathology_residual_tumor_" = "resid_tumor", "^hm_tumor_residual_disease_" = "residual",
  "^hm_alcohol_history_documented_" = "alcohol", "^hm_hemoglobin_result_" = "hb_q",
  "^hm_white_cell_count_result_" = "wbc_q", "^hm_visual_changes_" = "visual", "^hm_sensory_changes_" = "sensory",
  "^hm_motor_movement_changes_" = "motor", "^hm_tumor_location_" = "tumor_location",
  "^hm_albumin_result_specified_value$" = "albumin", "^hm_creatinine_value_in_mg_dl$" = "creatinine",
  "^hm_fetoprotein_outcome_value$" = "afp", "^hm_fibrosis_ishak_score_" = "fibrosis",
  "^hm_anatomic_neoplasm_subdivision_" = "location", "^hm_weight$" = "log_weight",
  #fields added above in PUB_SPEC
  "^hm_breast_carcinoma_surgical_procedure_name_" = "pub_surgery_breast",
  "^hm_lymphnode_neck_dissection_" = "pub_neck_dissection", "^hm_laterality_" = "pub_laterality",
  "^hm_surgery_performed_type_" = "pub_surgery_type", "^hm_radiation_therapy_" = "pub_radiation",
  "^hm_surgical_approach_" = "pub_surgical_approach",
  "^hm_total_aor_lnr$" = "pub_aor_removed", "^hm_total_pelv_lnr$" = "pub_pelv_removed")
#no counterpart in the CDR or PanCanAtlas files (so no ours version for those sets): LAML white cell count,
#blast % and FAB; SARC metastatic diagnosis; SKCM specimen site (Primary / Regional lymph node)

#our variable for each original column (NA = no counterpart)
.to_ours <- function(cols) vapply(cols, function(cl) {
  hit <- which(vapply(names(OURS_RULES), function(rx) grepl(rx, cl), logical(1)))
  if (length(hit)) unname(OURS_RULES[hit[1]]) else NA_character_
}, character(1))

#raw tables read once and reused across cancers
.pub_cache <- new.env()
.pub_cdr <- function() {
  if (is.null(.pub_cache$cdr))
    .pub_cache$cdr <- as.data.frame(read_xlsx(CDR_PATH, sheet = "TCGA-CDR", guess_max = 11000,
                                              .name_repair = "unique_quiet"))
  .pub_cache$cdr
}
.pub_pan <- function() {
  if (is.null(.pub_cache$pan))
    .pub_cache$pan <- fread(PANCAN_PATH, sep = "\t", header = TRUE, na.strings = SENT,
                            select = c("bcr_patient_barcode", "acronym", "pathologic_T", "pathologic_N", "pathologic_M"),
                            data.table = FALSE, showProgress = FALSE)
  .pub_cache$pan
}

#Jardillier's T/N/M coding: first two characters of the pathologic field, missing stays missing
.jd_tnm <- function(x) { x <- toupper(trimws(as.character(mark_na(x)))); out <- substr(x, 1, 2); out[is.na(x)] <- NA; out }

#ours: fetch the variables with load_cancer() (VAR_SPEC + level fixes), then apply the rare-level rule
#(under MIN_LEVEL_N patients -> missing) to categorical variables that have no hand-written fix here
.ours_columns <- function(code, vars, ids) {
  old <- VAR_SPEC
  on.exit(assign("VAR_SPEC", old, envir = .GlobalEnv), add = TRUE)
  assign("VAR_SPEC", c(old, PUB_SPEC[setdiff(names(PUB_SPEC), names(old))]), envir = .GlobalEnv)
  cl <- load_cancer(code, verbose = FALSE, vars = vars, read_expr = FALSE)$clinical
  m <- match(ids, cl$patient_id)
  if (anyNA(m)) stop(code, ": ", sum(is.na(m)), " modelling patients missing from the ours rebuild")
  fixed <- names(LEVEL_FIXES[[code]])
  out <- list()
  for (v in vars) {
    x <- cl[[v]][m]
    if (!is.numeric(x) && !v %in% fixed) {
      tb <- table(x); rare <- names(tb)[tb < MIN_LEVEL_N]
      if (length(rare)) { x[x %in% rare] <- NA
      cat(sprintf("  [rare -> NA] %s$%s: %s\n", code, v, paste(sprintf("%s(%d)", rare, as.integer(tb[rare])), collapse = ", "))) }
    }
    out[[paste0("o_", v)]] <- x
  }
  as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE)
}

#one cancer: every original recipe, its ours version where possible, and the specs rows
build_published <- function(code, verbose = TRUE) {
  dir.create(DIR_PUBLISHED, showWarnings = FALSE, recursive = TRUE)
  ids <- readRDS(path_modelling(code))$clinical$patient_id
  out <- data.frame(patient_id = ids, stringsAsFactors = FALSE)
  specs <- list()
  cdr <- .pub_cdr(); cdr <- cdr[cdr$type == code, ]
  m <- match(ids, cdr$bcr_patient_barcode)
  if (anyNA(m)) stop(code, ": ", sum(is.na(m)), " modelling patients not in the CDR")
  
  #SurvBoard, original (all CDR fields, raw levels, missing -> category in build_clinical)
  out$sb_age        <- suppressWarnings(as.numeric(cdr$age_at_initial_pathologic_diagnosis[m]))
  out$sb_gender     <- mark_na(cdr$gender[m])
  out$sb_race       <- mark_na(cdr$race[m])
  out$sb_stage_path <- mark_na(cdr$ajcc_pathologic_tumor_stage[m])
  out$sb_stage_clin <- mark_na(cdr$clinical_stage[m])
  out$sb_histology  <- mark_na(cdr$histological_type[m])
  sb <- grep("^sb_", names(out), value = TRUE)
  specs$survboard <- sb[vapply(sb, function(v) sum(!is.na(out[[v]])) > 0, logical(1))]  #all-missing columns carry nothing
  
  #Jardillier, original: completeness judged on every patient of the cancer, as in their script
  pan <- .pub_pan(); pan <- pan[pan$acronym == code, ]
  mp <- match(ids, pan$bcr_patient_barcode)
  jd <- c("jd_age", "jd_gender")
  out$jd_age <- out$sb_age; out$jd_gender <- out$sb_gender
  if (length(unique(na.omit(mark_na(cdr$gender)))) < 2) jd <- setdiff(jd, "jd_gender")  #one gender: not used
  cand <- list(jd_grade = list(all = mark_na(cdr$histological_grade), mine = mark_na(cdr$histological_grade[m])),
               jd_T = list(all = .jd_tnm(pan$pathologic_T), mine = .jd_tnm(pan$pathologic_T[mp])),
               jd_N = list(all = .jd_tnm(pan$pathologic_N), mine = .jd_tnm(pan$pathologic_N[mp])),
               jd_M = list(all = .jd_tnm(pan$pathologic_M), mine = .jd_tnm(pan$pathologic_M[mp])))
  jd_note <- character(0)
  for (v in names(cand)) {
    comp <- if (length(cand[[v]]$all)) mean(!is.na(cand[[v]]$all)) else 0  #field absent for this cancer
    if (comp >= JD_MIN_COMPLETE) { out[[v]] <- cand[[v]]$mine; jd <- c(jd, v) }
    else jd_note <- c(jd_note, sprintf("%s %.0f%%", sub("jd_", "", v), 100 * comp))
  }
  specs$jardillier <- jd
  
  #Herrmann, original: their numeric clinical columns joined by barcode
  hf <- file.path(DIR_PUBLISHED, sprintf("herrmann_clinical_%s.csv", code))
  hm_note <- "no Herrmann data for this cancer"
  if (file.exists(hf)) {
    h <- read.csv(hf, stringsAsFactors = FALSE, check.names = FALSE)
    mh <- match(ids, h$patient_id); absent <- is.na(mh)
    hv <- setdiff(names(h), "patient_id")
    #patients Herrmann did not have: one shared absent flag and a fixed 0 in his columns, so the columns
    #do not each grow an identical missing indicator (collinear, breaks the unpenalised fit); 0 is a constant, no leakage
    for (v in hv) { x <- suppressWarnings(as.numeric(h[[v]][mh])); x[absent] <- 0; out[[v]] <- x }
    if (any(absent)) { out$hm_absent <- as.numeric(absent); hv <- c(hv, "hm_absent") }
    if (anyNA(out[hv])) stop(code, ": missing values inside Herrmann's own columns")
    specs$herrmann <- hv
    hm_note <- sprintf("%d of %d modelling patients are in Herrmann's data", sum(!absent), length(ids))
  }
  
  #ours: one variable per idea, all of them, or none when any column has no counterpart here
  lat <- names(CANCER_VARS[[code]])
  ours <- list(); why_not <- list()
  for (p in names(specs)) {
    cols <- setdiff(specs[[p]], "hm_absent")
    to <- .to_ours(cols)
    if (anyNA(to)) { why_not[[p]] <- paste("no counterpart for", paste(cols[is.na(to)], collapse = ", ")); next }
    to <- unique(to)
    ours[[p]] <- c(lat[lat %in% to], setdiff(to, lat))  #lattice variables first, in lattice label order
  }
  need <- unique(unlist(ours))
  if (length(need)) out <- cbind(out, .ours_columns(code, need, ids))
  
  write.csv(out, file.path(DIR_PUBLISHED, sprintf("%s_published_clinical.csv", code)), row.names = FALSE)
  row <- function(p, version, vars, in_lat, note = "")
    data.frame(cancer = code, paper = p, version = version, variables = paste(vars, collapse = "+"),
               n_vars = length(vars), in_lattice = in_lat, note = note)
  rows <- do.call(rbind, c(
    lapply(names(specs), function(p) row(p, "original", specs[[p]], FALSE)),
    lapply(names(ours), function(p) {
      in_lat <- all(ours[[p]] %in% lat)
      row(p, "ours", if (in_lat) ours[[p]] else paste0("o_", ours[[p]]), in_lat)  #lattice label, or the o_ columns to fit
    }),
    lapply(names(why_not), function(p) row(p, "ours", character(0), FALSE, why_not[[p]]))))
  if (verbose) {
    cat(sprintf("\n=== %s === %d patients\n", code, length(ids)))
    for (p in names(specs)) cat(sprintf("  %-10s original : %s\n", p, paste(specs[[p]], collapse = ", ")))
    if (length(jd_note)) cat(sprintf("  jardillier left out (< 90%% complete): %s\n", paste(jd_note, collapse = ", ")))
    cat(sprintf("  herrmann   %s\n", hm_note))
    n_noage <- sum(is.na(out$sb_age))
    if (n_noage) cat(sprintf("  note: %d patient(s) without age (SurvBoard drops them; here age is median-filled)\n", n_noage))
    for (p in names(ours)) cat(sprintf("  %-10s ours     : %s  (%s)\n", p, paste(ours[[p]], collapse = "+"),
                                       if (all(ours[[p]] %in% lat)) "lattice" else "fitted"))
    for (p in names(why_not)) cat(sprintf("  %-10s ours     : not built, %s\n", p, why_not[[p]]))
  }
  rows
}

#every cancer, one specs file for all
build_published_all <- function(codes = CANCERS) {
  rows <- list(); failed <- character(0)
  for (code in codes) {
    r <- tryCatch(build_published(code), error = function(e) {
      cat(sprintf("[FAILED]     %-5s published: %s\n", code, conditionMessage(e))); NULL })
    if (is.null(r)) failed <- c(failed, code) else rows[[code]] <- r
  }
  sp <- do.call(rbind, rows)
  write.csv(sp, file.path(DIR_PUBLISHED, "published_specs.csv"), row.names = FALSE)
  cat(sprintf("\n[published] specs for %d cancers written to %s%s\n", length(rows),
              file.path(DIR_PUBLISHED, "published_specs.csv"),
              if (length(failed)) paste0(" | failed: ", paste(failed, collapse = ", ")) else ""))
  invisible(sp)
}