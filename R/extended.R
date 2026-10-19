#Extended
#Clinical with every eligible variable: the main set plus every other eligible variable,
#fitted together as one specification (ext_full) and run with EXT_REPEATS repeats.
#Same learners, folds and missing-data handling as the main lattice.

BODY <- c("log_weight", "log_height")  #body size enters as one item

EXT_REPEATS  <- 10L  #repeats for the extended run, whatever the lattice size
EXT_ADD_DROP <- FALSE  #TRUE also runs add-one (main + one item) and leave-one-out specs

#EXT_VARS[[code]]: named items, each one or more VAR_SPEC variables that enter together.
EXT_VARS <- list(
  BLCA = list(lvi = "lvi", body_size = BODY, hist_non_mibc = "hist_non_mibc",
              pack_years = "pack_years"),
  BRCA = list(location = "location", margin = "margin", menopause = "menopause"),
  CESC = list(ecog = "ecog", body_size = BODY, ethnicity = "ethnicity",
              keratinizing = "keratinizing", menopause = "menopause",
              nodes_pos = "nodes_pos"),
  COAD = list(cea = "cea", mmr = "mmr", body_size = BODY, location = "location",
              resid_tumor = "resid_tumor", synchronous = "synchronous",
              lymph_inv = "lymph_inv", venous_inv = "venous_inv"),
  ESCA = list(barretts = "barretts", alcohol_amount = "alcohol_amount",
              alcohol_freq = "alcohol_freq", alcohol = "alcohol", body_size = BODY,
              pack_years = "pack_years", prior_esoph = "prior_esoph",
              reflux = "reflux", resid_tumor = "resid_tumor"),
  GBM  = list(kps = "kps"),
  HNSC = list(lvi = "lvi", alcohol = "alcohol", pack_years = "pack_years",
              perineural = "perineural"),
  KIRC = list(ethnicity = "ethnicity", hb_q = "hb_q", plt_q = "plt_q", wbc_q = "wbc_q"),
  KIRP = list(body_size = BODY, hb_q = "hb_q", plt_q = "plt_q", wbc_q = "wbc_q"),
  LGG  = list(kps = "kps", allergy_animal = "allergy_animal", asthma = "asthma",
              fam_hist = "fam_hist", allergy_food = "allergy_food",
              hay_fever = "hay_fever", headache = "headache",
              mental_status = "mental_status", motor = "motor",
              allergy_mould = "allergy_mould", preop_antiseizure = "preop_antiseizure",
              preop_steroids = "preop_steroids", presenting_symptom = "presenting_symptom",
              seizures = "seizures", sensory = "sensory",
              supratentorial = "supratentorial", symptom_duration = "symptom_duration",
              tumor_location = "tumor_location", visual = "visual"),
  LIHC = list(afp = "afp", aetiology = "aetiology", albumin = "albumin",
              body_size = BODY, creatinine = "creatinine", fibrosis = "fibrosis",
              liver_inflam = "liver_inflam", plt_count = "plt_count",
              prothrombin = "prothrombin", resid_tumor = "resid_tumor",
              child_pugh = "child_pugh"),
  LUAD = list(location = "location", pack_years = "pack_years"),
  LUSC = list(location = "location", pack_years = "pack_years",
              resid_tumor = "resid_tumor"),
  OV   = list(residual = "residual"),
  PAAD = list(alcohol_exposure = "alcohol_exposure", pancreatitis = "pancreatitis",
              diabetes = "diabetes", fam_hist = "fam_hist", location = "location",
              resid_tumor = "resid_tumor"),
  READ = list(cea = "cea", mmr = "mmr", location = "location",
              relatives_cancer = "relatives_cancer", resid_tumor = "resid_tumor",
              lymph_inv = "lymph_inv", venous_inv = "venous_inv"),
  SARC = list(margin = "margin", multifocal = "multifocal", resid_tumor = "resid_tumor",
              tumour_depth_cm = "tumour_depth_cm", tumour_width = "tumour_width"),
  SKCM = list(body_size = BODY, clark = "clark", unknown_primary = "unknown_primary"),
  STAD = list(barretts = "barretts", location = "location", reflux = "reflux",
              resid_tumor = "resid_tumor", fam_hist_stomach = "fam_hist_stomach"),
  UCEC = list(body_size = BODY, menopause = "menopause",
              peritoneal_wash = "peritoneal_wash", resid_tumor = "resid_tumor")
)
#LAML: no extended variables (only age, sex and race are recorded)

#stop if an extended item repeats a main-set variable
for (.c in names(EXT_VARS)) {
  .dup <- intersect(unlist(EXT_VARS[[.c]]), names(CANCER_VARS[[.c]]))
  if (length(.dup)) stop("EXT_VARS$", .c, " repeats main-set variable(s): ", paste(.dup, collapse = ", "))
}
rm(.c, .dup)

ext_key <- function(code) paste0(code, "_ext")

#all variables the extended modelling object must contain
extended_vars <- function(code) {
  main <- names(CANCER_VARS[[code]])
  unique(c(main, unlist(EXT_VARS[[code]], use.names = FALSE)))
}

#specifications for run_cancer(): ext_full, plus add/drop specs if EXT_ADD_DROP
extended_specs <- function(code) {
  items <- EXT_VARS[[code]]
  if (!length(items)) stop("no EXT_VARS for ", code)
  main  <- names(CANCER_VARS[[code]])
  ext   <- extended_vars(code)
  specs <- list(ext_full = ext)
  if (isTRUE(EXT_ADD_DROP) && length(items) > 1) {  #with one item, add == ext_full and drop == main_full
    specs <- c(list(main_full = main), specs)
    for (nm in names(items)) {
      specs[[paste0("add:",  nm)]] <- c(main, items[[nm]])
      specs[[paste0("drop:", nm)]] <- setdiff(ext, items[[nm]])
    }
  }
  specs
}

#build (or reuse) the extended modelling object, then run the specifications
run_extended <- function(code, force_build = FALSE, n_genes = TOP_N_GENES) {
  key <- ext_key(code)
  dat <- if (!force_build && file.exists(path_modelling(key))) readRDS(path_modelling(key)) else NULL
  #rebuild when absent, forced, or missing a variable EXT_VARS now asks for
  if (is.null(dat) || !all(extended_vars(code) %in% names(dat$clinical)))
    dat <- build_modelling_data(code, vars = extended_vars(code), key = key)
  run_cancer(dat, code, n_genes = n_genes, subsets = extended_specs(code), key = key,
             n_repeats = EXT_REPEATS)
}