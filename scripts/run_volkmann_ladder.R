#Volkmann ladder (Volkmann et al. 2019, BMC Med Res Methodol 19:162, Table 1), refitted on BRCA
#Original: M1 nodal status N0/N+, M2 N0-N3 + T, M3 + ER, M4 + age as a (3,3) fractional polynomial;
#single mice imputation with PR as helper, fitted inside every training split (outer and inner).
#Ours: the same variables with this project's preprocessing, read off the main lattice (no refit).
#Usage, from the project root: source("master.R"); source("scripts/run_volkmann_ladder.R"); run_volkmann()

VK_CODE <- "BRCA"
VK_KEY  <- "BRCA_volkmann"
VK_IMPUTE_VARS <- c("vk_N4", "vk_T", "vk_ER", "vk_age", "vk_PR")  #his imputation model: ladder variables + PR, no outcome
VK_SPECS <- list(omics_only = character(0),
                 volkmann_M1 = "vk_N2",
                 volkmann_M2 = c("vk_N4", "vk_T"),
                 volkmann_M3 = c("vk_N4", "vk_T", "vk_ER"),
                 volkmann_M4 = c("vk_N4", "vk_T", "vk_ER", "vk_age"))
#the same sets in lattice labels (CANCER_VARS order), for the ours version
VK_OURS <- list(volkmann_M1 = "N", volkmann_M2 = c("T", "N"),
                volkmann_M3 = c("T", "N", "ER"), volkmann_M4 = c("age", "T", "N", "ER"))

#his codings, from the raw fields: pathologic T and N cut to two characters, X and indeterminate -> missing
.vk_tn <- function(x) { x <- toupper(substr(lc(x), 1, 2)); x[!grepl("^[TN][0-4]$", x)] <- NA; x }
.vk_receptor <- function(x) { low <- lc(x); out <- rep(NA_character_, length(low))
  out[low == "positive"] <- "positive"; out[low == "negative"] <- "negative"; out }

#his columns for the modelling cohort, matched by patient
vk_columns <- function(ids) {
  cdr <- as.data.frame(readxl::read_xlsx(CDR_PATH, sheet = "TCGA-CDR", guess_max = 11000,
                                         .name_repair = "unique_quiet"))
  cdr <- cdr[cdr$type == VK_CODE, ]
  pan <- data.table::fread(PANCAN_PATH, sep = "\t", header = TRUE, na.strings = SENT, data.table = FALSE,
                           showProgress = FALSE,
                           select = c("bcr_patient_barcode", "pathologic_T", "pathologic_N",
                                      "breast_carcinoma_estrogen_receptor_status",
                                      "breast_carcinoma_progesterone_receptor_status"))
  mc <- match(ids, cdr$bcr_patient_barcode); mp <- match(ids, pan$bcr_patient_barcode)
  if (anyNA(mc) || anyNA(mp)) stop("modelling patients missing from the CDR or PanCanAtlas file")
  T4 <- .vk_tn(pan$pathologic_T[mp]); T4[T4 %in% c("T0", "T1")] <- "T0_1"  #his T0_1 group
  N4 <- .vk_tn(pan$pathologic_N[mp])
  data.frame(patient_id = ids,
             vk_age = suppressWarnings(as.numeric(cdr$age_at_initial_pathologic_diagnosis[mc])),
             vk_T = T4, vk_N4 = N4,
             vk_N2 = ifelse(is.na(N4), NA, ifelse(N4 == "N0", "N0", "N+")),
             vk_ER = .vk_receptor(pan$breast_carcinoma_estrogen_receptor_status[mp]),
             vk_PR = .vk_receptor(pan$breast_carcinoma_progesterone_receptor_status[mp]),
             stringsAsFactors = FALSE)
}

#one split: mice fitted on the training rows only (ignore = test rows), then his age terms
#age/10 keeps the cubes small; the scale does not change the fit, since x^3 and x^3*log(x) span the same space for any scale
.vk_impute <- function(tr, te) {
  d <- rbind(tr[, VK_IMPUTE_VARS], te[, VK_IMPUTE_VARS])
  for (v in setdiff(VK_IMPUTE_VARS, "vk_age")) d[[v]] <- factor(d[[v]])
  ign <- c(rep(FALSE, nrow(tr)), rep(TRUE, nrow(te)))
  seed <- SEED + sum(strtoi(substr(tr$patient_id, 9, 12), 36L), na.rm = TRUE) %% 100000L  #fixed per split
  old <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = .GlobalEnv), add = TRUE)
  d <- if (anyNA(d[!ign, ]) || anyNA(d[ign, ]))
    mice::complete(mice::mice(d, m = 1, ignore = ign, seed = seed, printFlag = FALSE), 1) else d
  for (v in setdiff(VK_IMPUTE_VARS, "vk_age")) d[[v]] <- as.character(d[[v]])
  d$vk_N2 <- ifelse(d$vk_N4 == "N0", "N0", "N+")
  a <- d$vk_age / 10
  d$vk_age_fp1 <- a^3; d$vk_age_fp2 <- a^3 * log(a)
  list(train = d[!ign, , drop = FALSE], test = d[ign, , drop = FALSE])
}

#build_clinical for his columns: impute this split, swap age for its two terms, then the usual encoding
.build_clinical_base <- NULL
build_clinical_volkmann <- function(clinical_subset, tr, te) {
  if (!any(startsWith(clinical_subset, "vk_"))) return(.build_clinical_base(clinical_subset, tr, te))
  imp <- .vk_impute(tr, te)
  keep <- c(setdiff(clinical_subset, "vk_age"), if ("vk_age" %in% clinical_subset) c("vk_age_fp1", "vk_age_fp2"))
  .build_clinical_base(keep, imp$train, imp$test)
}

#the ours version: the four sets read off the main BRCA lattice summary, fused and offset arms
vk_lookup_ours <- function() {
  lat <- names(CANCER_VARS[[VK_CODE]])
  lab <- vapply(VK_OURS, function(v) paste(lat[lat %in% v], collapse = "+"), character(1))
  s <- read.csv(path_summary(VK_CODE), stringsAsFactors = FALSE)
  s <- s[s$arm %in% c("fused", "offset"), ]
  hit <- s$variables %in% lab
  if (!all(lab %in% s$variables)) warning("not in the lattice summary: ", paste(lab[!lab %in% s$variables], collapse = "; "))
  res <- cbind(spec = names(lab)[match(s$variables[hit], lab)], s[hit, ])
  write.csv(res, file.path(DIR_RESULTS, paste0(VK_KEY, "_ours.csv")), row.names = FALSE)
  invisible(res)
}

#fit the original ladder on the lattice's folds (same patients, outer and inner folds, repeats, arms)
run_volkmann <- function(force = FALSE) {
  if (!requireNamespace("mice", quietly = TRUE) || packageVersion("mice") < "3.12.0")
    stop("mice >= 3.12.0 is needed (for ignore =): renv::install('mice'), then renv::snapshot()")
  out <- path_results(VK_KEY)
  if (file.exists(out) && !force) {
    cat(sprintf("[volkmann] skipped (%s exists; force = TRUE to refit)\n", basename(out)))
  } else {
    if (force) {
      ck <- list.files(DIR_CHECKPOINTS, pattern = paste0("^", VK_KEY, "_repeat_[0-9]+\\.rds$"), full.names = TRUE)
      unlink(c(ck, file.path(DIR_CHECKPOINTS, paste0(VK_KEY, "_settings.txt")), out,
               sub("_results\\.csv$", "_genes.csv.gz", out), sub("_results\\.csv$", "_preds.rds", out),
               path_manifest(VK_KEY)))
    }
    dat <- readRDS(path_modelling(VK_CODE))
    vk <- vk_columns(dat$clinical$patient_id)
    for (v in setdiff(names(vk), "patient_id")) dat$clinical[[v]] <- vk[[v]]

    #what goes in: levels and missing counts of his columns
    for (v in c("vk_N2", "vk_N4", "vk_T", "vk_ER", "vk_PR"))
      cat(sprintf("  %-6s %s\n", v, paste(sprintf("%s=%d", names(table(vk[[v]], useNA = "ifany")),
                                                   table(vk[[v]], useNA = "ifany")), collapse = "  ")))
    cat(sprintf("  vk_age missing=%d\n", sum(is.na(vk$vk_age))))

    #his encoding replaces build_clinical for this run only; sequential, because workers re-source the library
    .build_clinical_base <<- build_clinical
    old_workers <- N_WORKERS
    build_clinical <<- build_clinical_volkmann; N_WORKERS <<- 1L
    on.exit({ build_clinical <<- .build_clinical_base; N_WORKERS <<- old_workers }, add = TRUE)

    res <- run_cancer(dat, VK_CODE, subsets = VK_SPECS, key = VK_KEY)

    #the offset arm has no penalized fallback in models.R: report its failures so we can decide
    n_off <- sum(!is.na(res$offset_error)); n_fb <- sum(!is.na(res$fused_fallback))
    cat(sprintf("[volkmann] fused cells rescued by penalized: %d | offset cells failed: %d (of %d)\n",
                n_fb, n_off, nrow(res)))
    if (n_off) warning("offset cells failed: these have no fallback yet, check before using the offset results")

    if (!exists("summarise_cancer", mode = "function")) source(file.path(DIR_SCRIPTS, "02_summarise.R"))
    summarise_cancer(VK_KEY, rescue = FALSE)
  }
  vk_lookup_ours()
  invisible(TRUE)
}
