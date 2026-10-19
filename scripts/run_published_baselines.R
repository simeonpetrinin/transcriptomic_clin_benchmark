#Published baselines: fits
#Fits each published clinical set (from build_published_baselines.R) on exactly the lattice's folds:
#same patients, same outer and inner folds, same repeats, same arms (clinical-only, IPF-LASSO fused, offset).
#Output key <CODE>_published (results, summary, specs): spec "<paper>" = the paper's original preprocessing,
#"<paper>_ours" = the same variables with this project's preprocessing. An ours set made only of main-lattice
#variables is not refitted: it is copied from the main lattice summary into results/<CODE>_published_ours.csv.
#Usage, from the project root: source("master.R"); source("scripts/build_published_baselines.R")
#                              source("scripts/run_published_baselines.R"); full_run_published()

published_key <- function(code) paste0(code, "_published")

#the specifications to fit for one cancer: omics-only reference, every original set, and the ours sets
#that are not lattice subsets
published_specs <- function(code) {
  sp <- read.csv(file.path(DIR_PUBLISHED, "published_specs.csv"), stringsAsFactors = FALSE)
  sp <- sp[sp$cancer == code & sp$n_vars > 0 & !sp$in_lattice, ]
  if (!any(sp$version == "original")) stop("no published specifications for ", code, ", run build_published_all() first")
  nm <- ifelse(sp$version == "original", sp$paper, paste0(sp$paper, "_ours"))
  c(list(omics_only = character(0)), setNames(strsplit(sp$variables, "+", fixed = TRUE), nm))
}

#the modelling object with the published columns added, matched by patient and kept in row order
.with_published <- function(code) {
  dat <- readRDS(path_modelling(code))
  pub <- read.csv(file.path(DIR_PUBLISHED, sprintf("%s_published_clinical.csv", code)), stringsAsFactors = FALSE)
  m <- match(dat$clinical$patient_id, pub$patient_id)
  if (anyNA(m)) stop(code, ": the published table is missing ", sum(is.na(m)), " modelling patients, rebuild it")
  for (v in setdiff(names(pub), "patient_id")) {
    if (v %in% names(dat$clinical)) stop(code, ": column ", v, " already in the modelling data")
    dat$clinical[[v]] <- pub[[v]][m]
  }
  dat
}

#one cancer: fit the original and non-lattice ours sets, summarise, and look up the lattice ours sets
run_published <- function(code, force = FALSE, n_genes = TOP_N_GENES) {
  key <- published_key(code); out <- path_results(key)
  if (file.exists(out) && !force) {
    cat(sprintf("[published] %-5s skipped (%s exists)\n", code, basename(out)))
  } else {
    if (force) {  #clear this key's checkpoints and outputs, as run_robustness() does
      ck <- list.files(DIR_CHECKPOINTS, pattern = paste0("^", key, "_repeat_[0-9]+\\.rds$"), full.names = TRUE)
      unlink(c(ck, file.path(DIR_CHECKPOINTS, paste0(key, "_settings.txt")), out,
               sub("_results\\.csv$", "_genes.csv.gz", out), sub("_results\\.csv$", "_preds.rds", out), path_manifest(key)))
    }
    sp <- published_specs(code)
    cat(sprintf("[published] %-5s %s | %d repeats\n", code,
                paste(sprintf("%s (%d)", names(sp), lengths(sp)), collapse = ", "), n_repeats_for(code)))
    run_cancer(.with_published(code), code, n_genes = n_genes, subsets = sp, key = key)
    if (!exists("summarise_cancer", mode = "function")) source(file.path(DIR_SCRIPTS, "02_summarise.R"))
    summarise_cancer(key, rescue = FALSE)  #these fits were never rescued
  }
  lookup_ours(code)
  invisible(TRUE)
}

#ours sets made only of lattice variables: copy their rows from the main summary
lookup_ours <- function(code) {
  sp <- read.csv(file.path(DIR_PUBLISHED, "published_specs.csv"), stringsAsFactors = FALSE)
  sp <- sp[sp$cancer == code & sp$version == "ours" & sp$n_vars > 0 & sp$in_lattice, ]
  if (!nrow(sp)) return(invisible(NULL))
  s <- read.csv(path_summary(code), stringsAsFactors = FALSE)
  s <- s[s$arm == "fused", ]
  hit <- match(sp$variables, s$variables)
  if (anyNA(hit)) warning(code, ": not in the lattice: ", paste(sp$variables[is.na(hit)], collapse = "; "), call. = FALSE)
  res <- cbind(sp[, c("cancer", "paper", "version", "variables", "n_vars")],
               s[hit, c("subset_id", "c_ref", "c_model", "delta_c", "delta_c_sd", "ibs_ref", "ibs_model", "delta_ibs", "delta_ibs_sd")])
  write.csv(res, file.path(DIR_RESULTS, sprintf("%s_published_ours.csv", code)), row.names = FALSE)
  invisible(res)
}

#every cancer in turn; one failing does not stop the rest
full_run_published <- function(codes = CANCERS, force = FALSE) {
  t0 <- proc.time()[["elapsed"]]; failed <- character(0)
  for (code in codes) {
    ok <- tryCatch({
      if (exists("timed_stage", mode = "function")) timed_stage(code, "published", run_published(code, force))
      else run_published(code, force)
      TRUE
    }, error = function(e) { cat(sprintf("[FAILED]     %-5s published: %s\n", code, conditionMessage(e))); FALSE })
    if (!ok) failed <- c(failed, code)
    invisible(gc())
  }
  cat(sprintf("\n[published] done in %.1f min | %d ok, %d failed%s\n", (proc.time()[["elapsed"]] - t0) / 60,
              length(codes) - length(failed), length(failed),
              if (length(failed)) paste0(": ", paste(failed, collapse = ", ")) else ""))
  invisible(list(failed = failed))
}