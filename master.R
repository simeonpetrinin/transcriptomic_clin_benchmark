#Master
#The single entry point. Stages: audit (once), prepare, run, summarise for the main
#lattice; full_run_extended() and full_run_robustness() for the secondary runs.
#Source it in a clean session from the project root; it defines everything, runs nothing.

#Setup: check the package versions against renv.lock before anything is loaded.
check_versions <- function(strict = TRUE){  #strict = TRUE stops on drift
  if (!requireNamespace("renv", quietly = TRUE)){
    warning("renv not installed -- package versions are NOT pinned");
    return(invisible(FALSE))
  }
  st <- tryCatch(renv::status(), error = function(e) NULL) 
  ok <- isTRUE(st$synchronized)  #st$synchronized is not TRUE when st is NULL
  if (!ok) {
    msg <- paste("Package library is out of sync with renv.lock.",
                 "Run renv::restore() to match the lockfile, or renv::snapshot()",
                 "to adopt what is installed -- then re-source master.R.",
                 "Do NOT change versions part-way through a run: the manifests",
                 "record one version per cohort and a mid-run switch makes them",
                 "disagree for a reason nobody would think to look for.")
    if (strict) stop(msg) else warning(msg)
  }
  invisible(ok)
}
check_versions()

source("paths.R")
source_library()
source(file.path(DIR_SCRIPTS, "run_robustness.R"))

#the cancers in scope, written by scripts/select_cancers.R (falls back to LUAD)
CANCERS <- tryCatch(
  readLines(file.path(DIR_RESULTS, "final_cancer_list.txt")),
  error = function(e) "LUAD"
)

ALL_STAGES <- c("prepare", "run", "summarise")

#Stage 1, audit (once): which cancers enter and which clinical variables pass.
run_audit <- function(force = FALSE){
  out <- file.path(DIR_RESULTS, "final_cancer_list.txt")
  if (file.exists(out) && !force) {
    CANCERS <<- readLines(out)
    cat("[audit]      skipped (final_cancer_list.txt exists; force = TRUE to redo)\n")
    return(invisible(FALSE))
  }
  cat("[audit]      running feasibility + clinical completeness audits ...\n")
  source(file.path(DIR_SCRIPTS, "tcga_feasibility_audit.R"))
  source(file.path(DIR_SCRIPTS, "select_cancers.R"))
  source(file.path(DIR_SCRIPTS, "audit_clinical.R"))
  CANCERS <<- readLines(out)
  invisible(TRUE)
}

#Stage 2, prepare: build one cancer's modelling object.
prepare_cancer <- function(code, force = FALSE){
  out <- path_modelling(code)
  if (file.exists(out) && !force) {
    have <- names(readRDS(out)$clinical)
    need <- names(CANCER_VARS[[code]])
    if (all(need %in% have)) {
      cat(sprintf("[prepare]    %-5s skipped (%s exists)\n", code, basename(out)))
      return(invisible(FALSE))
    }
    cat(sprintf("[prepare]    %-5s %s lacks %s, rebuilding\n", code, basename(out),
                paste(setdiff(need, have), collapse = ", ")))
  }
  cat(sprintf("[prepare]    %-5s building modelling object ...\n", code))
  dat <- build_modelling_data(code)  #defined in R/load_data.R
  cat(sprintf("[prepare]    %-5s saved %s (%d patients, %d genes)\n",
              code, basename(out), nrow(dat$clinical), ncol(dat$expr)))
  invisible(TRUE)
}

#Stage 3, run: repeated nested CV across the subset lattice
run_cv <- function(code, force = FALSE){
  out <- path_results(code)
  if (file.exists(out) && !force) {
    cat(sprintf("[run]        %-5s skipped (%s exists)\n", code, basename(out)))
    return(invisible(FALSE))
  }
  if (!file.exists(path_modelling(code)))
    stop("no modelling object for ", code, " -- run the 'prepare' stage first")
  
  if (force) {
    ckpts <- list.files(DIR_CHECKPOINTS,
                        pattern = paste0("^", code, "_repeat_[0-9]+\\.rds$"),
                        full.names = TRUE)
    settings <- file.path(DIR_CHECKPOINTS, paste0(code, "_settings.txt"))
    old_run <- c(out,
                 sub("_results\\.csv$", "_genes.csv.gz", out),
                 sub("_results\\.csv$", "_preds.rds", out),
                 path_manifest(code))
    unlink(c(ckpts, settings, old_run))
    cat(sprintf("[run]        %-5s force = TRUE: cleared old CV checkpoints/results\n", code))
  }
  
  #write the subset lattice for the record (inputs live in data/processed)
  idx <- build_subset_index(code, include_anchor = TRUE)
  write.csv(idx, path_subsets(code), row.names = FALSE)
  
  cat(sprintf("[run]        %-5s %d subsets x %d folds x %d repeats ...\n",
              code, nrow(idx), OUTER_FOLDS, n_repeats_for(code)))
  dat <- readRDS(path_modelling(code))
  run_cancer(dat, code)  #defined in R/cv.R
  invisible(TRUE)
}

#Stage 4, summarise: aggregate + check figure
summarise <- function(code, force = FALSE){
  out <- path_summary(code)
  if (file.exists(out) && !force){
    cat(sprintf("[summarise]  %-5s skipped (%s exists)\n", code, basename(out)))
    return(invisible(FALSE))
  }
  if (!file.exists(path_results(code)))
    stop("no results for ", code, " -- run the 'run' stage first")
  cat(sprintf("[summarise]  %-5s aggregating + drawing ...\n", code))
  source(file.path(DIR_SCRIPTS, "02_summarise.R"))
  summarise_cancer(code)
  invisible(TRUE)
}

#Stage timing: results/timing/stage_log.csv, one row per cohort and stage (appended)
timed_stage <- function(code, stage, expr) {
  t0  <- proc.time()[["elapsed"]]; start <- Sys.time()
  ok  <- FALSE
  on.exit({
    dir <- file.path(DIR_RESULTS, "timing"); dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    p   <- file.path(dir, "stage_log.csv")
    row <- data.frame(start = format(start, "%Y-%m-%d %H:%M:%S"), cancer = code, stage = stage,
                      wall_secs = round(proc.time()[["elapsed"]] - t0, 1), completed = ok,
                      host = tryCatch(as.character(Sys.info()[["nodename"]]), error = function(e) NA),
                      workers = if (exists("N_WORKERS")) N_WORKERS else NA)
    write.table(row, p, sep = ",", row.names = FALSE, col.names = !file.exists(p),
                append = file.exists(p))
  }, add = TRUE)
  val <- force(expr)
  ok  <- TRUE
  invisible(val)
}

#The main pipeline
run_pipeline <- function(codes = CANCERS, stages = ALL_STAGES, force = FALSE) {
  stopifnot(all(stages %in% ALL_STAGES))
  t0 <- proc.time()[["elapsed"]]
  
  cat("\n=========================================================\n")
  cat(sprintf(" pipeline: %d cancer(s) | stages: %s\n",
              length(codes), paste(stages, collapse = " -> ")))
  cat("=========================================================\n")
  
  failed <- character(0)
  for (code in codes) {
    cat(sprintf("\n--- %s ---\n", code))
    ok <- tryCatch({
      if ("prepare"   %in% stages) timed_stage(code, "prepare",   prepare_cancer(code, force))
      if ("run"       %in% stages) timed_stage(code, "run",       run_cv(code, force))
      if ("summarise" %in% stages) timed_stage(code, "summarise", summarise(code, force))
      TRUE
    }, error = function(e) {
      #one cancer failing must not kill a 21-cancer sweep
      cat(sprintf("[FAILED]     %-5s %s\n", code, conditionMessage(e)))
      FALSE
    })
    if (!ok) failed <- c(failed, code)
    #free this cancer's memory before the next one (the workers already exited)
    invisible(gc())
  }
  
  mins <- (proc.time()[["elapsed"]] - t0) / 60
  cat("\n=========================================================\n")
  cat(sprintf(" done in %.1f min | %d ok, %d failed\n",
              mins, length(codes) - length(failed), length(failed)))
  if (length(failed)) cat(" failed: ", paste(failed, collapse = ", "), "\n", sep = "")
  cat("=========================================================\n")
  invisible(list(failed = failed, minutes = mins))
}

#Extended run: every eligible clinical variable in one specification (R/extended.R)
full_run_extended <- function(codes = CANCERS, force = FALSE) {
  failed <- character(0)
  for (code in intersect(codes, names(EXT_VARS))) {
    key <- ext_key(code)
    if (file.exists(path_results(key)) && !force) {
      cat(sprintf("[extended]   %-5s skipped (%s exists)\n", code, basename(path_results(key)))); next
    }
    ok <- tryCatch({
      if (force) {
        ckpts <- list.files(DIR_CHECKPOINTS, pattern = paste0("^", key, "_repeat_[0-9]+\\.rds$"),
                            full.names = TRUE)
        unlink(c(ckpts, file.path(DIR_CHECKPOINTS, paste0(key, "_settings.txt")), path_results(key)))
      }
      timed_stage(code, "extended", run_extended(code, force_build = force))
      if (!exists("summarise_cancer", mode = "function"))
        source(file.path(DIR_SCRIPTS, "02_summarise.R"))
      summarise_cancer(key)
      TRUE
    }, error = function(e) {
      cat(sprintf("[FAILED]     %-5s extended: %s\n", code, conditionMessage(e))); FALSE })
    if (!ok) failed <- c(failed, code)
    invisible(gc())
  }
  skipped <- setdiff(codes, names(EXT_VARS))
  if (length(skipped)) cat("[extended]   no extended variables: ", paste(skipped, collapse = ", "), "\n", sep = "")
  invisible(list(failed = failed))
}

cat("master.R loaded.\n")
cat("  run_audit()                          # stage 1, once\n")
cat("  run_pipeline(\"LUAD\")                 # one cancer, end to end\n")
cat("  run_pipeline(CANCERS)                # the full sweep\n")
cat("  full_run_extended()                  # every eligible variable together, 10 repeats\n")
cat("  full_run_robustness()                # sensitivity arms on omics, age, main set + 10 sampled\n")
cat(sprintf("  CANCERS = %d cancer(s): %s\n", length(CANCERS),
            paste(head(CANCERS, 5), collapse = ", ")))