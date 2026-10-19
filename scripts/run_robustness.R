#Robustness
#Refits a fixed set of specifications per cancer (omics only, age alone, the full main set
#and a size-stratified random sample of the lattice) with the sensitivity arms switched on
#beside the fused and offset arms. Output key <CODE>_robust; run with full_run_robustness().
#The specifications actually run are read from results/robustness_specs.csv when that file exists, so a rerun uses
#exactly the reported subsets; otherwise they are drawn with the nested sampler below.

ROBUST_ARMS <- c("forest", "coxboost", "ridge", "score")

#the arm switches in cv.R that each name turns on
ROBUST_FLAGS <- c(forest = "RUN_FOREST_ARM", coxboost = "RUN_COXBOOST_ARM",
                  ridge = "RUN_RIDGE_ARM", score = "RUN_SCORE_ARM")

ROBUST_N_SAMPLED <- 10L  #random lattice subsets per cancer, beside the three fixed specs
ROBUST_SPECS_FILE <- file.path(DIR_RESULTS, "robustness_specs.csv")  #the subsets used in the paper
ROBUST_SEED      <- 20260928L  #fixed seed for the sample, set before any robustness result

#slots per subset size, added one at a time to the size with the fewest so far (ties: nearest the middle size,
#then the smaller one), capped by what each size holds; n + 1 is always n plus one slot, so allocations are nested
.robust_alloc <- function(avail, n) {
  S <- length(avail); alloc <- integer(S); mid <- (S + 1) / 2
  while (sum(alloc) < min(n, sum(avail))) {
    room <- which(alloc < avail)
    cand <- room[alloc[room] == min(alloc[room])]
    j <- cand[order(abs(cand - mid), cand)][1]
    alloc[j] <- alloc[j] + 1L
  }
  alloc
}

#the lattice subsets the sample is drawn from: all but age alone and the full set
robust_pool <- function(code) {
  vars <- names(CANCER_VARS[[code]])
  Filter(function(x) length(x) < length(vars) && !identical(x, "age"), all_subsets(vars))
}

#the reported specs for one cancer, read back from ROBUST_SPECS_FILE (NULL if the file or cancer is absent)
.saved_robust_specs <- function(code, path = ROBUST_SPECS_FILE) {
  if (!file.exists(path)) return(NULL)
  s <- read.csv(path, stringsAsFactors = FALSE)
  s <- s[s$cancer == code, , drop = FALSE]
  if (!nrow(s)) return(NULL)
  out <- lapply(s$variables, function(v) if (v == ANCHOR_LABEL) character(0) else strsplit(v, "+", fixed = TRUE)[[1]])
  setNames(out, s$spec)
}

#fixed specs (omics only, age alone, full main set) plus the stratified random sample;
#the saved specs are used unless redraw = TRUE
robust_specs <- function(code, n = ROBUST_N_SAMPLED, seed = ROBUST_SEED + sum(utf8ToInt(code)),
                         pool = robust_pool(code), redraw = FALSE) {
  if (!redraw) {
    saved <- .saved_robust_specs(code)
    if (!is.null(saved)) return(saved)
  }
  vars  <- names(CANCER_VARS[[code]])
  fixed <- list(omics_only = character(0), age_only = "age", main_full = vars)
  if (length(pool) <= n) picked <- pool else {
    sz    <- sort(unique(lengths(pool)))
    avail <- vapply(sz, function(k) sum(lengths(pool) == k), integer(1))
    alloc <- .robust_alloc(avail, n)
    set.seed(seed)
    ord <- lapply(avail, sample.int)  #one fixed random order per size, drawn the same way whatever n is
    picked <- unlist(lapply(seq_along(sz), function(i) {
      cand <- pool[lengths(pool) == sz[i]]
      cand[sort(ord[[i]][seq_len(alloc[i])])]  #the first alloc[i] of that order
    }), recursive = FALSE)
  }
  names(picked) <- sprintf("sample_%03d", seq_along(picked))
  c(fixed, picked)
}

#write every cancer's robustness specs to one file, before anything is run (draws afresh)
write_robust_specs <- function(codes = CANCERS) {
  out <- do.call(rbind, lapply(codes, function(code) {
    sp <- robust_specs(code, redraw = TRUE)
    data.frame(cancer = code, spec = names(sp), n_vars = lengths(sp),
               variables = vapply(sp, function(x)
                 if (length(x)) paste(x, collapse = "+") else ANCHOR_LABEL, character(1)),
               row.names = NULL)
  }))
  f <- file.path(DIR_RESULTS, "robustness_specs.csv")
  write.csv(out, f, row.names = FALSE)
  cat(sprintf("[robustness] %d specs for %d cancers written to %s\n", nrow(out), length(codes), f))
  invisible(out)
}

robust_key <- function(code) paste0(code, "_robust")

#fallback set if all cancers are too slow: the n cancers with most deaths and the n
robustness_fallback_codes <- function(codes = CANCERS, n = 5) {
  deaths <- vapply(codes, function(code) {
    f <- path_modelling(code)
    if (!file.exists(f)) return(NA_real_)
    cl <- readRDS(f)$clinical
    sum(cl$OS == 1 & !is.na(cl$OS.time) & cl$OS.time > 0, na.rm = TRUE)
  }, numeric(1))
  if (anyNA(deaths)) stop("no modelling object for: ", paste(codes[is.na(deaths)], collapse = ", "))
  ord <- names(sort(deaths))
  unique(c(tail(ord, n), head(ord, n)))
}

#one cancer: switch the chosen arms on, run the two specifications, switch them back
run_robustness <- function(code, arms = ROBUST_ARMS, force = FALSE, n_genes = TOP_N_GENES) {
  arms <- match.arg(arms, ROBUST_ARMS, several.ok = TRUE)
  key  <- robust_key(code)
  out  <- path_results(key)
  if (file.exists(out) && !force) {
    cat(sprintf("[robustness] %-5s skipped (%s exists)\n", code, basename(out)))
    return(invisible(FALSE))
  }
  if (!file.exists(path_modelling(code)))
    stop("no modelling object for ", code, ", run the 'prepare' stage first")
  if (force) {
    ckpts <- list.files(DIR_CHECKPOINTS, pattern = paste0("^", key, "_repeat_[0-9]+\\.rds$"),
                        full.names = TRUE)
    unlink(c(ckpts, file.path(DIR_CHECKPOINTS, paste0(key, "_settings.txt")), out,
             sub("_results\\.csv$", "_genes.csv.gz", out),
             sub("_results\\.csv$", "_preds.rds", out), path_manifest(key)))
  }
  #switch the arms on in the global settings (exported to the workers), restore on exit
  flags <- ROBUST_FLAGS[arms]
  old <- mget(unname(ROBUST_FLAGS), envir = .GlobalEnv)
  on.exit(for (f in names(old)) assign(f, old[[f]], envir = .GlobalEnv), add = TRUE)
  for (f in unname(ROBUST_FLAGS)) assign(f, f %in% flags, envir = .GlobalEnv)
  cat(sprintf("[robustness] %-5s arms: %s | %d specs | %d repeats\n",
              code, paste(arms, collapse = ", "), length(robust_specs(code)), n_repeats_for(code)))
  dat <- readRDS(path_modelling(code))
  run_cancer(dat, code, n_genes = n_genes, subsets = robust_specs(code), key = key)
  if (!exists("summarise_cancer", mode = "function"))
    source(file.path(DIR_SCRIPTS, "02_summarise.R"))
  summarise_cancer(key)
  invisible(TRUE)
}

#every cancer in turn
full_run_robustness <- function(codes = CANCERS, arms = ROBUST_ARMS, force = FALSE) {
  t0 <- proc.time()[["elapsed"]]
  failed <- character(0)
  for (code in codes) {
    ok <- tryCatch({
      if (exists("timed_stage", mode = "function"))
        timed_stage(code, "robustness", run_robustness(code, arms, force))
      else run_robustness(code, arms, force)
      TRUE
    }, error = function(e) {
      cat(sprintf("[FAILED]     %-5s robustness: %s\n", code, conditionMessage(e))); FALSE })
    if (!ok) failed <- c(failed, code)
    invisible(gc())
  }
  cat(sprintf("\n[robustness] done in %.1f h | %d ok, %d failed%s\n",
              (proc.time()[["elapsed"]] - t0) / 3600, length(codes) - length(failed),
              length(failed), if (length(failed)) paste0(": ", paste(failed, collapse = ", ")) else ""))
  invisible(list(failed = failed))
}
#check that the sampler is nested: every sample of size n sits inside the sample of size n + 1 (TRUE expected for every cancer)
check_robust_nesting <- function(codes = CANCERS, ns = 10:200) {
  ok <- vapply(codes, function(code) {
    pool <- robust_pool(code)
    key <- function(n) {
      sp <- robust_specs(code, n, pool = pool, redraw = TRUE)
      vapply(sp[grepl("^sample_", names(sp))], paste, character(1), collapse = "+")
    }
    prev <- key(ns[1]); nested <- TRUE
    for (n in ns[-1]) { cur <- key(n); nested <- nested && all(prev %in% cur); prev <- cur }
    nested
  }, logical(1))
  print(ok)
  invisible(ok)
}