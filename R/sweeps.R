#Sweeps
#run_sweeps() runs the gene, PFI and missing-data sweeps, together or apart, for
#chosen cancers. Each writes its own keyed files and checkpoints; defines only.
SWEEP_NAMES <- c("genes", "pfi", "missing")

run_sweeps <- function(codes = names(CANCER_VARS), which = SWEEP_NAMES) {
  which <- match.arg(which, SWEEP_NAMES, several.ok = TRUE)
  t0 <- proc.time()[["elapsed"]]
  status <- c()
  if ("genes"   %in% which) status <- c(status, sweep_genes(codes))
  if ("pfi"     %in% which) status <- c(status, sweep_pfi(codes))
  if ("missing" %in% which) status <- c(status, sweep_missing(codes))
  cat(sprintf("\n=== sweeps done in %.1f h | %d ok, %d failed ===\n",
              (proc.time()[["elapsed"]] - t0) / 3600,
              sum(status == "ok"), sum(status != "ok")))
  bad <- status[status != "ok"]
  if (length(bad)) for (k in names(bad)) cat(sprintf("  %s: %s\n", k, bad[[k]]))
  invisible(status)
}