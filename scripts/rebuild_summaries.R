#Rebuild summaries with the rescue
#Writes the subset summaries and per-repeat files with the rescued fused cells kept, for the
#cancers that were rescued or never summarised, plus "_norescue" sensitivity versions.
#Run after source("master.R"): rebuild_summaries()

RESCUED    <- c("BRCA", "CESC", "COAD", "KIRC", "KIRP", "READ", "SARC", "UCEC")  #have a <CODE>_rescue.csv
NO_SUMMARY <- c("BLCA", "ESCA")  #never summarised, nothing to rescue

#one key at a time; one failure must not stop the rest
rebuild_summaries <- function(main = c(NO_SUMMARY, RESCUED), rescued = RESCUED, robust = TRUE) {
  source(file.path(DIR_SCRIPTS, "02_summarise.R"))
  jobs <- c(lapply(main, function(k) list(key = k, rescue = TRUE)),
            lapply(rescued, function(k) list(key = k, rescue = FALSE)))  #sensitivity: failed cells dropped
  if (robust) jobs <- c(jobs,
                        lapply(paste0(rescued, "_robust"), function(k) list(key = k, rescue = TRUE)),
                        lapply(paste0(rescued, "_robust"), function(k) list(key = k, rescue = FALSE)))
  failed <- character(0)
  for (j in jobs) {
    lab <- paste0(j$key, if (j$rescue) "" else " (no rescue)")
    cat(sprintf("\n######## %s ########\n", lab))
    ok <- tryCatch({ summarise_cancer(j$key, rescue = j$rescue); TRUE },
                   error = function(e) { cat("[FAILED]", lab, conditionMessage(e), "\n"); FALSE })
    if (!ok) failed <- c(failed, lab)
    invisible(gc())
  }
  cat(sprintf("\nDone: %d ok, %d failed%s\n", length(jobs) - length(failed), length(failed),
              if (length(failed)) paste0(" (", paste(failed, collapse = ", "), ")") else ""))
  invisible(failed)
}