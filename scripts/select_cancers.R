#Select cancers
#Keeps cancers with at least 100 patients, 10 deaths and a 5% event ratio, and writes
#results/final_cancer_list.txt and cancer_selection_reasons.csv.

#Reads : results/tcga_feasibility_audit.csv   (generated from tcga_feasibility_audit.R)

if (!exists("DIR_RESULTS")) source("paths.R")  #standalone-safe

#selected thresholds
MIN_SAMPLES      <- 100  #n_analysis: patients with BOTH clinical + tumour RNA
MIN_EVENTS       <- 10  #OS deaths, PFI events
MIN_RATIO        <- 0.05  #OS event ratio, PFI event ratio

a <- read.csv(file.path(DIR_RESULTS, "tcga_feasibility_audit.csv"),
              stringsAsFactors = FALSE)

#survival gate (all three required)
a$pass_samples  <- a$n_analysis  >= MIN_SAMPLES
a$pass_events   <- a$events_OS   >= MIN_EVENTS
a$pass_ratio    <- a$event_ratio >= MIN_RATIO
a$survival_pass <- a$pass_samples & a$pass_events & a$pass_ratio

#reasons for the excluded
reason <- function(r) {
  f <- c()
  if (!r$pass_samples) f <- c(f, sprintf("n<%d(%d)",  MIN_SAMPLES, r$n_analysis))
  if (!r$pass_events)  f <- c(f, sprintf("ev<%d(%d)", MIN_EVENTS,  r$events_OS))
  if (!r$pass_ratio)   f <- c(f, sprintf("ratio<%.2f(%.3f)", MIN_RATIO, r$event_ratio))
  paste(f, collapse = "; ")
}
a$reason <- vapply(seq_len(nrow(a)), function(i) reason(a[i, ]), character(1))

#save to .csv (survival_pass is written out under the column name "selected")
scored <- a[, c("cancer","n_analysis","events_OS","event_ratio",
                "survival_pass","reason")]
names(scored)[names(scored) == "survival_pass"] <- "selected"
write.csv(scored, file.path(DIR_RESULTS, "cancer_selection_reasons.csv"), row.names = FALSE)
final <- a$cancer[a$survival_pass]
writeLines(final, file.path(DIR_RESULTS, "final_cancer_list.txt"))

#report in console
cat(sprintf("GATE: n>=%d, events>=%d, ratio>=%.2f\n\n",
            MIN_SAMPLES, MIN_EVENTS, MIN_RATIO))
cat("== FINAL SELECTED (", length(final), ") ==\n", sep = "")
sel <- a[a$survival_pass, c("cancer","n_analysis","events_OS","event_ratio")]
print(sel, row.names = FALSE)
cat("\n== EXCLUDED (", sum(!a$survival_pass), ") ==\n", sep = "")
print(a[!a$survival_pass, c("cancer","n_analysis","events_OS","reason")], row.names = FALSE)
cat("\nFinal list saved to ", file.path(DIR_RESULTS, "final_cancer_list.txt"), "\n", sep = "")