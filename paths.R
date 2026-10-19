#Paths
#Every directory and file path in the project, the shared settings from settings.txt
#(N_WORKERS), and source_library(), which sources the R/ library. Assumes the working
#directory is the project root.

DIR_R          <- "R"
DIR_SCRIPTS    <- "scripts"
DIR_RAW        <- file.path("data", "raw")  #never modified, never regenerated
DIR_PROCESSED  <- file.path("data", "processed")  #built once per cancer
DIR_CHECKPOINTS<- file.path("data", "checkpoints")  #transient, safe to delete
DIR_RESULTS    <- "results"  #scientific artifacts
DIR_FIGURES    <- "figures"  #manuscript-ready deliverables

#raw input files
CDR_PATH    <- file.path(DIR_RAW, "TCGA-CDR.xlsx")
PANCAN_PATH <- file.path(DIR_RAW, "clinical_PANCAN_patient_with_followup.tsv")
RNA_PATH    <- file.path(DIR_RAW,
                         "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz")

#per-cancer derived paths (helpers, so no script builds a filename)
path_modelling <- function(code) file.path(DIR_PROCESSED, paste0(code, "_modelling.rds"))
path_subsets   <- function(code) file.path(DIR_PROCESSED, paste0(code, "_subsets.csv"))
path_results   <- function(code) file.path(DIR_RESULTS,   paste0(code, "_results.csv"))
path_manifest  <- function(code) file.path(DIR_RESULTS,   paste0(code, "_manifest.txt"))
path_summary   <- function(code) file.path(DIR_RESULTS,   paste0(code, "_subset_summary.csv"))
path_timing    <- function(code) file.path(DIR_RESULTS,   paste0(code, "_timing_summary.csv"))
path_figure    <- function(code, what) file.path(DIR_FIGURES, sprintf("%s_%s.png", code, what))
path_checkpoint<- function(code, r) file.path(DIR_CHECKPOINTS,
                                              sprintf("%s_repeat_%02d.rds", code, r))

#create output dirs if absent (raw/processed must already exist)
for (d in c(DIR_CHECKPOINTS, DIR_RESULTS, DIR_FIGURES)) {
  if(!dir.exists(d)){
    dir.create(d, recursive=TRUE)  #Will output serious warnings, avoids 'directory already exists'
  }
}

#shared settings: settings.txt holds "name = value" lines, so the core count is set in one place
SETTINGS_PATH <- "settings.txt"
read_settings <- function(path = SETTINGS_PATH) {
  if (!file.exists(path)) stop("settings file not found: ", path)
  l <- trimws(sub("#.*$", "", readLines(path, warn = FALSE)))
  l <- l[nzchar(l)]
  if (!all(grepl("=", l, fixed = TRUE))) stop("every line of ", path, " must be name = value")
  setNames(as.list(trimws(sub("^[^=]*=", "", l))), trimws(sub("=.*$", "", l)))
}
SETTINGS  <- read_settings()
N_WORKERS <- suppressWarnings(as.integer(SETTINGS$workers))
if (length(N_WORKERS) != 1 || is.na(N_WORKERS) || N_WORKERS < 1)
  stop("settings.txt: 'workers' must be a whole number of at least 1")
if (N_WORKERS > parallel::detectCores(logical = TRUE))
  warning("settings.txt asks for ", N_WORKERS, " workers but this machine has ",
          parallel::detectCores(logical = TRUE), " logical cores")

#source the whole library
source_library <- function(){
  for (f in c("features.R", "metrics.R", "enumerate.R", "models.R",
              "load_data.R", "cv.R", "extended.R",
              "sweep_genes.R", "sweep_pfi.R", "sweep_missing.R", "sweeps.R")){
    source(file.path(DIR_R, f))
  }
  invisible(TRUE)
}