#Enumerate
#CANCER_VARS lists each cancer's main-set clinical variables: standard variables and
#AJCC 8th tier 1 items (BCLC 2022 for LIHC) that pass the eligibility filters. This
#file also builds the subset lattice from that list.

ANCHOR_LABEL <- "(omics_only)"  #the empty-clinical (omics-only) anchor

CANCER_VARS <- list(
  
  #12 variables, 4,095 subsets
  ESCA = list(age = "numeric", sex = "categorical", race = "categorical",
              stage = "categorical", T = "categorical", N = "categorical",
              M = "categorical", histology = "categorical", grade = "categorical",
              smoking = "categorical", esoph_loc = "categorical", nodes_pos = "numeric"),
  
  #11 variables, 2,047 subsets
  HNSC = list(age = "numeric", sex = "categorical", race = "categorical",
              stage = "categorical", T = "categorical", N = "categorical",
              grade = "categorical", smoking = "categorical",
              hn_site = "categorical", ece = "categorical", nodes_pos = "numeric"),
  #PAAD: size defines T1-T3 and the positive-node count defines N (1b)
  PAAD = list(age = "numeric", sex = "categorical", race = "categorical",
              stage = "categorical", T = "categorical", N = "categorical",
              histology = "categorical", grade = "categorical", smoking = "categorical",
              max_size = "numeric", nodes_pos = "numeric"),
  
  #10 variables, 1,023 subsets
  #BRCA: ER, PR and HER2 are required for stage grouping (1a)
  BRCA = list(age = "numeric", race = "categorical",
              stage = "categorical", T = "categorical", N = "categorical",
              histology = "categorical", ER = "categorical", PR = "categorical",
              HER2 = "categorical", nodes_pos = "numeric"),
  #COAD: tumour deposits (N1c) and the positive-node count (N1/N2) are part of N (1b)
  COAD = list(age = "numeric", sex = "categorical", race = "categorical",
              stage = "categorical", T = "categorical", N = "categorical",
              M = "categorical", histology = "categorical",
              tumour_deposits = "categorical", nodes_pos = "numeric"),
  #STAD: N is defined by the number of positive nodes (1b).
  STAD = list(age = "numeric", sex = "categorical", race = "categorical",
              stage = "categorical", T = "categorical", N = "categorical",
              M = "categorical", histology = "categorical", grade = "categorical",
              nodes_pos = "numeric"),
  
  #9 variables, 511 subsets
  #BLCA: grade fails the 95% rule on the modelling cohort
  BLCA = list(age = "numeric", sex = "categorical", race = "categorical",
              stage = "categorical", T = "categorical", N = "categorical",
              smoking = "categorical", papillary = "categorical", nodes_pos = "numeric"),
  
  #8 variables, 255 subsets
  CESC = list(age = "numeric", race = "categorical", stage = "categorical",
              T = "categorical", N = "categorical", histology = "categorical",
              grade = "categorical", smoking = "categorical"),
  #KIRC and KIRP: nodes are rarely assessed in kidney surgery, so NX is a large level
  KIRC = list(age = "numeric", sex = "categorical", race = "categorical",
              stage = "categorical", T = "categorical", N = "categorical",
              M = "categorical", grade = "categorical"),
  #KIRP: type 1 vs type 2 is the histology subtype field (i).
  KIRP = list(age = "numeric", sex = "categorical", race = "categorical",
              stage = "categorical", T = "categorical", N = "categorical",
              smoking = "categorical", tumor_type = "categorical"),
  #LIHC: vascular invasion is part of T1/T2 (1b)
  LIHC = list(age = "numeric", sex = "categorical", race = "categorical",
              stage = "categorical", T = "categorical", grade = "categorical",
              vascular_inv = "categorical", ecog = "numeric"),
  READ = list(age = "numeric", sex = "categorical",
              stage = "categorical", T = "categorical", N = "categorical",
              M = "categorical", histology = "categorical", nodes_pos = "numeric"),
  #SKCM: Breslow thickness and ulceration define T (1b).
  SKCM = list(age = "numeric", sex = "categorical",
              stage = "categorical", T = "categorical", N = "categorical",
              M = "categorical", breslow = "numeric", ulceration = "categorical"),
  #UCEC: myometrial invasion (T1a/T1b) and pelvic/para-aortic positive nodes (N1/N2) (1b)
  UCEC = list(age = "numeric", race = "categorical", stage = "categorical",
              histology = "categorical", grade = "categorical",
              myo_inv = "numeric", pelv_pos = "numeric", aor_pos = "numeric"),
  
  #7 variables, 127 subsets
  LUAD = list(age = "numeric", sex = "categorical", race = "categorical",
              stage = "categorical", T = "categorical", N = "categorical",
              smoking = "categorical"),
  LUSC = list(age = "numeric", sex = "categorical", race = "categorical",
              stage = "categorical", T = "categorical", N = "categorical",
              smoking = "categorical"),
  #SARC: no stage or grade in TCGA; site, size and necrosis are AJCC items (1b/1c)
  SARC = list(age = "numeric", sex = "categorical", race = "categorical",
              histology = "categorical", site = "categorical",
              tumour_length = "numeric", necrosis = "categorical"),
  
  #5 variables, 31 subsets
  LGG  = list(age = "numeric", sex = "categorical", race = "categorical",
              histology = "categorical", grade = "categorical"),
  #OV: one vs both ovaries is part of the T1a/T1b definition (1b).
  OV   = list(age = "numeric", race = "categorical", stage = "categorical",
              grade = "categorical", bilateral = "categorical"),
  
  #3 variables, 7 subsets
  GBM  = list(age = "numeric", sex = "categorical", race = "categorical"),
  LAML = list(age = "numeric", sex = "categorical", race = "categorical")
)

#all non-empty subsets of a character vector (the power set)
all_subsets <- function(vars) {
  unlist(lapply(seq_along(vars), function(k) combn(vars, k, simplify = FALSE)),
         recursive = FALSE)
}

#the subset index for one cancer (anchor prepended as subset_id 0)
build_subset_index <- function(code, include_anchor = TRUE) {
  if (is.null(CANCER_VARS[[code]]))
    stop("no clinical variables defined for '", code, "' in CANCER_VARS")
  vars <- names(CANCER_VARS[[code]])
  subs <- all_subsets(vars)
  idx <- data.frame(subset_id = seq_along(subs),
                    n_vars    = vapply(subs, length, integer(1)),
                    variables = vapply(subs, paste, character(1), collapse = "+"),
                    stringsAsFactors = FALSE)
  if (include_anchor)
    idx <- rbind(data.frame(subset_id = 0L, n_vars = 0L, variables = ANCHOR_LABEL,
                            stringsAsFactors = FALSE), idx)
  idx
}

#the compute bill
lattice_sizes <- function(codes = names(CANCER_VARS)) {
  d <- data.frame(cancer  = codes,
                  k       = vapply(codes, function(c) length(CANCER_VARS[[c]]), integer(1)),
                  subsets = vapply(codes, function(c) 2^length(CANCER_VARS[[c]]) - 1, numeric(1)),
                  stringsAsFactors = FALSE)
  d <- d[order(-d$subsets), ]
  cat("=== lattice size per cancer (2^k - 1 subsets, + 1 anchor each) ===\n")
  print(d, row.names = FALSE)
  cat(sprintf("\nTOTAL clinical subsets: %s\n", format(sum(d$subsets), big.mark = ",")))
  cat(sprintf("\nTOTAL specifications including anchors: %s\n", format(sum(d$subsets) + length(codes), big.mark = ",")))
  invisible(d)
}