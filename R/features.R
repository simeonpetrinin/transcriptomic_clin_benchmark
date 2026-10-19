#Features
#Leakage-safe encoding: every transform is fitted on the training rows only and then
#applied to the held-out rows. build_clinical() encodes a clinical subset, build_rna() the
#RNA block (once per fold), fuse() joins them with the penalty vector.

TOP_N_GENES          <- 500  #RNA block: top-variance genes (fit on train)
RNA_MIN_COMPLETENESS <- 0.80  #gene must be >=80% complete in the TRAIN fold

#fast column variances
.col_vars <- function(x) {
  
  #Use matrixStats for speed and numerical stability
  
  if (!requireNamespace("matrixStats", quietly = TRUE)) {
    stop("Package 'matrixStats' is required. Run renv::restore().")
  }
  matrixStats::colVars(x)
}

#clinical encoding: missing as its own level or a missing column; rare missing groups imputed
CLINICAL_MISSING  <- "category"  #"category" (main) | "impute" (sensitivity)
MISSING_LEVEL     <- "MISSING"
RARE_MISSING_VARS <- character(0)

#the variables of `vars` missing for 1 to n_min - 1 patients of `clin`
rare_missing_vars <- function(clin, vars, n_min = 6) {
  vars <- intersect(vars, names(clin))
  n_na <- vapply(vars, function(v) sum(is.na(clin[[v]])), numeric(1))
  vars[n_na > 0 & n_na < n_min]
}

build_clinical <- function(clinical_subset, tr, te){
  if (!CLINICAL_MISSING %in% c("category", "impute"))
    stop("CLINICAL_MISSING must be 'category' or 'impute'")
  as_cat <- identical(CLINICAL_MISSING, "category")
  cols_tr <- list()
  cols_te <- list()
  for (v in clinical_subset){
    xtr <- tr[[v]]; xte <- te[[v]]
    if (is.null(xtr) || is.null(xte)) stop("Missing clinical column: ", v)
    #A variable can be entirely missing in a smaller inner-training split.
    if (all(is.na(xtr))) next
    na_tr <- is.na(xtr); na_te <- is.na(xte)
    cat_v <- as_cat && !(v %in% RARE_MISSING_VARS)  #rare missing group: impute, no missing column
    
    if (is.numeric(xtr)){  #numeric: median-fill + z-score (+ indicator)
      med <- median(xtr, na.rm = TRUE)
      xtr[na_tr] <- med
      xte[na_te] <- med
      m <- mean(xtr)
      s <- sd(xtr)
      if (is.na(s) || s == 0) s <- 1
      mt <- matrix((xtr - m)/s, ncol = 1, dimnames = list(NULL, v))
      me <- matrix((xte - m)/s, ncol = 1, dimnames = list(NULL, v))
      if (cat_v && any(na_tr)) {  #indicator learned only if TRAIN has missings
        nm <- paste0(v, "_missing")
        mt <- cbind(mt, matrix(as.integer(na_tr), ncol = 1, dimnames = list(NULL, nm)))
        me <- cbind(me, matrix(as.integer(na_te), ncol = 1, dimnames = list(NULL, nm)))
      }
      cols_tr[[v]] <- mt
      cols_te[[v]] <- me
      
    }
    else{  #categorical: MISSING level (or mode) + one-hot
      xtr <- as.character(xtr); xte <- as.character(xte)
      mo <- names(sort(table(xtr[!na_tr]), decreasing = TRUE))[1]  #train mode (observed)
      if (cat_v && any(na_tr)) {
        xtr[na_tr] <- MISSING_LEVEL
        xte[na_te] <- MISSING_LEVEL
      } else {
        xtr[na_tr] <- mo
        xte[na_te] <- mo  #test-only missings -> TRAIN mode
      }
      levs <- sort(unique(xtr))
      unseen <- !(xte %in% levs)
      if (any(unseen)) xte[unseen] <- mo  #unseen TEST levels -> TRAIN mode
      if (length(levs) < 2) next  #dropped (constant -> contributes nothing)
      ref    <- mo
      nonref <- setdiff(levs, ref)
      mk <- function(x) {
        m <- sapply(nonref, function(l) as.integer(x == l))
        if (is.null(dim(m))) m <- matrix(m, nrow = length(x))  #ensure result is matrix
        colnames(m) <- paste0(v, "_", nonref)
        m
      }
      cols_tr[[v]] <- mk(xtr)
      cols_te[[v]] <- mk(xte)
    }
  }
  Xtr <- if (length(cols_tr)) do.call(cbind, cols_tr) else matrix(0, nrow(tr), 0)
  Xte <- if (length(cols_te)) do.call(cbind, cols_te) else matrix(0, nrow(te), 0)
  #record which variables were used or dropped in this split
  used <- names(cols_tr)
  list(train = Xtr, test = Xte,
       used = used, dropped = setdiff(clinical_subset, used))
}


#RNA block, fitted on training rows only: completeness, median fill, top-variance genes
build_rna <- function(tr_expr, te_expr, n_genes = TOP_N_GENES) {
  complete <- colMeans(!is.na(tr_expr))
  keep     <- which(complete >= RNA_MIN_COMPLETENESS)
  if (!length(keep))
    stop("No RNA genes pass the training completeness threshold.")
  
  Xtr0 <- tr_expr[, keep, drop = FALSE]
  Xte0 <- te_expr[, keep, drop = FALSE]
  
  med <- matrixStats::colMedians(Xtr0, na.rm = TRUE)
  miss_tr <- which(is.na(Xtr0), arr.ind = TRUE)
  miss_te <- which(is.na(Xte0), arr.ind = TRUE)
  if (nrow(miss_tr)) Xtr0[miss_tr] <- med[miss_tr[, 2]]
  if (nrow(miss_te)) Xte0[miss_te] <- med[miss_te[, 2]]
  
  v   <- .col_vars(Xtr0)  #variance per gene on TRAIN
  sel <- order(v, decreasing = TRUE)[seq_len(min(n_genes, ncol(Xtr0)))]
  Xtr <- Xtr0[, sel, drop = FALSE]
  Xte <- Xte0[, sel, drop = FALSE]
  mu  <- colMeans(Xtr)
  s <- apply(Xtr, 2, sd)
  s[is.na(s) | s == 0] <- 1
  list(train = scale(Xtr, center = mu, scale = s),
       test  = scale(Xte, center = mu, scale = s))
}

#RNA preprocessing cache for the shared inner folds
build_inner_rna <- function(tr_expr, foldid, n_genes = TOP_N_GENES) {
  if (length(foldid) != nrow(tr_expr) || anyNA(foldid) || any(foldid != as.integer(foldid)))
    stop("Inner fold IDs must identify every outer-training RNA row.")
  ids <- sort(unique(foldid))
  if (length(ids) < 3L || !identical(as.integer(ids), seq_along(ids)))
    stop("Inner fold IDs must be consecutive integers starting at 1 (at least 3 folds).")
  lapply(ids, function(k) {
    tri <- which(foldid != k); vai <- which(foldid == k)
    rn <- tryCatch(
      build_rna(tr_expr[tri, , drop = FALSE], tr_expr[vai, , drop = FALSE], n_genes),
      error = function(e) list(error = conditionMessage(e)))
    list(fold = k, train_idx = tri, valid_idx = vai, rna = rn)
  })
}

#fuse clinical + RNA and build the penalty-factor vector
fuse <- function(clin, rna, clin_penalty = 0) {
  list(
    train   = cbind(clin$train, rna$train),
    test    = cbind(clin$test,  rna$test),
    penalty = c(rep(clin_penalty, ncol(clin$train)), rep(1, ncol(rna$train))),
    block   = c(rep("clinical", ncol(clin$train)), rep("rna", ncol(rna$train)))
  )
}