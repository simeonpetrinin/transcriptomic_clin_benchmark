#Published-baseline search
#Collects candidate studies for the published-baseline figure, for screening against the inclusion rule
#(>= 2 TCGA cancers; clinical-only and clinical + RNA models compared per cancer; overall survival; named variables).
#Three routes, each logged with its date and count for the Methods:
#  1. forward citations of Herrmann et al. 2021, Zhao et al. 2015 and Yuan et al. 2014 (OpenAlex)
#  2. keyword search of titles and abstracts, 2011-2026 (OpenAlex)
#  3. the same keyword search in PubMed (E-utilities), abstracts then filled in from OpenAlex by PMID
#Writes data/processed/published/search/search_hits.csv (one row per study, all routes that found it)
#and search_log.txt. Needs only jsonlite and an internet connection; nothing in the pipeline is touched.
#Usage, from the project root: source("paths.R"); source("scripts/search_published_baselines.R"); run_search()

if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Package 'jsonlite' is missing.")

SEARCH_DIR <- file.path(DIR_PROCESSED, "published", "search")
OA_MAILTO  <- ""  #optional: an e-mail address here puts the requests in OpenAlex's faster "polite pool"
YEARS      <- c(2011, 2026)
SEEDS <- c(herrmann_2021 = "10.1093/bib/bbaa167",  #Brief Bioinform 22:bbaa167
           zhao_2015     = "10.1093/bib/bbu003",   #Brief Bioinform 16:291
           yuan_2014     = "10.1038/nbt.2940")     #Nat Biotechnol 32:644
#keyword search: TCGA x survival x clinical x RNA/omics x added value, titles and abstracts
PUBMED_QUERY <- paste(
  '("TCGA"[tiab] OR "The Cancer Genome Atlas"[tiab])',
  'AND (survival[tiab] OR prognos*[tiab])',
  'AND (clinical[tiab])',
  'AND ("gene expression"[tiab] OR transcriptom*[tiab] OR mRNA[tiab] OR "RNA-seq"[tiab] OR multi-omics[tiab] OR multiomics[tiab] OR "multi-omic"[tiab])',
  'AND ("added value"[tiab] OR "additional value"[tiab] OR incremental[tiab] OR "clinical model"[tiab] OR "clinical variables"[tiab] OR "clinical covariates"[tiab] OR "clinical data alone"[tiab] OR "clinical features"[tiab])')
OA_QUERY <- paste(
  '("TCGA" OR "The Cancer Genome Atlas") AND (survival OR prognosis OR prognostic) AND clinical',
  'AND ("gene expression" OR transcriptome OR transcriptomic OR mRNA OR "RNA-seq" OR "multi-omics" OR multiomics)',
  'AND ("added value" OR "additional value" OR incremental OR "clinical model" OR "clinical variables"',
  'OR "clinical covariates" OR "clinical data alone" OR "clinical features")')
OA_CAP <- 5000  #safety cap per route

.enc <- function(x) utils::URLencode(x, reserved = TRUE)
.get <- function(url) {
  for (i in 1:3) {  #retry a busy server twice
    x <- tryCatch(jsonlite::fromJSON(url, simplifyVector = FALSE), error = function(e) e)
    if (!inherits(x, "error")) return(x)
    Sys.sleep(2 * i)
  }
  stop("request failed: ", url, "\n  ", conditionMessage(x))
}
.mail <- function() if (nzchar(OA_MAILTO)) paste0("&mailto=", .enc(OA_MAILTO)) else ""

#OpenAlex stores abstracts as word -> positions; rebuild the text
.abstract <- function(inv) {
  if (is.null(inv) || !length(inv)) return(NA_character_)
  pos <- unlist(lapply(names(inv), function(w) setNames(unlist(inv[[w]]), rep(w, length(inv[[w]])))))
  paste(names(sort(pos)), collapse = " ")
}

#one OpenAlex work -> one row
.row <- function(w, route) {
  ids <- w$ids
  data.frame(route = route, openalex_id = sub("^https://openalex.org/", "", w$id),
             doi = if (is.null(w$doi)) NA_character_ else sub("^https://doi.org/", "", w$doi),
             pmid = if (is.null(ids$pmid)) NA_character_ else sub("^https://pubmed.ncbi.nlm.nih.gov/", "", ids$pmid),
             title = if (is.null(w$title)) NA_character_ else w$title,
             year = if (is.null(w$publication_year)) NA_integer_ else w$publication_year,
             venue = if (is.null(w$primary_location$source$display_name)) NA_character_ else w$primary_location$source$display_name,
             type = if (is.null(w$type)) NA_character_ else w$type,
             abstract = .abstract(w$abstract_inverted_index), stringsAsFactors = FALSE)
}

#every page of an OpenAlex works query (cursor paging)
.oa_all <- function(filter_or_search, route) {
  rows <- list(); cursor <- "*"; n <- 0
  repeat {
    x <- .get(sprintf("https://api.openalex.org/works?%s&per-page=200&cursor=%s%s", filter_or_search, .enc(cursor), .mail()))
    if (!length(x$results)) break
    rows <- c(rows, lapply(x$results, .row, route = route)); n <- n + length(x$results)
    cursor <- x$meta$next_cursor
    if (is.null(cursor) || n >= OA_CAP) break
    Sys.sleep(0.15)
  }
  if (length(rows)) do.call(rbind, rows) else NULL
}

#route 1: forward citations of each seed paper
search_citations <- function() {
  do.call(rbind, lapply(names(SEEDS), function(s) {
    w <- .get(sprintf("https://api.openalex.org/works/doi:%s%s", SEEDS[[s]], sub("^&", "?", .mail())))
    id <- sub("^https://openalex.org/", "", w$id)
    cat(sprintf("[search] %-14s %s, cited by %d (OpenAlex)\n", s, id, w$cited_by_count))
    .oa_all(sprintf("filter=cites:%s,publication_year:%d-%d", id, YEARS[1], YEARS[2]), paste0("cites_", s))
  }))
}

#route 2: OpenAlex keyword search of titles and abstracts
search_openalex <- function() {
  .oa_all(sprintf("filter=title_and_abstract.search:%s,publication_year:%d-%d", .enc(OA_QUERY), YEARS[1], YEARS[2]),
          "openalex_keywords")
}

#route 3: PubMed keyword search; details and abstracts looked up in OpenAlex by PMID (50 per request)
search_pubmed <- function() {
  u <- sprintf("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmode=json&retmax=10000&datetype=pdat&mindate=%d&maxdate=%d&term=%s",
               YEARS[1], YEARS[2], .enc(PUBMED_QUERY))
  pm <- unlist(.get(u)$esearchresult$idlist)
  cat(sprintf("[search] PubMed: %d records\n", length(pm)))
  if (!length(pm)) return(NULL)
  chunks <- split(pm, ceiling(seq_along(pm) / 50))
  rows <- do.call(rbind, lapply(chunks, function(ch) {
    Sys.sleep(0.15)
    x <- .get(sprintf("https://api.openalex.org/works?filter=pmid:%s&per-page=50%s", paste(ch, collapse = "|"), .mail()))
    if (length(x$results)) do.call(rbind, lapply(x$results, .row, route = "pubmed_keywords")) else NULL
  }))
  miss <- setdiff(pm, rows$pmid)  #PubMed records OpenAlex does not have: kept with the PMID only
  if (length(miss)) rows <- rbind(rows, data.frame(route = "pubmed_keywords", openalex_id = NA, doi = NA, pmid = miss,
                                                   title = NA, year = NA, venue = NA, type = NA, abstract = NA))
  rows
}

#all three routes, merged to one row per study with every route that found it
run_search <- function() {
  dir.create(SEARCH_DIR, showWarnings = FALSE, recursive = TRUE)
  t0 <- format(Sys.time(), "%Y-%m-%d %H:%M")
  parts <- list(citations = search_citations(), openalex = search_openalex(), pubmed = search_pubmed())
  all <- do.call(rbind, parts)
  key <- ifelse(!is.na(all$openalex_id), all$openalex_id, paste0("pmid:", all$pmid))
  routes <- tapply(all$route, key, function(r) paste(sort(unique(r)), collapse = ";"))
  first <- all[!duplicated(key), ]; first$key <- key[!duplicated(key)]
  first$routes <- routes[first$key]; first$route <- NULL
  first <- first[order(-as.numeric(first$year), first$title), ]
  write.csv(first, file.path(SEARCH_DIR, "search_hits.csv"), row.names = FALSE)
  counts <- table(all$route)
  log <- c(sprintf("search run: %s", t0), sprintf("years: %d-%d", YEARS[1], YEARS[2]),
           sprintf("seeds: %s", paste(sprintf("%s (doi %s)", names(SEEDS), SEEDS), collapse = "; ")),
           sprintf("PubMed query: %s", PUBMED_QUERY), sprintf("OpenAlex query (title_and_abstract): %s", OA_QUERY),
           sprintf("records per route: %s", paste(sprintf("%s %d", names(counts), as.integer(counts)), collapse = ", ")),
           sprintf("unique studies after merging: %d", nrow(first)))
  writeLines(log, file.path(SEARCH_DIR, "search_log.txt"))
  cat(paste(log, collapse = "\n"), "\n")
  invisible(first)
}
