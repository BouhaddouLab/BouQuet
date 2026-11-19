# =============================================================================
# BouQuet: Bouhaddou Lab Quantitative proteomics Toolkit - all package functions
#
# Every function in this file is available after library(BouQuet).
# To add a function: write it in this file (or a new .R file in R/), then reinstall the package.
#
#   1. Setup              check_packages()
#   2. Databases          ptm_species_info(), db_dir(), db_needs_download(), db_download(), db_log(),
#                         get_uniprot_genemap(), get_msigdb(), get_psp_regulatory_sites()
#   3. MSstats helpers    Check.MsStats.DataProcess.Output(), Make.LFQ.Contrast.Matrix(), Save.Csv()
#   4. Gene IDs           get_orgdb(), convert_geneids()
#   5. Enrichment (D2)    do_enrichmap(), simplifyEnrichBySimilarUniverseMembership(),
#                         fixMsigdbGONames(), fixMsigdbPNames(), heatmapNumbered(),
#                         enrichHeatmapBestPerGroup()
#   6. Kinase activity (D3)  do_KinaseActivity_Network(), kinaseActivity(),
#                            loadKinaseDataOmniPath(), loadKinaseData()
#   7. Built-in tables    bouquet_file()
# =============================================================================


# =============================================================================
# 1. Setup
# =============================================================================

# Everything BouQuet needs is installed together with it. This only checks the optional extras
# (data for msigdbr, mouse/rat gene annotations, reading .parquet files).
check_packages <- function(){
  cran <- c("msigdbdf", "arrow")
  bioc <- c("org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db", "AnnotationDbi", "biomaRt")
  missing_cran <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
  missing_bioc <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_cran) + length(missing_bioc) == 0){
    message("All optional packages are installed.")
  } else {
    message("Optional packages not installed (only needed for some steps). To install them, run:")
    if (length(missing_cran)) message('  install.packages(c("', paste(missing_cran, collapse = '", "'), '"), repos = c("https://igordot.r-universe.dev", "https://cloud.r-project.org"))')
    if (length(missing_bioc)) message('  BiocManager::install(c("', paste(missing_bioc, collapse = '", "'), '"))')
  }
  invisible(c(missing_cran, missing_bioc))
}


# =============================================================================
# 2. Reference databases (shared download cache, version log)
# =============================================================================

# Fetches every reference database the pipeline uses.
# Nothing needs to live in a per-project "databases/" folder anymore.
#
# How it works
#   * Files are downloaded from the source web portals and stored in ONE shared
#     cache folder (by default a per-user folder on this computer), so ten
#     projects share one copy instead of keeping ten copies.
#   * A cached file is re-downloaded when it is older than `ptm.db_max_age_days`,
#     or on every run if `ptm.db_refresh = TRUE`.
#   * If a download fails (offline, server down) and a cached copy exists, the
#     cached copy is used with a warning.
#   * Every database used is written to results/database_versions.tsv, so each
#     analysis records exactly which versions and dates it was run with.
#
# Settings (set in the Setup chunk of the Rmd with options()):
#   ptm.db_dir          where downloads are stored (default: shared user cache;
#                       use tempdir() to keep nothing after R closes)
#   ptm.db_refresh      TRUE = re-download everything each run
#   ptm.db_max_age_days re-download cached files older than this
#   ptm.results_dir     where database_versions.tsv is written
#   ptm.psp_file        path to the PhosphoSitePlus Regulatory_sites file
#                       (manual download, see get_psp_regulatory_sites())
# =============================================================================

# ---- Species information ----------------------------------------------------

ptm_species_info <- function(species){
  species <- tolower(species)
  info <- list(
    human = list(taxid = 9606L,  proteome = "UP000005640", ensembl = "hsapiens",
                 msigdbr = "Homo sapiens"),
    mouse = list(taxid = 10090L, proteome = "UP000000589", ensembl = "mmusculus",
                 msigdbr = "Mus musculus"),
    rat   = list(taxid = 10116L, proteome = "UP000002494", ensembl = "rnorvegicus",
                 msigdbr = "Rattus norvegicus")
  )
  if (!species %in% names(info))
    stop("Species '", species, "' not supported. Choose one of: ", paste(names(info), collapse = ", "))
  info[[species]]
}

# ---- Cache folder and download helper ----------------------------------------

db_dir <- function(){
  d <- getOption("ptm.db_dir", tools::R_user_dir("BouQuet", which = "cache"))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  d
}

db_needs_download <- function(path, refresh = getOption("ptm.db_refresh", FALSE),
                              max_age_days = getOption("ptm.db_max_age_days", 30)){
  if (refresh || !file.exists(path)) return(TRUE)
  age <- as.numeric(difftime(Sys.time(), file.mtime(path), units = "days"))
  if (!is.null(max_age_days) && is.finite(max_age_days) && age > max_age_days){
    message(basename(path), " is ", round(age), " days old; downloading a fresh copy.")
    return(TRUE)
  }
  FALSE
}

db_download <- function(url, filename, ...){
  path <- file.path(db_dir(), filename)
  if (db_needs_download(path, ...)){
    message("Downloading ", filename, " ...")
    tmp <- paste0(path, ".part")
    old <- options(timeout = max(3600, getOption("timeout")))
    on.exit(options(old), add = TRUE)
    ok <- tryCatch({
      suppressWarnings(utils::download.file(url, tmp, mode = "wb", quiet = TRUE))
      file.exists(tmp) && file.size(tmp) > 0
    }, error = function(e){ message("  Download failed: ", conditionMessage(e)); FALSE })
    if (ok){
      file.rename(tmp, path)
    } else {
      unlink(tmp)
      if (file.exists(path)){
        warning("Could not download ", filename, "; using cached copy from ",
                format(file.mtime(path), "%Y-%m-%d"), call. = FALSE)
      } else {
        stop("Could not download ", filename, " from ", url,
             " and no cached copy exists. Check your internet connection.", call. = FALSE)
      }
    }
  }
  path
}

# Record which databases (and versions) were used in this analysis
db_log <- function(name, source, path = NA_character_, version = NA_character_){
  out_dir <- getOption("ptm.results_dir", "results")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  log_file <- file.path(out_dir, "database_versions.tsv")
  retrieved <- if (!is.na(path) && file.exists(path)) format(file.mtime(path), "%Y-%m-%d") else format(Sys.Date())
  entry <- data.table(database = name, source = source, version = version,
                      retrieved = retrieved, local_path = path, used_on = format(Sys.Date()))
  if (file.exists(log_file)){
    log <- fread(log_file, colClasses = "character")
    log <- rbind(log[database != name], entry, fill = TRUE)
  } else log <- entry
  fwrite(log, log_file, sep = "\t")
  invisible(entry)
}

# ---- UniProt ----------------------------------------------------------------

# UniProt accession -> primary gene symbol table (replaces the 38 MB idmapping file)
get_uniprot_genemap <- function(species = "human", ...){
  s <- ptm_species_info(species)
  url <- paste0("https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=tsv",
                "&fields=accession%2Cgene_primary&query=%28organism_id%3A", s$taxid, "%29")
  path <- db_download(url, paste0("uniprot_genemap_", s$taxid, ".tsv.gz"), ...)
  db_log(paste0("UniProt accession-to-gene map (", species, ")"), url, path)
  gm <- fread(path, header = TRUE, sep = "\t", quote = "")
  setnames(gm, c("uniprot", "gene"))
  gm[gene != ""]
}

# ---- MSigDB (via the msigdbr package; no GMT files needed) --------------------

# collection: "C5" + subcollection "GO:BP"  -> GO biological process
#             "C2" + subcollection "CP"     -> all canonical pathways (Reactome, KEGG, WikiPathways...)
# gmt_file:   optionally give a local .gmt file instead (e.g. to reproduce an older analysis)
get_msigdb <- function(collection = "C5", subcollection = "GO:BP", species = "human", gmt_file = NULL){
  if (!is.null(gmt_file)){
    gmt <- setDT(clusterProfiler::read.gmt(gmt_file))
    db_log(paste0("MSigDB ", collection, " ", subcollection), "local GMT file", gmt_file)
    return(gmt[, .(term = as.character(term), gene = as.character(gene))])
  }
  if (!requireNamespace("msigdbr", quietly = TRUE))
    stop("Please install msigdbr: install.packages('msigdbr')", call. = FALSE)
  sp <- ptm_species_info(species)$msigdbr
  new_api <- "collection" %in% names(formals(msigdbr::msigdbr))   # msigdbr >= 10 renamed the arguments
  m <- tryCatch({
    if (new_api) msigdbr::msigdbr(species = sp, collection = collection)
    else         msigdbr::msigdbr(species = sp, category  = collection)
  }, error = function(e){
    stop("msigdbr failed: ", conditionMessage(e),
         "\nIf it asks for 'msigdbdf', install it with:\n",
         "  install.packages('msigdbdf', repos = c('https://igordot.r-universe.dev', 'https://cloud.r-project.org'))",
         call. = FALSE)
  })
  setDT(m)
  subcol <- intersect(c("gs_subcollection", "gs_subcat"), names(m))[1]
  # "CP" keeps every CP:* sub-collection, matching the c2.cp.*.gmt file from the MSigDB website
  m <- m[m[[subcol]] == subcollection | startsWith(m[[subcol]], paste0(subcollection, ":"))]
  if (nrow(m) == 0) stop("No gene sets found for ", collection, " / ", subcollection, call. = FALSE)
  ver <- paste0("msigdbr ", utils::packageVersion("msigdbr"),
                if ("db_version" %in% names(m)) paste0(", MSigDB ", m$db_version[1]) else "")
  db_log(paste0("MSigDB ", collection, " ", subcollection), "msigdbr R package", NA_character_, ver)
  unique(m[, .(term = gs_name, gene = gene_symbol)])
}

# ---- PhosphoSitePlus (manual download: requires free registration + license) --

# Which file is used, in order:
#   1. the file set as ptm.psp_file (psp_file in Setup)
#   2. the newest "Regulatory_sites..." file in the project's resources/ folder
#   3. the newest one in the shared cache folder (db_dir())
#   4. the copy that comes with BouQuet (inst/extdata in the package)
# The file can be the original text download (.txt / .gz), a .tsv.gz table, or a .parquet copy.
get_psp_regulatory_sites <- function(path = getOption("ptm.psp_file", NULL),
                                     search_dirs = c("resources", db_dir(), system.file("extdata", package = "BouQuet"))){
  if (is.null(path)){
    for (d in search_dirs[nzchar(search_dirs) & dir.exists(search_dirs)]){
      hits <- list.files(d, pattern = "Regulatory_sites", full.names = TRUE)
      if (length(hits) > 0){ path <- hits[which.max(file.mtime(hits))]; break }
    }
  }
  if (is.null(path) || !file.exists(path))
    stop("PhosphoSitePlus Regulatory_sites file not found.\n",
         "PSP cannot be downloaded automatically (it needs a free login and non-commercial license).\n",
         "  1. Log in at https://www.phosphosite.org and open 'Downloads'\n",
         "  2. Download 'Regulatory_sites.gz'\n",
         "  3. Put it (gzipped is fine) in the project's resources/ folder or in: ", db_dir(), "\n",
         "     or set options(ptm.psp_file = 'path/to/Regulatory_sites') in the Setup chunk.", call. = FALSE)
  age <- as.numeric(difftime(Sys.time(), file.mtime(path), units = "days"))
  file_date <- regmatches(basename(path), regexpr("[0-9]{4}-[0-9]{2}-[0-9]{2}", basename(path)))   # date in the file name, if any
  if (length(file_date)) age <- as.numeric(Sys.Date() - as.Date(file_date))
  message("Using PhosphoSitePlus file: ", path)
  if (age > 180) message("Note: your PhosphoSitePlus file is ", round(age), " days old. Consider downloading a new one.")
  db_log("PhosphoSitePlus Regulatory_sites", "https://www.phosphosite.org (manual download)", path)
  # .parquet (smaller lab copy in resources/) needs the arrow package; the original text file is read with fread
  if (grepl("\\.parquet$", path)){
    if (!requireNamespace("arrow", quietly = TRUE))
      stop("Reading ", basename(path), " needs the 'arrow' package: install.packages('arrow')", call. = FALSE)
    return(as.data.table(arrow::read_parquet(path)))
  }
  fread(path, skip = "ACC_ID")   # skips the license text at the top of the file
}


# =============================================================================
# 3. MSstats helpers (contrast matrix, saving tables)
# =============================================================================

Check.MsStats.DataProcess.Output <- function (input.data.frame, type = ""){
  if (type == "") {
    if (!"ProteinLevelData" %in% names(input.data.frame) || !"FeatureLevelData" %in% names(input.data.frame) ){
      stop("Error: Input dataframe list does not include ProteinLevelData and FeatureLevelData items.\n\tThis variable should come directly from the `dataProcess` function.")
    }
  } else if (type == "PTM") {
    if (!"PTM" %in% names(input.data.frame) || !"PROTEIN" %in% names(input.data.frame) ){
      stop("Error: Input dataframe list does not include PTM and PROTEIN items.\n\tThis variable should come directly from the `dataSummarizationPTM` function.\n\tIf this is not a PTM experiment, you should be using the standard MsStats workflow.")
    }
  } else {
    stop("DEV ERROR: Wrong type qualifier for wrapper function. Users should not see this error")
  }
}

Make.LFQ.Contrast.Matrix <- function(input.data.frame, condition.vector, exact.match = TRUE){
  Check.MsStats.DataProcess.Output(input.data.frame)
  columnNames <- as.character(levels(input.data.frame$ProteinLevelData$GROUP))
  
  positives = c()
  negatives = c()
  for (i in columnNames){
    for (j in columnNames){
      if (i != j){
        positives <- c(positives,i)
        negatives <- c(negatives,j)
      }
    }
  }
  contrasts <- data.table (positives, negatives)
  contrasts[,name := paste(positives,negatives, sep="-")]
  
  contrastNames <- character(0)
  for (re in condition.vector){
    if (exact.match) {
      contrastNames <- union(contrastNames, grep(paste0("^", re, "$"), contrasts$name, value = TRUE))
    } else {
      contrastNames <- union(contrastNames, grep(re, contrasts$name, value = TRUE))
    }
    
    
  }
  print (contrastNames)
  contrasts <- contrasts[name %in% contrastNames]
  print (contrasts)
  
  contrastMat <- matrix(rep(0,length(columnNames) * nrow(contrasts)), nrow=nrow(contrasts), ncol=length(columnNames),
                        dimnames = list(contrasts$name, columnNames))
  
  for (i in seq_len(nrow(contrasts))){
    #rows are named according to the positive condition:
    contrastMat[contrasts$name[i],contrasts$negatives[i]] = -1
    contrastMat[contrasts$name[i],contrasts$positives[i]] =  1
  }
  return (contrastMat)  
}

Save.Csv <- function(export.data, filename, base.dir = getwd()) {
  data.filename <- filename
  full.path <- file.path(base.dir, data.filename)
  data.table::fwrite(export.data, full.path)
}


# =============================================================================
# 4. Gene IDs: UniProt -> gene symbol
# =============================================================================

# convert_geneids(data_in, from_type, to_type)
#   1) maps with the species' OrgDb, if that package is installed (chosen by get_orgdb())
#   2) maps everything else UNIPROT <-> SYMBOL with UniProt's own accession-to-gene table:
#      downloaded the first time it is needed, then reused from the cache folder
#      (download a fresh copy with db_refresh = TRUE in Setup)
#   If several matches exist for one ID, the first is kept.

# Map a species name to its Bioconductor OrgDb package and load it.
# OrgDb packages are optional: if it is not installed, NULL is returned and gene symbols come from
# UniProt's accession-to-gene table instead (downloaded once into the cache folder, then reused).
get_orgdb <- function(species){
  orgdb_pkgs <- c(
    human = "org.Hs.eg.db",     "homo sapiens"            = "org.Hs.eg.db",
    mouse = "org.Mm.eg.db",     "mus musculus"            = "org.Mm.eg.db",
    rat = "org.Rn.eg.db",       "rattus norvegicus"       = "org.Rn.eg.db",
    zebrafish = "org.Dr.eg.db", "danio rerio"             = "org.Dr.eg.db",
    fly = "org.Dm.eg.db",       "drosophila melanogaster" = "org.Dm.eg.db",
    worm = "org.Ce.eg.db",      "caenorhabditis elegans"  = "org.Ce.eg.db",
    chicken = "org.Gg.eg.db",   "gallus gallus"           = "org.Gg.eg.db",
    cow = "org.Bt.eg.db",       "bos taurus"              = "org.Bt.eg.db",
    pig = "org.Ss.eg.db",       "sus scrofa"              = "org.Ss.eg.db",
    dog = "org.Cf.eg.db",       "canis familiaris"        = "org.Cf.eg.db",
    rhesus = "org.Mmu.eg.db",   "macaca mulatta"          = "org.Mmu.eg.db"
  )
  pkg <- orgdb_pkgs[tolower(trimws(species))]
  if (is.na(pkg)){
    message("No OrgDb known for species '", species, "'; using UniProt mapping only.")
    return(NULL)
  }
  if (!requireNamespace(pkg, quietly = TRUE)){
    message(pkg, " is not installed; using UniProt's gene names (downloaded once, then reused from the cache).")
    return(NULL)
  }
  getExportedValue(pkg, pkg)
}


convert_geneids <- function(data_in, from_type, to_type, strip_version = TRUE,
                            orgdb = NULL, species = "human", ...){
  if (length(to_type) != 1) stop("to_type must be a single column name, e.g. 'SYMBOL'", call. = FALSE)
  ids_raw <- as.character(data_in)

  # Clean IDs for lookup (the originals are not changed)
  query <- ids_raw
  if (strip_version && from_type == "ENSEMBL") query <- sub("\\..*$", "", query)
  if (strip_version && from_type == "UNIPROT") query <- sub("-[0-9]+$", "", query)
  keys <- unique(query[!is.na(query) & query != ""])

  # 1) OrgDb (chosen from species unless supplied explicitly)
  result <- rep(NA_character_, length(ids_raw))
  if (is.null(orgdb)) orgdb <- get_orgdb(species)
  if (!is.null(orgdb) && length(keys) > 0){
    if (!requireNamespace("AnnotationDbi", quietly = TRUE))
      stop("Please install AnnotationDbi: BiocManager::install('AnnotationDbi')", call. = FALSE)
    mapped <- tryCatch(
      suppressMessages(AnnotationDbi::select(orgdb, keys = keys, keytype = from_type, columns = to_type)),
      error = function(e){
        # select() errors when *none* of the keys are valid (or the keytype is absent); fall through to UniProt
        message("OrgDb lookup found no matches (", conditionMessage(e), ")"); NULL
      })
    if (!is.null(mapped)){
      mapped <- mapped[!is.na(mapped[[to_type]]) & !duplicated(mapped[[from_type]]), ]
      result <- mapped[[to_type]][match(query, mapped[[from_type]])]
    }
  }
  n_orgdb <- sum(!is.na(result))

  # 2) Fill the gaps from UniProt (UNIPROT <-> SYMBOL only)
  if (any(is.na(result)) && all(c(from_type, to_type) %in% c("UNIPROT", "SYMBOL"))){
    gm <- get_uniprot_genemap(species, ...)
    gm[, gene := sub(";.*", "", gene)]
    fill <- if (from_type == "UNIPROT") gm$gene[match(query, gm$uniprot)] else gm$uniprot[match(query, gm$gene)]
    result[is.na(result)] <- fill[is.na(result)]
  }

  message(sprintf("Mapped %.1f%% of input IDs from %s to %s (%d / %d; %d from OrgDb, %d from UniProt).",
                  100 * mean(!is.na(result)), from_type, to_type, sum(!is.na(result)), length(result),
                  n_orgdb, sum(!is.na(result)) - n_orgdb))
  names(result) <- names(data_in)
  result
}


# =============================================================================
# 5. Pathway enrichment (D2)
# =============================================================================

do_enrichmap = function(M,cutHeight,topN,flag_cluster_col=TRUE,flag_p=FALSE,
                        path_to_databasefiles='',database=NULL,
                        databasefile = NULL, database_in=NULL,
                        column_order_values = NULL,minGSSize = 10, allGenes = NULL){

# Inputs explained:
# M table with two columns: "genes" and "group".
# cutHeight = 0.8
# topN = 3
# If flag_cluster_col, then cluster the columns ("groups") as well as the rows. 
# If flag_p, then use p-value rather than adjusted p-value.
# database_in should be in a dataframe format with two columns: "term" and "gene"
#   (get_msigdb() in databases.R returns exactly this, so no GMT files are needed)
  

# Get enrichment database
# if (database=='GOBP'){
# gmt <- clusterProfiler::read.gmt (paste(path_to_databasefiles,"c5.bp.v7.1.symbols.gmt",sep=""))
# setDT(gmt)
# }
# if (database=='Pathways'){
# gmt <- clusterProfiler::read.gmt (paste(path_to_databasefiles,"c2.cp.v7.5.1.symbols.gmt",sep=""))
# setDT(gmt)
# }
# if (!is.null(databasefile)){  
# gmt <- clusterProfiler::read.gmt (paste0(path_to_databasefiles,databasefile))
# setDT(gmt)
# }
if (!is.null(database_in)){
  gmt = database_in
} else {
  gmt <- clusterProfiler::read.gmt(paste0(path_to_databasefiles,databasefile))
  setDT(gmt)
}
  
# Background (universe): by default all genes in the gene-set database.
# Pass allGenes (e.g. all genes detected in your experiment) to use a custom background.
# (Previously the condition was inverted, so a supplied allGenes was silently overwritten.)
if (is.null(allGenes)){
  allGenes <- unique(gmt$gene)
}

# Prepare input data frame
M$genes=as.character(M$genes)
setDT(M)

# Number of input genes per comparison, shown as a bar plot above the heatmap
geneCounts = M[, .(n = uniqueN(genes)), by = group]
geneCounts = setNames(geneCounts$n, geneCounts$group)

# Perform Enrichment
groups = unique(M$group)
# forked parallelism is not available on Windows, so run serially there
numProcessors = if (.Platform$OS.type == "windows") 1 else max(1, floor(parallel::detectCores()/2))
message ("Computing enrichments on ", length(groups), " groups")
enrichList =  pbapply::pblapply(groups, 
                                 function(g){
                                   setDT(as.data.frame(clusterProfiler::enricher(unique(M[group==g,]$genes),
                                                                                 pAdjustMethod="fdr",
                                                                                 universe=allGenes,
                                                                                 TERM2GENE = gmt,
                                                                                 qvalueCutoff = 1.1,
                                                                                 pvalueCutoff = 1.1,
                                                                                 minGSSize = minGSSize)))
                                 },
                                 cl=numProcessors)
names(enrichList) <- groups

#Combine all results into a single table
enrichTable = rbindlist(enrichList, idcol="group")

# Simplification of overlapping GO terms
if (flag_p){ # If flag_p, then use p-value rather than adjusted p-value.
  enrichTable$p.adjust=enrichTable$pvalue
}

# Simplification and the heatmap need at least 2 enriched terms to cluster.
# With fewer, return the enrichment table only (no heatmap) instead of stopping with an hclust error.
sig_ids <- if ("p.adjust" %in% names(enrichTable)) unique(enrichTable$ID[which(enrichTable$p.adjust < 0.05)]) else character(0)
if (length(sig_ids) < 2){
  message("Only ", length(sig_ids), " enriched term(s) (adj.p < 0.05); returning the table without simplification or heatmap")
  return(list(enrichTable, NULL, NULL, NULL))
}

simp <- simplifyEnrichBySimilarUniverseMembership(enrichTable, 
                                                  gmt, 
                                                  groupColumn="group",
                                                  cutHeight = cutHeight, 
                                                  broadest=FALSE)
simpTable = setDT(simp[[1]])

# Draw Heatmaps
#if (flag_plot){
a <- enrichHeatmapBestPerGroup(simp[[1]], simp[[2]],
                               groupColumn = "group",
                               topN = topN,
                               cluster_columns=flag_cluster_col, 
                               reduceRedundantsAcrossGroup=FALSE,
                               column_order_values = column_order_values,
                               geneCounts = geneCounts) # passed on to heatmapNumbered() for the bar plot
simpTable_topN = a[[1]]
#}

return(list(enrichTable,a,simpTable,simpTable_topN))

}


simplifyEnrichBySimilarUniverseMembership <- function(enrichResultsTable, gmt, groupColumn="bait", cutHeight = 0.99, broadest=TRUE){
  if (length(unique(enrichResultsTable$ID)) < 2){
    message ("Nothing to simplify")
    return (list (enrichResultsTable, data.frame()))
  }
  setDT(enrichResultsTable); setDT(gmt)
  
  ##Select Significant Terms
  target_overrep_sig <- enrichResultsTable[p.adjust < 0.05,]#qvalue < 0.01]
  #target_overrep_sig <- enrichResultsTable[pvalue < 0.05,]#qvalue < 0.01]
  
  ##Prepare Significant GO Term Jaccard Similarity Matrix
  sig_go_terms <- unique(target_overrep_sig$ID)
  
  message ("Computing universal gene overlap between ", length(sig_go_terms), " significant GO terms from ", length(unique(enrichResultsTable[[groupColumn]])), " ", groupColumn, "(s)")
  
  gmt.subset <- gmt[term %in% sig_go_terms, .(term=factor(term), gene=factor(gene))]
  termByGeneMat <-  Matrix::sparseMatrix(as.integer(gmt.subset$term), as.integer(gmt.subset$gene), 
                                         dimnames=list(levels(gmt.subset$term), 
                                                       levels(gmt.subset$gene)))
  
  go_dist_mat <- dist(termByGeneMat, method="binary")
  hc <- hclust(go_dist_mat) # you can change the method here.
  
  clusters <- cutree(hc, h=cutHeight)
  clusters <- data.table (cluster = as.numeric(clusters), ID = attributes(clusters)$names )
  
  message ("GO terms clustered into ", max(clusters$cluster), " clusters")
  
  ## go gene set lengths
  gmt.setLengths <- gmt[,.(setSize = length(unique(gene))), by = term]
  clusters <- merge (clusters, gmt.setLengths, by.x="ID", by.y="term")
  
  ## local gene set lengths  -- this will only count those genes that make it into an enriched list
  id.gene.long <- enrichResultsTable[,.(ID, gene = unlist(strsplit(geneID, "/"))), by = seq_len(nrow(enrichResultsTable))]
  genesPerTerm <- id.gene.long[,.(localSetLength=length(unique(gene))), by = ID]
  clusters <- merge (clusters, genesPerTerm, by.x="ID", by.y="ID")
  #setorder(clusters, -localSetLength)  # the tie breaker
  
  clusterInfo <- merge (enrichResultsTable, clusters, by = "ID")
  
  clusterInfo[,maxSet := max (setSize), by = c("cluster", groupColumn)]
  #winners <- clusterInfo[,.SD[which(count == maxSet)], by = .(cluster, Bait)]  #keeps all ties...needs updating
  
  if (broadest){
    winners <- clusterInfo[p.adjust < 0.05,.SD[which.max(setSize),],by=c("cluster", groupColumn)]  #chooses the first in case of tie breakers
    message (length(unique(winners$ID)), " representative GO terms choosing the BROADEST significant term per GO-cluster per ", groupColumn)
  }else{
    winners <- clusterInfo[p.adjust < 0.05,.SD[which.min(p.adjust),],by=c("cluster", groupColumn)]  #chooses the first in case of tie breakers
    message (length(unique(winners$ID)), " representative GO terms choosing the MOST significant term per GO-cluster per ", groupColumn)
  }
  result <- enrichResultsTable[ID %in% winners$ID,]
  list(simplified = result, clusterInfo = clusterInfo)
}


fixMsigdbGONames <- function(names){
  names <- gsub("^GO_","",names)
  names <- gsub("_"," ",names)
  names <- tolower(names)
  patterns <- c("\\bdna|dna\\b", 
                "\\brna|rna\\b", 
                "\\bgtp|gtp\\b",
                "\\batp|atp\\b"
  )
  replacements <- c("DNA", 
                    "RNA", 
                    "GTP",
                    "ATP"
  ) 
  for (i in seq_along(patterns)){
    names <- gsub (patterns[i], replacements[i], names)
  }
  #capitalize first letter
  substr(names,1,1) <- toupper(substr(names,1,1))
  return(names)
}

fixMsigdbPNames <- function(names){
  names <- gsub("_"," ",names)
  names = tolower(names)
  return(names)
}


heatmapNumbered <- function (main.mat, counts.mat, negCols = NULL, title="",
                             borderMatrix = NULL, borderColFun = NULL,
                             borderMM = 2,
                             brewerPalette = "Blues",
                             show_column_dend = TRUE,
                             border = TRUE, column_order_values = column_order_values,
                             geneCounts = NULL, ...){
  # geneCounts: optional named vector (names = heatmap columns) of input genes per comparison,
  #             drawn as a bar plot above the heatmap
  Blues = colorRampPalette(RColorBrewer::brewer.pal(9, brewerPalette))
  colors <- Blues(100)
  #colors = circlize::colorRamp2 (breaks=seq(from=0, to = max(main.mat), length.out=101), colors =colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(101))
  #colors = circlize::colorRamp2(c(0, max(main.mat)), c("#EEEEEE", "blue"))
  
  heatmap_legend_param = list(legend_direction="horizontal", title = "global enrich\n-log10(adj.p)")
  
  if (!is.null(negCols)){
    colors <- circlize::colorRamp2 (breaks=seq(from=0, to = max(main.mat), length.out=101), colors =colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(101))
    main.mat[,negCols] = -main.mat[, negCols]
    heatmap_legend_param = c (heatmap_legend_param, list(at=c(-4,-2,0,2,4), labels = c(4,2,0,2,4)) )
  }
  
  
  ## Bar plot of input gene counts per comparison, above the heatmap
  top_ann = NULL
  if (!is.null(geneCounts)){
    n_genes = as.numeric(geneCounts[colnames(main.mat)])   # match counts to heatmap columns by name
    n_genes[is.na(n_genes)] = 0
    top_ann = HeatmapAnnotation(`Input genes` = anno_barplot(n_genes, add_numbers = TRUE, numbers_rot = 0,
                                                             gp = gpar(fill = "grey70", col = NA),
                                                             height = unit(2, "cm")),
                                annotation_name_side = "left")
  }

  ##Plot main figure heatmap
  hm <- ComplexHeatmap::Heatmap(as.matrix(main.mat), col = colors, border = border, rect_gp = gpar(col = "grey", lwd = 1), column_title = title,
                                top_annotation = top_ann,
                                column_names_rot = 90, row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10),
                                show_row_dend = TRUE, show_column_dend = show_column_dend, heatmap_legend_param = heatmap_legend_param,
                                row_names_max_width = max_text_width(rownames(main.mat), gp = gpar(fontsize = 12)),
                                column_order = column_order_values,
                                cell_fun = function(j, i, x, y, width, height, fill) {
                                  if (!is.null(borderMatrix) & !is.null(borderColFun)){
                                    lwd <- unit(borderMM,"mm")
                                    grid.rect(x, y, width, height, gp = gpar(fill = borderColFun(borderMatrix[i,j]), col = NA))
                                    grid.rect(x, y, width-lwd, height-lwd, gp = gpar(fill = fill, col = NA))
                                  }
                                  if (!is.na(counts.mat[i,j])){
                                    color <- ifelse(main.mat[i,j] < (-log10(0.05)), "darkgrey", "red")
                                    grid.text(sprintf("%s",counts.mat[i, j]), x, y, gp = gpar(fontsize=10, col=color))
                                  }
                                }, ...) #+1  # this makes it a list!
  
  
  if (!is.null(borderMatrix) & !is.null(borderColFun)){
    legendList <-  list (Legend(col_fun = colorFun, title= "viral enrich\n-log10(p)"))
    hm <- hm + 1 # this makes it a list so I can add annotations to it
  } else{
    legendList <- list()
  }
  
  draw(hm,heatmap_legend_side="top", annotation_legend_list = legendList)
  
  invisible (hm)
  
}


enrichHeatmapBestPerGroup <- function(simplifiedEnrichTable, fullEnrichTable, groupColumn="bait", topN = 1, title="", cols = NULL, 
                                      negCols = NULL, reduceRedundantsAcrossGroups=TRUE, column_order_values=column_order_values,...){
  setorder(simplifiedEnrichTable, p.adjust)
  bestTermPerBait <- simplifiedEnrichTable[p.adjust<0.05,.(ID=ID[1:topN]),by=groupColumn]
  if(reduceRedundantsAcrossGroups){  
    #reduce redundant based on clusters in fullEnrichTable
    countsByID <- fullEnrichTable[ID %in% bestTermPerBait$ID, .(geneCount  = length(unique(unlist(strsplit(geneID, "/"))))), by = .(ID, cluster)]
    # get the term with most genes across whole dataset per term-cluster
    setorder(countsByID, -geneCount)
    bestTerms <- countsByID[,.SD[1],by=cluster]$ID
  } else bestTerms <- unique(bestTermPerBait$ID)
  
  main.wide <- dcast (fullEnrichTable[ID %in% bestTerms], as.formula(paste("Description", groupColumn, sep="~")), value.var="p.adjust")
  for(col in unique(c(cols,negCols))){
    if (is.null(main.wide[[col]])) main.wide[[col]] <- NA
  }
  main.mat <- -log10(as.matrix(main.wide, rownames = "Description"))
  
  main.mat[is.na(main.mat)] <- 0
  #main.mat[main.mat > 5] <- 5
  if (all(grepl("^GO_", rownames(main.mat)))){
    rownames(main.mat) <- fixMsigdbGONames(rownames(main.mat))}
  else {
    rownames(main.mat) <- fixMsigdbPNames(rownames(main.mat))
  }
  
  counts.wide <- dcast (fullEnrichTable[ID %in% bestTerms], as.formula(paste("Description", groupColumn, sep="~")), value.var="Count")
  for(col in unique(c(cols,negCols))){
    if (is.null(counts.wide[[col]])) counts.wide[[col]] <- NA
  }
  counts.mat <- as.matrix(counts.wide, rownames = "Description")
  
  
  geneTable <- fullEnrichTable[ID %in% bestTerms, .(gene = unlist(strsplit(geneID, split="/"))),by = ID]
  geneTable[,cleanName := fixMsigdbGONames(ID)]
  
  if (!is.null(cols)){
    if (! all(cols %in% colnames(counts.mat) & cols %in% colnames(main.mat))){
      message ("Not all requested columns for heatmap found in data")
      message ("main.mat ", paste(colnames(main.mat), collapse=" "))
      message ("counts.mat ", paste(colnames(counts.mat), collapse=" "))
    }else{
      counts.mat<- counts.mat[,cols]
      main.mat<- main.mat[,cols]
    }
    
  }
  
  hm <- heatmapNumbered (main.mat, counts.mat, negCols, title, column_order_values=column_order_values, ...)
  
  invisible(list(geneTable = geneTable, main.mat = main.mat, counts.mat = counts.mat, hmList = hm))
}


# =============================================================================
# 6. Kinase activity from OmniPath (D3, phospho only)
# =============================================================================

# This script calculates kinase activities based on phosphoproteomics data.
# Required input:

# 1) log2FCData should have 4 columns: group, log2FC, pos, Protein. 

# group has the comparison being considered.
# log2FC has the log2FC per phosphorylation site.
# pos has the amino acid and position of the detected phosphorylation site
# Protein has the uniprot protein ID. 

# 2) kinaseDataSource should be either "Omnipath" or NULL.

# 3) refreshOmnipath (default TRUE): download the latest network from the OmniPath web
#    service (omnipathdb.org). FALSE = reuse the newest saved snapshot if there is one.

# 4) out_dir: folder for the output files (default: the results folder set in the Rmd).

# Output files (in out_dir):
#   kinaseSummaryScores.txt    - kinase activity scores per comparison
#   kinaseSubstrateData.txt    - your sites that mapped to the network (subset)
#   omnipath_kinaseSubstrateNetwork_<date>.txt - the FULL OmniPath network used,
#                                saved as a snapshot for reproducibility


do_KinaseActivity_Network = function(log2FCData, kinaseDataSource="Omnipath", species="HUMAN",
                                     refreshOmnipath = TRUE,
                                     out_dir = getOption("ptm.results_dir", ".")){
  
  # Load the kinase-substrate network ONCE (instead of once per comparison)
  # (if kinaseDataSource is NULL, kinaseActivity() falls back to kinaseDataFile as before)
  kinaseData <- NULL
  if (!is.null(kinaseDataSource)){
    # Snapshots are looked for in this project's results and in the shared database cache
    cache_dir <- if (exists("db_dir")) db_dir() else NULL
    kinaseData <- loadKinaseDataOmniPath(species = species, refresh = refreshOmnipath,
                                         snapshot_dirs = c(out_dir, cache_dir))
    
    if (!is.null(attr(kinaseData, "snapshot"))){
      # A saved snapshot is being used: record that, don't save it again as a new one
      if (exists("db_log")) db_log("OmniPath enzyme-substrate (saved snapshot)",
                                   "https://omnipathdb.org/enzsub", attr(kinaseData, "snapshot"))
    } else {
      # Save the full network snapshot used for this run (not just the mapped subset),
      # in the project results and in the shared cache (used as a fallback if a later download fails)
      dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
      networkFile <- file.path(out_dir, paste0("omnipath_kinaseSubstrateNetwork_", Sys.Date(), ".txt"))
      fwrite(kinaseData, file = networkFile, sep = '\t', quote = FALSE)
      if (!is.null(cache_dir)) fwrite(kinaseData, file = file.path(cache_dir, basename(networkFile)), sep = '\t', quote = FALSE)
      message("Saved OmniPath network snapshot to ", networkFile)
      if (exists("db_log")) db_log("OmniPath enzyme-substrate", "https://omnipathdb.org/enzsub", networkFile)
    }
  }
  
  # Run Kinase Activity for every label separately
  labels <- unique(log2FCData$group)
  kinActList <- lapply (labels, FUN=function(lab){
    kinaseActivity(log2FCData[group==lab,], kinaseDataSource=kinaseDataSource,
                   kinaseData=kinaseData, uniprot=TRUE)})
  names(kinActList) <- labels
  
  kinActFull.scores <- rbindlist(lapply(kinActList, FUN = function(x)x$scores), idcol="Label")
  kinActFull.mapped <- rbindlist(lapply(kinActList, FUN = function(x)x$kinaseMapped)) # Label is already in these tables
  
  # Save out data
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  write.table(kinActFull.scores, file=file.path(out_dir, "kinaseSummaryScores.txt"), row.names=F, quote=F, sep='\t')
  write.table(kinActFull.mapped, file=file.path(out_dir, "kinaseSubstrateData.txt"), row.names=F, quote=F, sep='\t')
  invisible(kinActFull.scores)
}

#-----------------------------------------------------------------------------------------------

kinaseActivity <- function (log2FCData, kinaseDataSource="Omnipath", kinaseData = NULL,
                            kinaseDataFile = bouquet_file("kinaseSiteList_BachmanGyoriSorger2019.csv.gz"),
                            plots=FALSE, outlierLog2FC = 10, requireAAMatch = TRUE, uniprot=FALSE, do.sea = FALSE,
                            sea.center = FALSE, sea.weightParam = 1){
  
  #kinaseDataSource should be "Omnipath" or NULL. If NULL, then kinaseDataFile will be used.
  
  # deal with columns
  requiredColumns <- c("log2FC", "pos", "Protein")
  if (requireAAMatch) requiredColumns <- c(requiredColumns, "aa")
  if (uniprot){
    requiredColumns <- c(requiredColumns, "uniprot")
  }else{
    requiredColumns <- c(requiredColumns, "Gene_Name")
  }
  
  missingColumns <- setdiff(requiredColumns, colnames(log2FCData))
  if (length(missingColumns) > 0) stop (length(missingColumns), " missing column(s) ", paste0(missingColumns, collapse=", "), " in log2FCData  table")
  
  if (uniprot){
    mergeXCols <- c("uniprot")
    mergeYCols <- c("TARGET_UP_ID")
  }else{
    mergeXCols <- c("Gene_Name")
    mergeYCols <- c("TARGET_GENE_NAME")
  }
  if(requireAAMatch){
    mergeXCols <- c(mergeXCols, "aa")
    mergeYCols <- c(mergeYCols, "TARGET_RES")
  }
  mergeXCols <- c(mergeXCols, "pos")
  mergeYCols <- c(mergeYCols, "TARGET_POS")
  
  # done with columns
  
  
  # # Select complete rows only, if selected above
  # if (flag_complete.across.contrasts){
  #   D.wide = dcast(D,Protein~Label,value.var="log2FC")
  #   D.wide = as.matrix(D.wide[,2:ncol(D.wide)],rownames=D.wide$Protein)
  #   D.wide[is.infinite(D.wide)]=NA
  #   inds = complete.cases(D.wide)
  #   message(paste0(round(sum(inds)/length(inds)*100,2),"% of phosphorylation sites were kept after removing incomplete cases. Make sure you selected the right contrasts2include as that will affect what is kept/removed."))
  #   D.complete.cases = D[as.logical(match(D$Protein,rownames(D.wide)[inds],nomatch=F)),]
  #   D = D.complete.cases
  # }
  
  
  if ("representative" %in% colnames(log2FCData) & any(log2FCData$representative == FALSE)){
    warning ("data.table includes a 'representative' column, but some are set to FALSE. Did you forget to pre-filter?")
  }
  
  # Only load the network here if it was not passed in already
  if (is.null(kinaseData)){
    if (is.null(kinaseDataSource)){
      kinaseData <- loadKinaseData (kinaseDataFile)
    } else {
      kinaseData <- loadKinaseDataOmniPath()
    }
  }
  
  # clean up data
  ## remove NA and infinite
  cleanData <- log2FCData[!is.na(log2FC) & is.finite(log2FC)]
  if (plots) boxplot(cleanData$log2FC)
  ## outliers
  if (!is.na(outlierLog2FC)){
    cleanData <- cleanData[abs(log2FC) < abs(outlierLog2FC)]
    if (plots) boxplot(cleanData$log2FC)
  }
  
  
  
  kinaseMapped <- merge (cleanData, kinaseData, by.x=mergeXCols, by.y=mergeYCols)
  
  if (uniprot) kinaseMapped[, Gene_Name := uniprot]
  if (! "aa" %in% colnames(kinaseMapped)) kinaseMapped[,aa := ""]
  
  if(plots){
    plot (density(cleanData[,log2FC]))
    lines (density(kinaseMapped[,log2FC]), col="red")
    #legend ...
  }
  
  # the denominator is SD of the full dataset, whether it maps to a kinase or not:
  fullSD <- sd(cleanData[,log2FC])
  #each kinase mean will be compared to global mean
  fullMean <- mean(cleanData[,log2FC])
  
  # Columns holding each substrate's UniProt ID and gene symbol (the one not used for matching comes from the network)
  upCol   <- if (uniprot) "uniprot" else "TARGET_UP_ID"
  geneCol <- if (uniprot) "TARGET_GENE_NAME" else "Gene_Name"

  # summarize per kinase, Z score based on standard error of the mean log2FC per kinase
  # substrate lists show each site with its log2FC, e.g. "P31749_S473(1.52)" and "AKT1_S473(1.52)"
  scores <- kinaseMapped[,.(Z = (unique(.SD)[,mean(log2FC)] - fullMean)*sqrt(length(unique(Protein)))/fullSD,
                            #N_unique_proteins = length(unique(Protein)),#.N,  # don't double count when two sites are in the same peptide
                            N = length(unique(paste(Gene_Name, paste0(aa, pos), sep="_"))), # Number of unique sites
                            meanLog2FC  = mean(log2FC),
                            uniprot_substrates  = paste(unique(paste0(get(upCol),   "_", aa, pos, "(", round(log2FC, 2), ")")), collapse=","),  # 'unique' because some sites will occur multiple times in different phospho-combos
                            genename_substrates = paste(unique(paste0(get(geneCol), "_", aa, pos, "(", round(log2FC, 2), ")")), collapse=",")),
                         by = CTRL_GENE_NAME,
                         .SDcols = c("Protein", "log2FC")]
  
  # p values.  2*pnorm... makes it two-tailed. This is different from KSEApp, but I think it is more appropriate here where we include both positive and negative effects
  scores[,pValue := 2*pnorm(-abs(Z), lower.tail= TRUE)]
  scores[, sigScore := -log10(pValue) * ifelse(Z < 0, -1, 1) ]
  scores[, fdr.BH := p.adjust(pValue, method = "fdr")]  #method=fdr is alias for "BH" Benjamini & Hochberg
  scores[,c("bgMean", "bgSD") := .(fullMean, fullSD)]
  
  
  if (do.sea == TRUE){
    seaScores <- kinaseActivity.sea(kinaseMapped, fullData = cleanData, center = sea.center, gseaWeightParam = sea.weightParam, nperm = 10000)
    setnames(seaScores, old = c("pval", "padj"), new = paste0(c("pval", "padj"), ".sea"))
    scores <- merge (scores, seaScores, by.x = c("CTRL_GENE_NAME"), by.y = c("pathway"),suffixes = c("", ".sea"), all = TRUE)
  }
  
  setorder(scores, pValue)
  
  if (plots){
    hist (scores$pValue, breaks=100)
    p <- ggplot(scores, aes(x=Z, y = -log10(pValue), col = sqrt(N))) + geom_point()
    print (p)
    p <- ggplot(scores, aes(x=Z, y = -log10(fdr.BH), col = sqrt(N))) + geom_point()
    print (p)
    p <- ggplot(scores, aes(x=meanLog2FC, y = -log10(pValue), col = sqrt(N))) + geom_point()
    print (p)
    p <- ggplot(scores, aes(x=meanLog2FC, y = -log10(fdr.BH), col = sqrt(N))) + geom_point() #+ scale_color_steps(low="white", high="navy", breaks=c(0,1,3,10,50,100,200)) + theme_classic()
    print (p)
    
    p <- BarplotKinaseActivities(scores, kinaseMapped, max_pValue = 0.05)
    print (p)
    p <- BarplotKinaseActivities(scores, kinaseMapped, max_pValue = 1.0)
    print (p)
    
  }
  return (list(scores = scores, kinaseMapped = kinaseMapped))
}

#-----------------------------------------------------------------------------------------------

# Downloads the enzyme-substrate network directly from the OmniPath web service
# (https://omnipathdb.org/enzsub, the same database browsed at https://explore.omnipathdb.org/),
# without the OmnipathR package.
#
# refresh = TRUE:  download a fresh copy
# refresh = FALSE: use the newest saved snapshot if there is one, otherwise download
# snapshot_dirs:   folders with saved networks (omnipath_kinaseSubstrateNetwork_<date>.txt).
#                  If the download fails (e.g. the server is down), the newest snapshot is used, with a warning.
loadKinaseDataOmniPath <- function(species = "HUMAN", removeNonKinases = TRUE, refresh = FALSE,
                                   snapshot_dirs = NULL, ...){
  species <- toupper(species)
  organismID <- unname(c(HUMAN = 9606L, RAT = 10116L, MOUSE = 10090L)[species])
  if(is.na(organismID))
    stop(species, " not a recognized organism")
  
  # Newest saved snapshot (if any)
  snapshot_dirs <- as.character(snapshot_dirs)                    # NULL (no folders) becomes an empty list of folders
  snapshot_dirs <- snapshot_dirs[nzchar(snapshot_dirs)]
  snaps <- unlist(lapply(snapshot_dirs[dir.exists(snapshot_dirs)], list.files,
                         pattern = "^omnipath_kinaseSubstrateNetwork_.*\\.txt$", full.names = TRUE))
  snap <- if (length(snaps)) snaps[which.max(file.mtime(snaps))] else NULL
  read_snapshot <- function(){
    kinaseData <- fread(snap)
    attr(kinaseData, "snapshot") <- snap
    kinaseData
  }
  if (!refresh && !is.null(snap)){
    message("Using saved OmniPath network (refresh = FALSE): ", snap)
    return(read_snapshot())
  }
  
  # Download from the OmniPath web service
  url <- paste0("https://omnipathdb.org/enzsub?genesymbols=yes&fields=sources&format=tsv&organisms=", organismID)
  enzsub <- tryCatch({
    tmp <- tempfile(fileext = ".tsv")
    utils::download.file(url, tmp, quiet = TRUE, mode = "wb")
    d <- fread(tmp)
    if (nrow(d) == 0) stop("the server returned no data")
    d
  }, error = function(e) e)
  
  if (inherits(enzsub, "error")){
    if (is.null(snap))
      stop("Could not download the OmniPath network from ", url, ", and no saved snapshot was found.\n",
           "This is usually a temporary server problem; try again later.\n",
           "Original error: ", conditionMessage(enzsub), call. = FALSE)
    warning("Could not download the OmniPath network (", conditionMessage(enzsub), ").\n",
            "Using the most recent saved snapshot instead: ", snap, call. = FALSE)
    return(read_snapshot())
  }
  setDT (enzsub)
  
  # Record what was used, for reproducibility
  message("OmniPath enzyme-substrate data: ", nrow(enzsub), " interactions from omnipathdb.org; retrieved ", Sys.Date())
  
  neededCols <- c("enzyme_genesymbol", "substrate_genesymbol", "residue_type",
                  "residue_offset", "substrate", "modification")
  missingCols <- setdiff(neededCols, colnames(enzsub))
  if (length(missingCols) > 0)
    stop("OmniPath output is missing expected column(s): ", paste(missingCols, collapse=", "),
         ". The OmniPath web service output format may have changed.")
  # for consistency with loadKinaseData I need columns: "CTRL_GENE_NAME", "TARGET_GENE_NAME", "TARGET_RES", "TARGET_POS", "TARGET_UP_ID"
  
  if (any(grepl("ProtMapper", enzsub$sources))){
    message ("OmniPath includes ProtMapper data which has many non-kinases.")
    
    if(removeNonKinases == TRUE) {
      if (!species %in% c("HUMAN", "MOUSE")) stop ("I don't yet know how to find kinases for ", species)
      orgdb <- suppressMessages(get_orgdb(tolower(species)))
      if (!is.null(orgdb)){
        # all genes annotated with GO:0016301 (kinase activity) or any of its child terms
        message ("Removing all enzymes that are not kinases according to GO:0016301 (", orgdb$packageName, ")")
        kinaseTable <- suppressMessages(AnnotationDbi::select(orgdb, keys = "GO:0016301", keytype = "GOALL", columns = "SYMBOL"))
        kinases <- unique(kinaseTable$SYMBOL)
      } else {
        # built-in human kinase list (mouse gene symbols matched case-insensitively)
        message ("Removing all enzymes that are not kinases according to BouQuet's built-in kinase list")
        human_kinases <- fread(bouquet_file("human_kinases.txt"))$Kinase
        kinases <- unique(enzsub$enzyme_genesymbol[toupper(enzsub$enzyme_genesymbol) %in% human_kinases])
      }
      #label and filter into a new table
      enzsub[, is_kinase := enzyme_genesymbol %in% kinases]
      enzsub <- enzsub[is_kinase == TRUE]
    }
  }
  
  setnames(enzsub, 
           old = c("enzyme_genesymbol", "substrate_genesymbol", "residue_type", "residue_offset", "substrate"), 
           new = c("CTRL_GENE_NAME", "TARGET_GENE_NAME", "TARGET_RES", "TARGET_POS", "TARGET_UP_ID"))
  
  return (enzsub[modification == "phosphorylation"])
}

#-----------------------------------------------------------------------------------------------

# Used when kinaseDataSource = NULL (static Bachman, Gyori & Sorger 2019 kinase-substrate list,
# shipped with BouQuet because it is small and no longer updated).
# (Was previously only defined in do_kinase_enrichment.R, so kinaseDataSource = NULL failed.)
loadKinaseData <- function (kinaseDataFile = bouquet_file("kinaseSiteList_BachmanGyoriSorger2019.csv.gz")){
  kinaseSubstrates <- fread (kinaseDataFile)
  kinaseSubstrates <- kinaseSubstrates[CTRL_IS_KINASE == TRUE]
  # for CTRL_NS == FPLX the CTRL_ID looks like a kinase name but CTRL_GENE_NAME is blank
  kinaseSubstrates[CTRL_NS == "FPLX" & CTRL_GENE_NAME == "", CTRL_GENE_NAME := CTRL_ID]
  return (kinaseSubstrates)
}


# =============================================================================
# 7. Built-in tables (shipped with BouQuet in inst/extdata)
# =============================================================================

# Full path to a table that comes with BouQuet, e.g. bouquet_file("human_kinases.txt")
# Built-in tables: human_kinases.txt, kinaseSiteList_BachmanGyoriSorger2019.csv.gz,
#                  Regulatory_sites_<date>.tsv.gz (PhosphoSitePlus; see inst/extdata/README.md)
bouquet_file <- function(name){
  path <- system.file("extdata", name, package = "BouQuet")
  if (path == "") stop("'", name, "' is not one of BouQuet's built-in tables. Available: ",
                       paste(list.files(system.file("extdata", package = "BouQuet")), collapse = ", "), call. = FALSE)
  path
}
