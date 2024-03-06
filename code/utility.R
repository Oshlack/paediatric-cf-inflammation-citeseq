
add_gene_information <- function(sce){
  require(AnnotationDbi)
  require(EnsDb.Hsapiens.v86)
  require(Homo.sapiens)
  require(msigdbr)
  
  # Add chromosome location so we can filter on mitochondrial genes.
  location <- mapIds(
    x = EnsDb.Hsapiens.v86,
    # NOTE: Need to remove gene version number prior to lookup.
    keys = rownames(sce),
    keytype = "GENEID",
    column = "SEQNAME")
  rowData(sce)$CHR <- location
  # Additional gene metadata from ENSEMBL and NCBI
  # NOTE: These columns were customised for this project.
  ensdb_columns <- c(
    "GENEBIOTYPE", "GENENAME", "GENESEQSTART", "GENESEQEND", "SEQNAME", "SYMBOL")
  names(ensdb_columns) <- paste0("ENSEMBL.", ensdb_columns)
  stopifnot(all(ensdb_columns %in% columns(EnsDb.Hsapiens.v86)))
  ensdb_df <- DataFrame(
    lapply(ensdb_columns, function(column) {
      mapIds(
        x = EnsDb.Hsapiens.v86,
        keys = rownames(sce),
        keytype = "GENEID",
        column = column,
        multiVals = "CharacterList")
    }),
    row.names = rownames(sce))
  # NOTE: Can't look up GENEID column with GENEID key, so have to add manually.
  ensdb_df$ENSEMBL.GENEID <- rownames(sce)
  # NOTE: Homo.sapiens combines org.Hs.eg.db and
  #       TxDb.Hsapiens.UCSC.hg19.knownGene (as well as others) and therefore
  #       uses entrez gene and RefSeq based data.
  
  # NOTE: These columns were customised for this project.
  ncbi_columns <- c(
    # From TxDB: None required
    # From OrgDB
    "ALIAS", "ENTREZID", "GENENAME", "REFSEQ", "SYMBOL")
  names(ncbi_columns) <- paste0("NCBI.", ncbi_columns)
  stopifnot(all(ncbi_columns %in% columns(Homo.sapiens)))
  ncbi_df <- DataFrame(
    lapply(ncbi_columns, function(column) {
      mapIds(
        x = Homo.sapiens,
        keys = rownames(sce),
        keytype = "ENSEMBL",
        column = column,
        multiVals = "CharacterList")
    }),
    row.names = rownames(sce))
  rowData(sce) <- cbind(rowData(sce), ensdb_df, ncbi_df)
  
  # Some useful gene sets
  ribo_set <- grep("^RP(S|L)", rownames(sce), value = TRUE)
  # NOTE: A more curated approach for identifying ribosomal protein genes
  # http://bioconductor.org/books/3.14/OSCA.advanced/more-hvgs.html#feature-selection-positive
  c2_sets <- msigdbr(species = "Homo sapiens", category = "C2")
  ribo_set <- union(
    ribo_set,
    c2_sets[c2_sets$gs_name == "KEGG_RIBOSOME", ]$human_ensembl_gene)
  rowData(sce)$is_ribo <- rownames(sce) %in% ribo_set
  rowData(sce)$is_sex <- any(rowData(sce)$ENSEMBL.SEQNAME %in% c("X", "Y"))
  rowData(sce)$is_mito <- rowData(sce)$CHR == "MT"
  rowData(sce)$is_pseudogene <- any(grepl("pseudogene", rowData(sce)$ENSEMBL.GENEBIOTYPE))
  
  sce
}
