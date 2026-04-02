#' Add gene-level annotation and feature flags to a SingleCellExperiment object
#'
#' Enriches the rowData of a SingleCellExperiment with gene metadata from
#' Ensembl (via EnsDb.Hsapiens.v86) and NCBI (via Homo.sapiens), and appends
#' logical flag columns indicating mitochondrial, ribosomal, sex-chromosome,
#' and pseudogene status.
#'
#' @param sce A SingleCellExperiment object whose rownames are Ensembl gene IDs
#'   (GENEID format, e.g. "ENSG00000000003").
#'
#' @return The input \code{sce} with additional columns appended to
#'   \code{rowData}: Ensembl and NCBI annotation fields, plus logical flags
#'   \code{is_ribo}, \code{is_sex}, \code{is_mito}, and \code{is_pseudogene}.
#'
#' @details
#' Annotation columns added (prefixed \code{ENSEMBL.}) come from
#' EnsDb.Hsapiens.v86 and include biotype, gene name, genomic coordinates,
#' chromosome, and symbol. NCBI columns (prefixed \code{NCBI.}) are drawn from
#' Homo.sapiens (org.Hs.eg.db + TxDb.Hsapiens.UCSC.hg19.knownGene) and
#' include aliases, Entrez ID, RefSeq accessions, gene name, and symbol.
#' Ribosomal gene membership combines a name-pattern grep with the curated
#' MSigDB C2 KEGG_RIBOSOME gene set.
#'
#' @examples
#' \dontrun{
#'   sce <- add_gene_information(sce)
#' }
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

#' Parse a GMT file into a named list of gene sets
#'
#' Reads a Gene Matrix Transposed (GMT) file and returns a named list where
#' each element corresponds to one gene set. The gene set name (column 1) is
#' used as the list element name, and the remaining columns (after the URL in
#' column 2) contain the Entrez gene IDs for that set.
#'
#' @param file_path Character string; path to a .gmt file.
#'
#' @return A named list of character vectors. Each name is a gene set
#'   identifier and each vector contains the Entrez IDs belonging to that set.
#'
#' @examples
#' \dontrun{
#'   gene_sets <- convert_gmt_to_list("c2.all.v2023.1.Hs.entrez.gmt")
#' }
convert_gmt_to_list <- function(file_path){
  # Read the file content
  lines <- readLines(file_path)
  
  # Pre-allocate the list with the number of lines in the file
  gene_sets <- vector("list", length(lines))
  
  # Loop over each line to process it
  for (i in seq_along(lines)) {
    # Split the line by tabs
    elements <- strsplit(lines[i], "\t")[[1]]
    
    # The first element is the name of the gene set
    gene_set_name <- elements[1]
    
    # The rest are the Entrez IDs (after the URL)
    entrez_ids <- elements[-(1:2)]
    
    # Store the gene set name and the vector of Entrez IDs
    names(gene_sets)[i] <- gene_set_name
    gene_sets[[i]] <- entrez_ids
  }
  
  gene_sets
}

#' Plot multi-dimensional scaling across successive PC pairs, coloured by a sample factor
#'
#' Computes an MDS layout (via \code{limma::plotMDS} on log-CPM values) for
#' four pairs of leading dimensions (1-2, 2-3, 3-4, 4-5) and assembles them
#' into a 2×2 patchwork, with a shared colour legend for the specified grouping
#' variable.
#'
#' @param data An \code{edgeR} DGEList object containing raw counts and a
#'   \code{samples} data frame with a \code{sample.id} column.
#' @param factor Character string; the name of a column in \code{data$samples}
#'   to use for point colour. Passed to \code{eval(parse(text = ...))} so
#'   quoted column names are accepted (e.g. \code{"Group"}).
#' @param lab Character string; legend title for the colour aesthetic.
#'
#' @return A patchwork object containing four ggplot panels.
#'
#' @examples
#' \dontrun{
#'   mds_by_factor(dge, factor = "Group", lab = "Treatment group")
#' }
mds_by_factor <- function(data, factor, lab){
  dims <- list(c(1,2), c(2:3), c(3,4), c(4,5))
  p <- vector("list", length(dims))
  
  for(i in 1:length(dims)){
    
    mds <- limma::plotMDS(edgeR::cpm(data, 
                                     log = TRUE), 
                          gene.selection = "common",
                          plot = FALSE, dim.plot = dims[[i]])
    
    data.frame(x = mds$x, 
               y = mds$y,
               sample = rownames(mds$distance.matrix.squared)) %>%
      left_join(rownames_to_column(data$samples, var = "sample")) -> dat
    
    p[[i]] <- ggplot(dat, aes(x = x, y = y, 
                              colour = eval(parse(text=(factor))))) +
      geom_point(size = 3) +
      ggrepel::geom_text_repel(aes(label = sample.id),
                               size = 2) +
      labs(x = glue("Principal Component {dims[[i]][1]}"),
           y = glue("Principal Component {dims[[i]][2]}"),
           colour = lab) +
      theme(legend.direction = "horizontal",
            legend.text = element_text(size = 8),
            legend.title = element_text(size = 9),
            axis.text = element_text(size = 8),
            axis.title = element_text(size = 9)) -> p[[i]]
  }
  
  wrap_plots(p, ncol = 2) + 
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
}

#' Volcano plot for differential expression results
#'
#' Draws a volcano plot (-log10 p-value vs log2 fold-change) for the top
#' table returned by edgeR, colouring points by DE direction (up/down/not
#' significant) and labelling genes that pass the FDR threshold.
#'
#' @param top Data frame of DE results as returned by \code{edgeR::topTags}
#'   (rownames are gene symbols).
#' @param cutoff Numeric; FDR threshold for significance (e.g. \code{0.05}).
#' @param dt Matrix returned by \code{edgeR::decideTests}; used to classify
#'   genes as up-regulated (+1), down-regulated (-1), or not significant (0).
#' @param pval_col Character string; name of the p-value column in \code{top}
#'   to plot on the y-axis (e.g. \code{"PValue"}).
#' @param fdr_col Character string; name of the adjusted p-value column used
#'   for the significance colour and gene labelling (e.g. \code{"FDR"}).
#' @param pal Character vector of colours for DE status levels (up, down,
#'   not significant).
#'
#' @return A ggplot object.
top_deg_volcano <- function(top, cutoff, dt, pval_col, fdr_col, pal){
  
  top %>% 
    mutate(sig = ifelse(!!sym(fdr_col) <= cutoff, glue("<= {cutoff}"), 
                        glue("> {cutoff}")))  %>%
    rownames_to_column(var = "SYMBOL") %>%
    left_join(dt[,1] %>% 
                data.frame %>%
                rownames_to_column(var = "SYMBOL") %>%
                dplyr::rename(status = 2)) %>%
    mutate(status = ifelse(status == 1, "Up",
                           ifelse(status == -1, "Down",
                                  "NotSig"))) %>%
    mutate(status = as.factor(status)) %>%
    mutate(status = fct_relevel(status, "NotSig", after = Inf)) %>%
    ggplot(aes(x = logFC, y = -log10(!!sym(pval_col)), color = status)) +
    geom_point(alpha = 0.75) +
    ggrepel::geom_text_repel(data = function(x) subset(x, eval(parse(text = fdr_col)) < cutoff), 
                             aes(x = logFC, y = -log10(eval(parse(text = pval_col))), 
                                 label = SYMBOL), 
                             size = 2, colour = "black", max.overlaps = 15) +
    labs(x = expression(~log[2]~"(Fold Change)"), 
         y = expression(~-log[10]~"(P-value)"),
         colour = glue("FDR < {cutoff}")) +
    scale_colour_manual(values = pal) +
    theme_classic() +
    theme(legend.position = "bottom")
}

#' Strip chart of normalised expression for top differentially expressed genes
#'
#' Plots per-sample log-CPM expression (both raw and TMM-normalised) as a
#' jittered strip chart for the top DE genes from a given contrast, faceted
#' by gene. A horizontal bar indicates the group mean. Only groups involved
#' in the contrast are shown.
#'
#' @param raw_counts DGEList or count matrix of raw counts.
#' @param norm_counts DGEList or count matrix of normalised counts (e.g.
#'   post-\code{calcNormFactors}).
#' @param group_info Data frame mapping sample identifiers to group labels;
#'   must contain a \code{Group} column and a column matching the sample names
#'   in the count matrices.
#' @param contr Single-column contrast matrix (genes × 1) defining the
#'   comparison; groups with non-zero weights are displayed.
#' @param top Data frame of DE results (rownames are gene identifiers); genes
#'   are taken from the first \code{num} rows that also pass \code{cutoff}.
#' @param num Integer; maximum number of genes to plot (default \code{9}).
#' @param severity Logical; if \code{TRUE} severity suffixes (\code{.M},
#'   \code{.S}) are retained in group labels, otherwise they are stripped.
#'
#' @return A ggplot object.
top_deg_stripchart <- function(raw_counts, norm_counts, group_info, contr, top, num = 9,
                               severity = FALSE){
  # plot up to top X DGE
  grps <- names(contr[,1])[abs(contr[,1]) > 0]
  
  edgeR::cpm(raw_counts, log = TRUE) %>% 
    data.frame %>%
    rownames_to_column(var = "gene") %>%
    pivot_longer(-gene, 
                 names_to = "sample", 
                 values_to = "raw") %>%
    inner_join(edgeR::cpm(norm_counts, log = TRUE) %>% 
                 data.frame %>%
                 rownames_to_column(var = "gene") %>%
                 pivot_longer(-gene, 
                              names_to = "sample", 
                              values_to = "norm")) %>%
    left_join(group_info) %>%
    left_join(top %>% 
                rownames_to_column(var = "gene")) %>%
    dplyr::filter(Group %in% grps,
                  gene %in% rownames(top)[1:min(num, max(which(top$FDR < cutoff)))]) %>%
    mutate(Group = case_when(severity ~ Group,
                             !severity ~ str_remove_all(Group, "\\.M$|\\.S$"))) %>% 
    ggplot(aes(x = Group,
               y = norm,
               colour = Group)) +
    geom_jitter(stat = "identity",
                width = 0.15,
                size = 1.25) +
    stat_summary(geom = "point",
                 fun.y = "mean",
                 col = "black",
                 shape = "_",
                 size = 14) +
    geom_jitter(aes(x = Group,
                    y = raw), stat = "identity",
                width = 0.15,
                size = 1.25, 
                alpha = 0.2) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     vjust = 0.5),
          legend.position = "bottom",
          legend.direction = "horizontal",
          strip.text = element_text(size = 7),
          axis.text.y = element_text(size = 6)) +
    labs(x = "Group", y = "log2 CPM") +
    facet_wrap(~gene, scales = "free_y") + 
    scale_color_brewer(palette = "Set2") +
    ggtitle(colnames(contr))
}

#' Dot plot of top CAMERA competitive gene set results
#'
#' Visualises the top \code{num} gene sets from each entry in a named list of
#' CAMERA results as a dot plot, with point size proportional to gene set size
#' and x-axis showing -log10(FDR). A dashed line marks the FDR = 0.05
#' threshold. Results are faceted by comparison type.
#'
#' @param results_list Named list of CAMERA result data frames (as returned by
#'   \code{limma::camera} or \code{limma::cameraPR}), one per contrast or
#'   gene set collection.
#' @param num Integer; maximum number of top sets to show per facet
#'   (default \code{10}).
#'
#' @return A ggplot object.
top_camera_sets <- function(results_list, num = 10){
  
  pal <- c(paletteer::paletteer_d("RColorBrewer::Set1")[2:1], "grey") 
  
  lapply(seq_along(results_list), function(i){
    results_list[[i]] %>%
      data.frame %>%
      dplyr::slice(1:min(num, n())) %>%
      rownames_to_column(var = "Set") %>%
      mutate(Type = names(results_list)[i],
             Rank = 1:min(num, n()))
  }) %>%
    bind_rows  %>%
    mutate(Set = str_wrap(str_replace_all(Set, "_", " "), width = 75),
           Set = str_remove_all(Set, "GO |REACTOME |HALLMARK |WP ")) %>%
    ggplot(aes(x = -log10(FDR), y = fct_reorder(Set, -Rank),
               colour = Direction)) +
    geom_point(aes(size = NGenes)) +
    facet_wrap(~Type, ncol = 1, scales = "free_y") +
    geom_vline(xintercept = -log10(0.05),
               linetype = "dashed")  +
    scale_colour_manual(values = pal) +
    labs(y = "Gene set", size = "Set size") +
    theme_classic(base_size = 10) +
    ggtitle("Camera gene set analysis")
}

#' Dot plot of top CAMERA results, faceted by cell type and comparison
#'
#' Similar to \code{top_camera_sets} but labels each facet with both the
#' comparison name and the current cell type (via the \code{cell} variable in
#' the calling environment). Gene set source prefixes (GO, REACTOME, etc.) are
#' replaced with compact bracketed abbreviations. Y-axis positions are set by
#' rank so labels remain readable when sets differ across facets.
#'
#' @param results_list Named list of CAMERA result data frames.
#' @param num Integer; maximum number of top sets per facet (default \code{10}).
#' @param wrap_width Integer; character width at which gene set names are
#'   wrapped (default \code{75}).
#' @param labeller Character string or labeller function passed to
#'   \code{facet_wrap} (default \code{"label_value"}).
#'
#' @return A ggplot object.
top_camera_sets_by_cell <- function(results_list, num = 10, wrap_width = 75,
                                    labeller = "label_value"){
  
  pal <- c(paletteer::paletteer_d("RColorBrewer::Set1")[2:1], "grey") 
  
  lapply(seq_along(results_list), function(i){
    results_list[[i]] %>%
      data.frame %>%
      dplyr::slice(1:min(num, n())) %>%
      rownames_to_column(var = "Set") %>%
      mutate(Type = glue("{names(results_list)[i]}; {cell}"))
  }) %>%
    bind_rows  %>%
    mutate(Set = str_replace_all(Set, "pulmonary_fibrosis_ctd", 
                                 "(CTD) PULMONARY FIBROSIS"),
           Set = str_wrap(str_replace_all(Set, "_", " "), width = wrap_width),
           Set = str_replace_all(Set, "GO", "(GO)"),
           Set = str_replace_all(Set, "REACTOME", "(R)"),
           Set = str_replace_all(Set, "WP", "(WP)"),
           Set = str_replace_all(Set, "HALLMARK", "(H)"),
           Rank = n():1) %>%
    # wrap in curly brackets so we can access the augmented dataset multiple times
    {
      ggplot(., aes(x = -log10(FDR), y = Rank,
                    colour = Direction)) +
        geom_point(aes(size = NGenes)) +
        facet_wrap(~Type, ncol = 1, scales = "free_y",
                   labeller = labeller) +
        geom_vline(xintercept = -log10(0.05),
                   linetype = "dashed")  +
        scale_colour_manual(values = pal) +
        scale_y_continuous(
          breaks = .$Rank,
          labels = .$Set,
          expand = c(0,.4)) +
        scale_size(range = c(1, 4)) +
        labs(y = "Gene set", size = "Set size") +
        theme_classic(base_size = 10) +
        ggtitle("Camera gene set analysis")
    }
}


#' Dot plot of top over-representation analysis (ORA) results
#'
#' Visualises the top \code{num} gene sets per contrast from a named list of
#' ORA results. Point colour encodes the gene ratio (overlap / universe size)
#' on a plasma viridis scale; point size encodes gene set size. A dashed line
#' marks FDR = 0.05. Results are faceted by comparison type.
#'
#' @param results_list Named list of ORA result data frames (e.g. as returned
#'   by \code{gene_set_test_ora}).
#' @param num Integer; maximum number of top sets to display per facet
#'   (default \code{10}).
#'
#' @return A ggplot object.
top_ora_sets <- function(results_list, num = 10){
  
  lapply(seq_along(results_list), function(i){
    results_list[[i]] %>%
      data.frame %>%
      dplyr::slice(1:min(num, n())) %>%
      rownames_to_column(var = "Set") %>%
      mutate(Type = names(results_list)[i],
             Rank = 1:min(num, n()))
  }) %>%
    bind_rows  %>%
    mutate(Set = str_wrap(str_replace_all(Set, "_", " "), width = 75),
           Set = str_remove_all(Set, "GO |REACTOME |HALLMARK |WP ")) %>%
    ggplot(aes(x = -log10(FDR), y = fct_reorder(Set, -Rank),
               colour = GR)) +
    geom_point(aes(size = N)) +
    facet_wrap(~Type, ncol = 1, scales = "free_y") +
    geom_vline(xintercept = -log10(0.05),
               linetype = "dashed")  +
    scale_colour_viridis_c(option = "plasma") +
    labs(y = "Gene set",
         colour = "Gene ratio",
         size = "Set size") +
    theme_classic(base_size = 10) +
    ggtitle("Over-representation gene set analysis")
}

#' Dot plot of top ORA results, faceted by cell type and comparison
#'
#' Analogous to \code{top_camera_sets_by_cell} but for ORA output. Gene ratio
#' (DE / N) is mapped to a continuous plasma colour scale. Facet labels
#' combine the comparison name and the cell type from the calling environment.
#'
#' @param results_list Named list of ORA result data frames.
#' @param num Integer; maximum number of top sets per facet (default \code{10}).
#' @param wrap_width Integer; character width for wrapping gene set names
#'   (default \code{75}).
#' @param labeller Character string or labeller function for \code{facet_wrap}
#'   (default \code{"label_value"}).
#'
#' @return A ggplot object.
top_ora_sets_by_cell <- function(results_list, num = 10, wrap_width = 75,
                                 labeller = "label_value"){
  
  lapply(seq_along(results_list), function(i){
    results_list[[i]] %>%
      data.frame %>%
      dplyr::slice(1:min(num, n())) %>%
      rownames_to_column(var = "Set") %>%
      mutate(Type = glue("{names(results_list)[i]}; {cell}"))
  }) %>%
    bind_rows %>%
    mutate(Set = str_replace_all(Set, "pulmonary_fibrosis_ctd", 
                                 "(CTD) PULMONARY FIBROSIS"),
           Set = str_wrap(str_replace_all(Set, "_", " "), width = wrap_width),
           Set = str_replace_all(Set, "GO", "(GO)"),
           Set = str_replace_all(Set, "REACTOME", "(R)"),
           Set = str_replace_all(Set, "WP", "(WP)"),
           Set = str_replace_all(Set, "HALLMARK", "(H)"),
           Rank = n():1) %>%
    # wrap in curly brackets so we can access the augmented dataset multiple times
    {
      ggplot(., aes(x = -log10(FDR), y = Rank,
                    colour = DE/N)) +
        geom_point(aes(size = N)) +
        facet_wrap(~Type, ncol = 1, scales = "free_y",
                   labeller = labeller) +
        geom_vline(xintercept = -log10(0.05),
                   linetype = "dashed")  +
        scale_colour_viridis_c(option = "plasma") +
        scale_y_continuous(
          breaks = .$Rank,
          labels = .$Set,
          expand = c(0,.4)) +
        scale_size(range = c(1, 4)) +
        labs(y = "Gene set",
             colour = "Gene ratio",
             size = "Set size") +
        theme_classic(base_size = 10) +
        ggtitle("Over-representation gene set analysis")
    }
}

#' Run over-representation analysis (ORA) across multiple gene set collections
#'
#' For each gene set collection in \code{gene_sets_list}, tests for enrichment
#' of \code{deg} (significant DE genes) against the full gene universe using
#' \code{missMethyl::gsaseq} / \code{topGSA}. Annotates the top 50 results
#' with the gene symbols of the overlapping DE genes, computes the gene ratio
#' (k/n), and writes a CSV for each collection to \code{cellDir}.
#'
#' @param gene_sets_list Named list of gene set collections; each element is
#'   itself a named list mapping gene set names to vectors of Entrez IDs.
#' @param deg Character vector of significant DE gene identifiers (Ensembl or
#'   the key type accepted by \code{gns}).
#' @param gns Named character vector mapping gene identifiers to Entrez IDs
#'   (used to convert \code{deg} and the full universe for ORA).
#' @param contr Single-column contrast matrix; used only to derive the output
#'   file name via \code{colnames(contr)}.
#' @param cellDir Character string; directory path where CSV result files are
#'   written.
#'
#' @return A named list (one element per collection) of ORA result data frames,
#'   each augmented with \code{DEG.GENES} (comma-separated symbols), \code{k}
#'   (overlap size), \code{n} (total DE genes), and \code{GR} (gene ratio).
gene_set_test_ora <- function(gene_sets_list, deg, gns, contr, cellDir){
  
  # GeneRatio = k/n
  # 
  # k is the overlap between your genes-of-interest and the geneset
  # n is the number of all unique genes-of-interest
  
  ora_list <- lapply(seq_along(gene_sets_list), function(i){
    topGSA(gsaseq(unname(gns[deg]),
                  universe = unname(gns),
                  collection = gene_sets_list[[i]],
                  plot.bias = FALSE),  
           n = 50) %>%
      data.frame -> tmp
    
    sig_genes <- unname(gns[deg])
    unlist(lapply(1:nrow(tmp), function(j){
      # Get gene symbols of significant genes
      # code adapted from Belinda Phipson missMethyl::gsameth 06/12/2024
      sig_genes_entrez <- sig_genes[sig_genes %in% gene_sets_list[[i]][[rownames(tmp)[j]]]]
      sig_genes_symbol <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, 
                                                                 keys = sig_genes_entrez,
                                                                 columns = "SYMBOL"))
      paste(sig_genes_symbol$SYMBOL,collapse=",")
    })) -> tmp$DEG.GENES
    
    tmp$k <- sapply(strsplit(tmp$DEG.GENES, ","), length)
    tmp$n <- length(deg)
    tmp$GR <- tmp$k/tmp$n
    
    data.table::fwrite(tmp %>%
                         data.frame %>%
                         rownames_to_column(var = "Set"),
                       file = file.path(cellDir, glue("ORA.{names(gene_sets_list[i])}.{colnames(contr)}.csv")),
                       sep = ",", quote = "auto")
    
    tmp
  })
  names(ora_list) <- names(gene_sets_list)
  
  ora_list
}

#' Run CAMERA competitive gene set testing across multiple collections
#'
#' For each gene set collection in \code{gene_sets_list}, maps gene sets to
#' indices in the fitted model via \code{limma::ids2indices} and runs a
#' pre-ranked CAMERA test (\code{limma::cameraPR}) using the supplied signed
#' likelihood ratio statistic. The top 50 results are annotated with DE gene
#' symbols and written to CSV.
#'
#' @param gene_sets_list Named list of gene set collections; each element maps
#'   gene set names to vectors of Entrez IDs.
#' @param gns Named character vector mapping gene identifiers (rownames of the
#'   fit) to Entrez IDs.
#' @param lrt Object returned by \code{edgeR::glmLRT}; rownames define the
#'   gene universe.
#' @param statistic Numeric vector of pre-ranked statistics (one per gene),
#'   e.g. signed LR statistics: \code{sign(logFC) * sqrt(LR)}.
#' @param contr Single-column contrast matrix; column name used in output
#'   filenames.
#' @param cellDir Character string; directory where CSV results are written.
#' @param deg Character vector of significant DE gene identifiers.
#'
#' @return A named list (one element per collection) of CAMERA result data
#'   frames augmented with a \code{DEG.GENES} column of comma-separated gene
#'   symbols.
gene_set_test_camera <- function(gene_sets_list, gns, lrt, statistic, contr, cellDir, deg){
  
  cam_list <- lapply(seq_along(gene_sets_list), function(i){
    id <- ids2indices(gene_sets_list[[i]], unname(gns[rownames(lrt)]))
    tmp <- cameraPR(statistic, id)
    tmp %>%
      data.frame %>%
      slice_head(n = 50) -> tmp
    
    sig_genes <- unname(gns[deg])
    unlist(lapply(1:nrow(tmp), function(j){
      # Get gene symbols of significant genes
      # code adapted from Belinda Phipson missMethyl::gsameth 06/12/2024
      sig_genes_entrez <- sig_genes[sig_genes %in% gene_sets_list[[i]][[rownames(tmp)[j]]]]
      sig_genes_symbol <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, 
                                                                 keys = sig_genes_entrez,
                                                                 columns = "SYMBOL"))
      paste(sig_genes_symbol$SYMBOL,collapse=",")
    })) -> tmp$DEG.GENES
    
    data.table::fwrite(tmp %>%
                         data.frame %>%
                         rownames_to_column(var = "Set"),
                       file = file.path(cellDir, glue("CAM.{names(gene_sets_list[i])}.{colnames(contr)}.csv")),
                       sep = ",", quote = "auto")
    
    tmp
  })
  names(cam_list) <- names(gene_sets_list)
  
  cam_list
}

#' Assemble a per-contrast DE summary panel (volcano + stripchart + ORA + CAMERA)
#'
#' Iterates over each contrast in \code{contr}, runs LRT differential
#' expression, saves the full results table to CSV, and—for contrasts with at
#' least one significant gene—performs ORA and CAMERA gene set tests, then
#' assembles a composite patchwork panel containing a volcano plot, strip
#' chart, ORA dot plot, and CAMERA dot plot.
#'
#' @param contr Contrast matrix (genes × contrasts) defining all comparisons.
#' @param cutoff Numeric FDR threshold for significance (e.g. \code{0.05}).
#' @param cellDir Character string; output directory for CSV result files.
#' @param gene_sets_list Named list of gene set collections passed to
#'   \code{gene_set_test_ora} and \code{gene_set_test_camera}.
#' @param gns Named Entrez ID vector for gene set mapping.
#' @param raw_counts Raw count DGEList for strip chart plotting.
#' @param norm_counts Normalised count DGEList for strip chart plotting.
#' @param group_info Sample-level metadata data frame.
#' @param layout patchwork layout specification string for panel arrangement.
#' @param pal Colour palette vector for volcano/strip chart DE direction.
#' @param severity Logical vector (length == \code{ncol(contr)}); passed to
#'   \code{top_deg_stripchart} for each contrast.
#' @param pval_col Character; p-value column name for the volcano y-axis
#'   (default \code{"PValue"}).
#' @param fdr_col Character; FDR column name for significance colouring
#'   (default \code{"FDR"}).
#'
#' @return A list of patchwork objects (one per contrast); elements for
#'   contrasts with no significant genes are \code{NULL}.
plot_ruv_results_summary <- function(contr, cutoff, cellDir, gene_sets_list, gns,
                                     raw_counts, norm_counts, group_info, layout, 
                                     pal, severity,
                                     pval_col = "PValue",
                                     fdr_col = "FDR"){
  p <- vector("list", ncol(contr))
  
  for(i in 1:(ncol(contr))){
    lrt <- glmLRT(fit, contrast = contr[,i])
    topTags(lrt, n = Inf) %>%
      data.frame -> top
    
    if(sum(top$FDR < cutoff) > 0){
      # top DGE results
      data.table::fwrite(top %>%
                           rownames_to_column(var = "gene"), 
                         file = file.path(cellDir, glue("{colnames(contr)[i]}.csv")),
                         sep = ",", quote = "auto")
      
      deg <- rownames(top)[top$FDR < cutoff]
      ora_list <- gene_set_test_ora(gene_sets_list, deg, gns, 
                                    contr[, i, drop = FALSE], cellDir)
      
      # run camera competitive gene set test  
      # use signed likelihood ratio test statistic as recommended by GS here: https://support.bioconductor.org/p/112937/
      statistic <- sign(lrt$table$logFC) * sqrt(lrt$table$LR)
      cam_list <- gene_set_test_camera(gene_sets_list, gns, lrt, statistic, 
                                       contr[, i, drop = FALSE], cellDir, deg)
      
      top_deg_volcano(top, cutoff, dt[[i]], pval_col, fdr_col, pal) -> p1
      top_deg_stripchart(raw_counts, 
                         norm_counts, 
                         group_info, 
                         contr = contr[, i, drop = FALSE], 
                         top = top, 
                         severity = severity[i],
                         num = 9) -> p2
      
      # over-representation analysis top 10 plots
      top_ora_sets(ora_list, num = 10) -> p3
      # camera top 10 plots
      top_camera_sets(cam_list, num = 10) -> p4
      
      p[[i]] <- wrap_elements(p1 + p2) + 
        wrap_elements(p3) + 
        wrap_elements(p4) +
        plot_layout(design = layout)
    }
  }
  
  p
}

#' Load custom gene lists from a CSV and map gene symbols to Entrez IDs
#'
#' Reads a wide-format CSV where each column represents a named gene set
#' (header = set name, rows = gene symbols). Pivots to long format, looks up
#' Entrez IDs via \code{org.Hs.eg.db}, and returns a named list of Entrez ID
#' vectors, one per gene set.
#'
#' @param file Character string; path to a comma-delimited CSV file with gene
#'   symbols arranged in columns (one column per gene set, column headers are
#'   set names).
#'
#' @return A named list where each element is a character vector of Entrez IDs
#'   for the corresponding gene set.
#'
#' @examples
#' \dontrun{
#'   custom_sets <- create_custom_gene_lists_from_file("gene_lists.csv")
#' }
create_custom_gene_lists_from_file <- function(file = file){
  
  read.csv2(file = file,
            header = TRUE, sep = ",") %>% 
    pivot_longer(cols = everything()) %>%
    mutate(value = str_trim(value),
           entrez = unname(unlist(AnnotationDbi::mapIds(org.Hs.eg.db, 
                                                        keys = value, 
                                                        column = c("ENTREZID"),
                                                        keytype = "SYMBOL",
                                                        multiVals = "first"))[value])) %>%
    dplyr::rename("symbol" = "value",
                  "source" = "name") -> tmp
  
  split(tmp$entrez, tmp$source)
}

#' UMAP coloured by cell type annotation with repelled cluster labels
#'
#' Draws a Seurat \code{DimPlot} UMAP coloured by the specified annotation
#' level, removes the legend, and overlays bold repelled cluster labels using
#' \code{LabelClusters}. Axis elements are suppressed for a clean plot.
#'
#' @param seu A Seurat object with a UMAP reduction computed.
#' @param ann_level Character string; name of the metadata column to use for
#'   colouring and labelling (e.g. \code{"cell_type_l2"}).
#' @param cluster_pal Character string; a paletteer palette identifier (e.g.
#'   \code{"ggsci::default_igv"}) used to colour clusters.
#' @param direction Integer; palette direction, \code{1} (default) or
#'   \code{-1} to reverse.
#'
#' @return A ggplot object with cluster labels overlaid.
draw_umap_with_labels <- function(seu, ann_level, cluster_pal, direction = 1){
  DimPlot(seu, 
          group.by = ann_level, label = F, repel = T,
          label.size = 3) +
    scale_color_paletteer_d(cluster_pal, direction = direction) +
    NoLegend() -> p1
  
  LabelClusters(p1, id = ann_level, repel = TRUE,
                size = 3, box = TRUE, fontfamily = "arial",
                fontface = "bold") +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          plot.title = element_blank())
  
}

#' Seurat dot plot of top marker genes, faceted by cluster with coloured strips
#'
#' Selects the top \code{num} marker genes per cluster (ranked by
#' \code{avg_log2FC} when available), draws a Seurat \code{DotPlot} with a
#' two-colour gradient, and facets columns by cluster using
#' \code{ggh4x::facet_grid2} with palette-matched strip backgrounds. Y-axis
#' labels and strip text can be remapped via \code{lab_map}.
#'
#' @param seu A Seurat object.
#' @param markers Data frame of marker genes as returned by
#'   \code{Seurat::FindAllMarkers}, containing at minimum \code{cluster} and
#'   \code{gene} columns.
#' @param ann_level Character string; metadata column used as the
#'   \code{group.by} axis in \code{DotPlot}.
#' @param cluster_pal Character string; paletteer palette for strip and point
#'   colours.
#' @param lab_map Named character vector mapping raw cluster/cell-type labels
#'   to display labels; used for both y-axis tick labels and facet strip text.
#' @param direction Integer; palette direction (default \code{1}).
#' @param num Integer; number of top markers to select per cluster
#'   (default \code{5}).
#' @param assay Character string; Seurat assay to use (default \code{"SCT"}).
#' @param strip.text.size Numeric; font size for facet strip labels
#'   (default \code{8}).
#' @param strip.text.blank Logical; if \code{TRUE} strip labels are blank
#'   (colours only), useful when labels are shown elsewhere (default
#'   \code{FALSE}).
#' @param strip.alpha Numeric in [0,1]; alpha applied to strip fill colours
#'   (default \code{0.75}).
#' @param dot.scale Numeric; passed to \code{Seurat::DotPlot} to control dot
#'   size scaling (default \code{4}).
#'
#' @return A ggplot object.
draw_marker_gene_dotplot <- function(seu, markers, ann_level, cluster_pal, 
                                     lab_map, direction = 1, num = 5, assay = "SCT",
                                     strip.text.size = 8, strip.text.blank = FALSE,
                                     strip.alpha = 0.75, dot.scale = 4){
  
  # --- top markers per cluster ---
  markers_top <- markers %>%
    distinct(cluster, gene, .keep_all = TRUE) %>%
    group_by(cluster) %>%
    {
      if ("avg_log2FC" %in% names(.)) {
        arrange(., desc(avg_log2FC), .by_group = TRUE)
      } else {
        .  # no sorting
      }
    } %>%
    slice_head(n = num) %>%
    ungroup()
  
  
  # features in cluster blocks (keeps your ranking)
  features_vec <- markers_top %>%
    {
      if ("avg_log2FC" %in% names(.)) {
        arrange(., cluster, desc(avg_log2FC), .by_group = TRUE)
      } else {
        .  # no sorting
      }
    } %>%
    #arrange(cluster, desc(avg_log2FC)) %>%
    pull(gene) %>%
    unique()
  
  # feature -> cluster lookup
  feat2clust <- markers_top %>%
    dplyr::select(gene, cluster) %>%
    distinct()
  
  cluster_labels <- unique(feat2clust$cluster)
  cluster_lvls <- as.character(cluster_labels[order(cluster_labels)])
  
  # build strip colours in that order
  pal_vec_raw <- paletteer::paletteer_d(cluster_pal, direction = direction)
  pal_vec <- setNames(rep_len(pal_vec_raw, length(cluster_lvls)), cluster_lvls)
  
  strip <- ggh4x::strip_themed(
    background_x = ggh4x::elem_list_rect(fill = scales::alpha(unname(pal_vec[cluster_lvls]), strip.alpha)),
    text_x       = ggh4x::elem_list_text(colour = "black", face = "bold")
  )
  
  # ---- DotPlot then inject faceting var ----
  DefaultAssay(seu) <- assay
  p <- Seurat::DotPlot(seu, features = features_vec, group.by = ann_level,
                       cols = c("#e0ecf4", "#88419d"), dot.scale = dot.scale)
  
  p$data <- p$data %>%
    left_join(feat2clust, by = c("features.plot" = "gene")) %>%
    mutate(cluster = factor(cluster, levels = cluster_lvls))
  
  # build a labeller 
  facet_labeller <- if (strip.text.blank) {
    # Safe blank labeller — returns empty strings (no element_blank() needed)
    labeller(cluster = as_labeller(function(x) rep("", length(x))))
  } else {
    # Your original labeller using lab_map
    labeller(cluster = as_labeller(lab_map, default = label_value))
  }
  
  p <- p +
    ggh4x::facet_grid2(
      ~ cluster,
      scales   = "free_x",
      space    = "free_x",
      strip    = strip,
      labeller = facet_labeller
    ) +
    scale_y_discrete(labels = function(x) lab_map[x]) +  # y-axis short names
    labs(x = NULL, y = NULL) +
    theme_bw(base_family = "Arial") +
    theme(
      axis.text.x   = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
      axis.text.y   = element_text(size = 8),
      legend.text   = element_text(size = 8),
      legend.title  = element_text(size = 9),
      legend.position = "bottom",
      strip.text.x  = element_text(size = strip.text.size, 
                                   margin = margin(2, 2, 2, 2)),
      strip.background = element_rect(colour = NA),
      panel.spacing = unit(3, "pt"),
      axis.ticks    = element_blank(),
      axis.line     = element_blank(),
      panel.border = element_blank()
    )
  
  p
  
}

#' Stacked bar chart of cell type proportions per sample, faceted by condition
#'
#' Computes arcsine-transformed cell type proportions via
#' \code{propeller::getTransformedProps}, then draws a stacked bar chart with
#' samples on the x-axis, coloured by cell type, and faceted by condition
#' group using \code{ggh4x::facet_grid2} with palette-matched coloured strip
#' backgrounds. A colour legend for condition groups is generated via an
#' invisible point layer.
#'
#' @param seu A Seurat object with \code{sample.id} and \code{Group} metadata
#'   columns.
#' @param ann_level Character string; metadata column containing cell type
#'   labels.
#' @param cluster_pal Character string; paletteer palette for cell type fill
#'   colours.
#' @param strip_colours Named character vector mapping group names to hex
#'   colours for facet strip backgrounds; names must match levels of
#'   \code{seu$Group}.
#' @param direction Integer; palette direction for \code{cluster_pal}
#'   (default \code{1}).
#'
#' @return A ggplot object.
draw_cell_type_proportions_barplot <- function(seu, ann_level, cluster_pal,
                                               strip_colours, direction = 1){
  
  props <- getTransformedProps(clusters = seu@meta.data[,ann_level],
                               sample = seu$sample.id, transform="asin")
  
  strip <- ggh4x::strip_themed(
    background_x = ggh4x::elem_list_rect(
      fill = unname(strip_colours)
    ))
  
  facet_labeller <- labeller(Group = as_labeller(function(x) rep("", length(x))))
  
  # Keep Group order consistent with your strip_colours vector
  group_lvls <- intersect(names(strip_colours), unique(seu@meta.data$Group))
  
  # ---- Legend dummy data: one (sample, y=0) row per Group ----
  legend_df <- seu@meta.data %>%
    dplyr::group_by(Group) %>%
    dplyr::summarise(
      sample = dplyr::first(sample.id),  # pick any sample present in this Group
      Freq   = 0,
      .groups = "drop"
    ) %>%
    dplyr::mutate(Group = factor(Group, levels = group_lvls))
  
  props$Proportions %>%
    data.frame %>%
    inner_join(seu@meta.data %>%
                 dplyr::select(sample.id,
                               Group),
               by = c("sample" = "sample.id")) %>%
    distinct() %>%
    ggplot(aes(x = sample, y = Freq, fill = clusters)) +
    geom_bar(stat = "identity", color = "black", size = 0.1) +
    # ---- legend-only layer for facet colours ----
  #invisible points create a legend for `colour = Group`
  geom_point(data = legend_df,
             aes(x = sample, y = Freq, colour = Group),
             size = 0, alpha = 0, inherit.aes = FALSE) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45,
                                     vjust = 1,
                                     hjust = 1,
                                     size = 8),
          axis.title = element_text(size = 9),
          strip.text = element_text(size = 8),
          strip.background = element_rect(colour = NA),
          panel.spacing = unit(2, "points", data = NULL),
          plot.margin = unit(c(0, 0, 0, 0), "lines"),
          legend.text = element_text(size = 8),
          legend.key.size = unit(0.8, "lines"),
          legend.title = element_text(size = 9)) +
    labs(y = "Cell type proportion", fill = "Cell type", x = "Sample") +
    scale_fill_paletteer_d(cluster_pal, direction = direction) +
    # facet-colour legend (groups) — matches strip colours
    scale_colour_manual(
      name   = "Condition",
      values = strip_colours[unique(seu$Group)],
      drop   = FALSE
    ) +
    guides(
      colour = guide_legend(
        override.aes = list(
          shape = 15,        # filled square
          size = 5,          # legend key size
          alpha = 1         # fully opaque
        ),
        order = 1
      ),
      fill   = guide_legend(order = 2)
    ) +
    ggh4x::facet_grid2(~Group, 
                       scales = "free_x", 
                       space = "free_x",
                       strip = strip,
                       labeller = facet_labeller) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.02)))
}

#' Heatmap of top DE genes for a given contrast
#'
#' Filters log-CPM expression to the relevant sample groups and top 50 DE
#' genes (by FDR rank), then draws a row-scaled heatmap using
#' \code{tidyHeatmap::heatmap} with a red-white-blue diverging palette.
#' Columns are not clustered; the group axis is determined automatically from
#' the contrast name.
#'
#' @param top Data frame of DE results (rownames are gene IDs); only the first
#'   \code{min(50, nrow(top))} genes are shown.
#' @param comparison Character string; contrast name used to parse which sample
#'   groups to include (e.g. \code{"IVA-CF"}).
#' @param counts DGEList or count matrix from which log-CPM values are derived.
#' @param sample_data Data frame (or DataFrame) of sample metadata containing
#'   at minimum \code{sample.id}, \code{Group}, \code{Group_severity},
#'   \code{Severity}, and \code{Age} columns.
#'
#' @return A \code{tidyHeatmap} / ComplexHeatmap object, or \code{NULL}
#'   invisibly if \code{top} has zero rows.
top_deg_heatmap <- function(top, comparison, counts, sample_data){
  
  if(nrow(top) > 0) {
    groups <- unlist(str_split(str_remove_all(comparison, "-?0.5\\*|-?1\\*"), " "))
    if(any(str_detect(groups, "IVA"))) groups <- unique(str_remove_all(groups, "\\.M|\\.S"))
    group_column <- ifelse(all(groups %in% ySub$samples$Group),
                           "Group",
                           "Group_severity")
    
    edgeR::cpm(counts, log = TRUE) %>% 
      data.frame %>%
      rownames_to_column(var = "Gene") %>%
      pivot_longer(-Gene, 
                   names_to = "Sample", 
                   values_to = "adj. cpm scaled") %>%
      left_join(sample_data %>% 
                  data.frame %>%
                  dplyr::select(sample.id, 
                                Group,
                                Group_severity,
                                Severity,
                                Age),
                by = c("Sample" = "sample.id")) %>%
      dplyr::filter(!!sym(group_column) %in% groups) %>%
      dplyr::filter(Gene %in% rownames(top)[1:min(50, nrow(top))]) %>%
      group_by(!!sym(group_column)) -> dat
    
    tidyHeatmap::heatmap(
      .data = dat,
      .row = Gene,  
      .column = Sample,
      .value = `adj. cpm scaled`,
      scale = "row",
      cluster_columns = FALSE,
      row_names_gp = gpar(fontsize = 8),
      column_names_gp = gpar(fontsize = 8),
      column_title_gp = gpar(fontsize = 10),
      row_title_gp = gpar(fontsize = 10),
      heatmap_legend_param = list(labels_gp = gpar(fontsize = 7),
                                  title_gp = gpar(fontsize = 8, 
                                                  fontface = 'bold')),
      #palette_value = c("red", "white", "blue"))
      palette_value = circlize::colorRamp2(
        seq(-4, 4, length.out = 11),
        RColorBrewer::brewer.pal(11, "RdBu")))
  }
}

#' Aggregate significant DE results across cell types from pre-fitted edgeR objects
#'
#' Loads a list of RDS files (each containing a fitted edgeR model and contrast
#' matrix), runs \code{glmTreat} for the named contrast, and row-binds the
#' significant results across all cell types into a single data frame. Cell
#' type identity is parsed from the filename prefix. Results are filtered to
#' \code{FDR < cutoff}, joined to cell-frequency weights, and sorted by
#' descending cell frequency then log fold-change.
#'
#' @param files Character vector of paths to RDS files. Each file must contain
#'   a list with elements \code{$fit} (a fitted DGEGLM) and \code{$contr} (a
#'   contrast matrix); the cell type label is extracted from the filename
#'   prefix (everything before the first \code{.}).
#' @param cont_name Character string; name of the contrast column in
#'   \code{$contr} to test.
#' @param cell_freq Data frame with a \code{cell} column (matching filename
#'   prefixes) and an \code{n} column used to weight the sort order.
#' @param treat_lfc Numeric; minimum absolute log2 fold-change passed to
#'   \code{glmTreat} (default \code{0}, equivalent to \code{glmLRT}).
#' @param cutoff Numeric; FDR threshold for filtering significant genes
#'   (default \code{0.05}).
#' @param suffix Character string; filename suffix used in documentation only
#'   (default \code{".all_samples.fit.rds"}).
#'
#' @return A data frame with columns \code{gene}, \code{sig} (DE direction),
#'   \code{cell}, \code{logFC}, and \code{FDR}, restricted to significant
#'   results and sorted by cell frequency then \code{logFC}.
get_deg_data <- function(files, cont_name, cell_freq, treat_lfc = 0, 
                         cutoff = 0.05, suffix = ".all_samples.fit.rds"){
  
  bind_rows(lapply(files, function(f){
    
    deg_results <- readRDS(f)
    lrt <- glmTreat(deg_results$fit, 
                    contrast = deg_results$contr[,cont_name],
                    lfc = treat_lfc)
    top <- as.data.frame(topTags(lrt, n = Inf))
    cell <- str_extract(basename(f), "^[^.]+")
    dt <- decideTests(lrt, p.value = 0.05)
    
    as.data.frame(dt) %>% 
      rename_with(.cols = 1, ~"sig") %>%
      rownames_to_column(var = "gene") %>%
      mutate(cell = cell) %>%
      left_join(top %>%
                  rownames_to_column(var = "gene") %>%
                  dplyr::select(gene, logFC, FDR))
  })) %>%
    dplyr::filter(FDR < cutoff) %>%
    left_join(cell_freq) %>%
    arrange(-n, logFC) %>%
    dplyr::select(-n)
  
}  


#' Volcano plot overlaying standard LRT and glmTreat significance for a single cell type
#'
#' Loads a pre-fitted edgeR object for a given cell type, runs both
#' \code{glmLRT} and \code{glmTreat} for the contrast named by the
#' \code{cont_name} variable in the calling environment, and draws a volcano
#' plot that distinguishes four gene categories: not significant (grey);
#' significant by LRT only (coloured outline); significant by both LRT and
#' \code{glmTreat} (filled, labelled). This layered display makes the
#' additional stringency of the fold-change threshold immediately visible.
#'
#' @param cell Character string; cell type label used to construct the RDS
#'   file path as \code{here("data", "intermediate_objects", paste0(cell, suffix))}.
#' @param suffix Character string; filename suffix appended to \code{cell}
#'   (e.g. \code{".all_samples.fit.rds"}).
#' @param cutoff Numeric; FDR threshold for the significance dashed line and
#'   \code{decideTests} calls (default \code{0.05}).
#' @param lfc_cutoff Numeric; minimum absolute log2 fold-change for
#'   \code{glmTreat}; \code{0} reduces to standard LRT (default \code{0}).
#'
#' @return A ggplot object.
#'
#' @note Relies on \code{cont_name} being defined in the calling environment.
draw_treat_volcano_plot <- function(cell, suffix, cutoff = 0.05, lfc_cutoff = 0){
  f <- here("data",
            "intermediate_objects",
            glue("{cell}{suffix}"))
  
  deg_results <- readRDS(f)
  lrt <- glmLRT(deg_results$fit, 
                contrast = deg_results$contr[, cont_name])
  treat_lrt <- glmTreat(deg_results$fit, 
                        contrast = deg_results$contr[, cont_name],
                        lfc = lfc_cutoff)
  dat <- data.frame(gene = rownames(lrt),
                    logCPM = lrt$table$logCPM, 
                    logFC = lrt$table$logFC,
                    p = p.adjust(lrt$table$PValue, method = "BH"),
                    p_treat = p.adjust(treat_lrt$table$PValue, method = "BH"),
                    dt = as.vector(decideTests(lrt, p.value = cutoff)),
                    dt_treat = as.vector(decideTests(treat_lrt, p.value = cutoff)))
  
  pal_dt <- paletteer::paletteer_d("RColorBrewer::Set1") 
  
  ggplot(dat, aes(x = logFC, y = -log10(p))) +
    geom_point(colour = "lightgrey", shape = 1) +
    geom_point(data = dat[dat$dt > 0,], colour = pal_dt[1], shape = 1) +
    geom_point(data = dat[dat$dt < 0,], colour = pal_dt[2], shape = 1) +
    geom_point(data = dat[dat$dt_treat > 0,], colour = pal_dt[1]) +
    geom_point(data = dat[dat$dt_treat < 0,], colour = pal_dt[2]) +
    ggrepel::geom_text_repel(data = dat[dat$dt_treat > 0,], 
                             aes(label = gene), size = 2,
                             colour = pal_dt[1]) +
    ggrepel::geom_text_repel(data = dat[dat$dt_treat < 0,], 
                             aes(label = gene), size = 2,
                             colour = pal_dt[2]) +
    geom_hline(yintercept = -log10(cutoff), linetype = "dashed") +
    theme_classic() +
    labs(x = "log2(FC)",
         y = "-log10(FDR)")
  
}