
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

top_deg_volcano <- function(top, cutoff, dt, pval_col, fdr_col, pal){
  
  # pal <- ifelse(sum(summary(dt)[c(1,3),]) < 2,
  #               ifelse(summary(dt)[1,] == 1, pal[-2], pal[-1]),
  #               pal)
  
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

top_deg_stripchart <- function(raw_counts, norm_counts, group_info, contr, top, num = 9){
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
    dplyr::filter(Group %in% grps,
                  gene %in% rownames(top)[1:min(num, max(which(top$FDR < cutoff)))]) %>%
    mutate(Group = ifelse(str_detect(Group, str_remove(grps[1], "CF.")),
                          grps[1], 
                          grps[2])) %>%
    ggplot(aes(x = Group,
               y = norm,
               colour = Group)) +
    geom_boxplot(outlier.shape = NA, colour = "grey") +
    geom_jitter(stat = "identity",
                width = 0.15,
                size = 1.25) +
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

top_camera_sets <- function(results_list, num = 10){
  
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
    labs(y = "Gene set") +
    theme_classic(base_size = 10) +
    ggtitle("Camera gene set analysis")
}

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
               colour = DE/N*100)) +
    geom_point(aes(size = N)) +
    facet_wrap(~Type, ncol = 1, scales = "free_y") +
    geom_vline(xintercept = -log10(0.05),
               linetype = "dashed")  +
    scale_colour_viridis_c(option = "plasma") +
    labs(y = "Gene set",
         colour = "% DEGs in set",
         size = "Set size") +
    theme_classic(base_size = 10) +
    ggtitle("Over-representation gene set analysis")
}

gene_set_test_ora <- function(gene_sets_list, deg, gns, contr, cellDir){

  ora_list <- lapply(seq_along(gene_sets_list), function(i){
    topGSA(gsaseq(unname(gns[deg]),
                  universe = unname(gns),
                  collection = gene_sets_list[[i]],
                  plot.bias = FALSE)) -> tmp
    
    write.table(tmp %>%
                  data.frame %>%
                  rownames_to_column(var = "Set"),
                file = file.path(cellDir, glue("ORA.{names(gene_sets_list[i])}.{colnames(contr)}.csv")),
                sep = ",", quote = F, col.names = NA)
    
    tmp
  })
  names(ora_list) <- names(gene_sets_list)
  
  ora_list
}

gene_set_test_camera <- function(gene_sets_list, gns, lrt, statistic, contr, cellDir){
  
  cam_list <- lapply(seq_along(gene_sets_list), function(i){
    id <- ids2indices(gene_sets_list[[i]], unname(gns[rownames(lrt)]))
    tmp <- cameraPR(statistic, id)
    
    write.table(tmp %>%
                  data.frame %>%
                  rownames_to_column(var = "Set"),
                file = file.path(cellDir, glue("CAM.{names(gene_sets_list[i])}.{colnames(contr)}.csv")),
                sep = ",", quote = F, col.names = NA)
    
    tmp
  })
  names(cam_list) <- names(gene_sets_list)
  
  cam_list
}

plot_ruv_results_summary <- function(contr, cutoff, cellDir, gene_sets_list, gns,
                                     raw_counts, norm_counts, group_info, layout, 
                                     pal,
                                     pval_col = "PValue",
                                     fdr_col = "FDR"){
  p <- vector("list", ncol(contr))
  
  for(i in 1:(ncol(contr))){
    lrt <- glmLRT(fit, contrast = contr[,i])
    topTags(lrt, n = Inf) %>%
      data.frame -> top
    
    if(sum(top$FDR < cutoff) > 0){
      # top DGE results
      write.table(top, 
                  file = file.path(cellDir, glue("{colnames(contr)[i]}.csv")),
                  sep = ",", quote = F, col.names = NA)
      
      deg <- rownames(top)[top$FDR < cutoff]
      ora_list <- gene_set_test_ora(gene_sets_list, deg, gns, 
                                    contr[, i, drop = FALSE], cellDir)
      
      # run camera competitive gene set test  
      # use signed likelihood ratio test statistic as recommended by GS here: https://support.bioconductor.org/p/112937/
      statistic <- sign(lrt$table$logFC) * sqrt(lrt$table$LR)
      cam_list <- gene_set_test_camera(gene_sets_list, gns, lrt, statistic, 
                                       contr[, i, drop = FALSE], cellDir)
      
      top_deg_volcano(top, cutoff, dt[[i]], pval_col, fdr_col, pal) -> p1
      top_deg_stripchart(raw_counts, 
                         norm_counts, 
                         group_info, 
                         contr = contr[, i, drop = FALSE], 
                         top = top, 
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

draw_umap_with_labels <- function(seu, ann_level, cluster_pal, direction = 1){
  DimPlot(seu, 
          group.by = ann_level, label = F, repel = T,
          label.size = 3) +
    scale_color_paletteer_d(cluster_pal, direction = direction) +
    NoLegend() -> p1
  
  LabelClusters(p1, id = ann_level, repel = TRUE,
                size = 3, box = TRUE, fontfamily = "arial") +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          plot.title = element_blank())
  
}

draw_marker_gene_dotplot <- function(seu, rds_path, ann_level, cluster_pal, direction = 1){
  markers <- readRDS(rds_path)
  
  markers %>%
    group_by(cluster) %>%
    slice_head(n = 10) %>%
    mutate(cluster = as.character(cluster)) %>%
    ungroup() %>%
    dplyr::arrange(cluster, .by_group = FALSE) -> markers
  
  d <- duplicated(markers$gene)
  markers[!d,] %>%
    group_by(cluster) %>%
    slice_head(n = 5) -> top
  
  pal <- setNames(paletteer::paletteer_d(cluster_pal, direction = direction),
                  unique(markers$cluster))
  cell_type_cols <- pal[top$cluster]
  
  DefaultAssay(seu) <- "RNA"
  strip <- strip_themed(background_x = elem_list_rect(fill = unique(cell_type_cols)))
  DotPlot(seu,
          features = top$gene,
          group.by = ann_level,
          cols = c("azure1", "blueviolet"),
          dot.scale = 3,
          assay = "SCT") +
    FontSize(x.text = 9, y.text = 9) +
    labs(y = element_blank(), x = element_blank()) +
    facet_grid2(~top$cluster,
                scales = "free_x",
                space = "free_x",
                 strip = strip) +
    theme(axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     vjust = 0.5,
                                     size = 8),
          axis.text.y = element_text(size = 8),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 9),
          legend.position = "bottom",
          strip.text = element_text(size = 0),
          text = element_text(family = "arial"),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          panel.spacing = unit(2, "points"))
}

draw_cell_type_proportions_barplot <- function(seu, ann_level, cluster_pal,
                                               direction = 1){
  props <- getTransformedProps(clusters = seu@meta.data[,ann_level],
                               sample = seu$sample.id, transform="asin")
  
  props$Proportions %>%
    data.frame %>%
    inner_join(seu@meta.data %>%
                 dplyr::select(sample.id,
                               Group),
               by = c("sample" = "sample.id")) %>%
    distinct() %>%
    ggplot(aes(x = sample, y = Freq, fill = clusters)) +
    geom_bar(stat = "identity", color = "black", size = 0.1) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45,
                                     vjust = 1,
                                     hjust = 1,
                                     size = 8),
          axis.title = element_text(size = 9),
          strip.text = element_text(size = 8),
          #strip.background = element_blank(),
          panel.spacing = unit(2, "points", data = NULL),
          plot.margin = unit(c(0, 0, 0, 0), "lines"),
          legend.text = element_text(size = 8),
          legend.key.size = unit(1, "lines"),
          legend.title = element_text(size = 9)) +
    labs(y = "Cell type proportion", fill = "Cell type", x = "Sample") +
    scale_fill_paletteer_d(cluster_pal, direction = direction) +
    facet_grid(~Group, scales = "free_x", space = "free_x")
}
