# Create SingleCellExperiment object from 10X data and identify empty droplets
# using DropletUtils.
# Peter Hickey
# 2021-11-24
# Modified by Jovana Maksimovic
# 2024-02-19

# Setup ------------------------------------------------------------------------

library(DropletUtils)
library(here)
library(BiocParallel)

batchDirs <- list.files(here("data"), pattern = "C133_Neeland_batch", 
                        full.names = TRUE) 

for (dir in batchDirs) {
  batch <- basename(dir)
  
  emptyDir <- here("data", batch, "data", "emptyDrops")
  if(!dir.exists(emptyDir)) dir.create(emptyDir, recursive = TRUE)
  
  # Load SingleCellExperiment object ----------------------------------------
  
  sce_raw <- readRDS(here("data",
                          batch,
                          "data",
                          "SCEs",
                          paste0(batch,".CellRanger.SCE.rds")))
  
  # Identify empty droplets ------------------------------------------------------
  
  capture_names <- unique(sce_raw$Sample)
  capture_names <- setNames(capture_names, capture_names)
  
  set.seed(100)
  list_of_empties <- lapply(capture_names, function(cn) {
    message(cn)
    emptyDrops(counts(sce_raw)[, sce_raw$Sample == cn])
  })
  
  # Check if more permutations are needed; see
  # https://osca.bioconductor.org/quality-control.html#testing-for-empty-droplets
  more_permutations_needed <- sapply(list_of_empties, function(e) {
    table(
      Sig = e$FDR <= 0.001,
      Limited = e$Limited)[1, 2] > 0
  })
  stopifnot(all(!more_permutations_needed))
  
  # Save outputs -----------------------------------------------------------------
  
  for (cn in capture_names) {
    message(cn)
    empties <- list_of_empties[[cn]]
    saveRDS(
      object = empties,
      file = here(
        "data",
        batch,
        "data",
        "emptyDrops",
        paste0(cn, ".emptyDrops.rds")))
    
    writeLines(
      text = sce_raw[["Barcode"]][sce_raw$Sample == cn][which(empties$FDR <= 0.001)],
      con = here(
        "data",
        batch,
        "data",
        "emptyDrops",
        paste0(cn, ".barcodes.txt")))
  }
}