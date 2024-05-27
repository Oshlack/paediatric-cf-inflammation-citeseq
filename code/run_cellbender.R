#!/usr/bin/env Rscript

captures <- c(list.files(here::here("data",
                              "C133_Neeland_batch0",
                              "data",
                              "190930_A00152_0150_BHTYCMDSXX",
                              "GE"), full.names = TRUE),
              list.files(here::here("data",
                                    "C133_Neeland"), full.names = TRUE))

files <- file.path(captures,
                "outs",
                "multi",
                "count",
                "raw_feature_bc_matrix.h5")

jobDir <- here::here()
outDirs <- dirname(files)

for (i in 1:length(files)){
  # Start writing to this file
  jobFile <- glue::glue("{jobDir}/cn_{basename(captures[i])}.job")
  sink(file = jobFile)
  
  # the basic job submission script is a bash script
  cat("#!/bin/bash\n")
  cat(glue::glue("#SBATCH --job-name={basename(captures)}.job"), "\n")
  cat(glue::glue("#SBATCH --output={outDirs[i]}/{basename(captures[i])}.out"), "\n")
  cat(glue::glue("#SBATCH --error={outDirs[i]}/{basename(captures[i])}.err"), "\n")
  cat("#SBATCH --time=24:00:00\n")
  cat("#SBATCH --mem=16384\n")
  cat("#\n")
  cat("source activate cellbender\n")
  cat("module load gcc/9.2.0\n")
  cat("#\n")
  cat(glue::glue("cellbender remove-background \
    --input {files[i]} \
    --output {outDirs[i]}/output.h5"), 
      "\n")
  
  # Close the sink!
  sink()
  
  # Submit to run on cluster
  system(glue::glue("sbatch {jobFile}"))
}