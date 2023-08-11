###############################################
## Marginal (SQ/CP) and Marginal(SQ) applied ##
## to UK Biobank                             ##
###############################################
#setwd("~/Dropbox/postdoc/ajbass/lit/analysis/ukb")

# Load packages + data
library(parallel)
library(tidyverse)
library(genio)
library(lit)
path <- "/Users/ajbass/UKB_data/imputed_genotypes/chunks/c22/chr_22_p1"
adjustment <- read_tsv("./data/obesity_covariates.tab")
phenotype <- read_tsv("./data/obesity_phenotypes.tab")
fam_file <- genio::read_fam(file = path)

# Make sure order in file is the same as bed file (should be)
tmp <- data.frame(IID = as.numeric(fam_file$id))
phenotype <- tmp %>% left_join(phenotype)
adjustment <- tmp %>% left_join(adjustment)
adjustment <- as.matrix(adjustment[, 3:24])
phenotype <- as.matrix(phenotype[, -c(1)])
PC <- svd(as.matrix(adjustment[, 3:22]))$u
rm(tmp, fam_file)

run_marginal <- function(file) {
  library(tidyverse)
  library(lit)
  out <- marginal_plink(y = phenotype,
                        file = file,
                        pop_struct = PC,
                        verbose = FALSE)

  chr = str_split(file, pattern = "/")[[1]][9]
  save(out,
       file = paste0("./data/ukb_marginal/", chr, ".RData"))
  rm(out)
  return("completed")
}

## Run cluster
fnames <- list.files("/Users/ajbass/UKB_data/imputed_genotypes/chunks/", recursive = TRUE, full.names = TRUE)
fnames <- unique(sapply(fnames,FUN = function(i) str_split(i, "\\.")[[1]][1]))

cl <- makeCluster(12, type = "PSOCK")
clusterExport(cl, varlist = c("run_marginal", "phenotype", "PC", "fnames"))
out <- parLapply(cl, 1:length(fnames), function(i) {
  return(run_marginal(file = fnames[i]))
})
stopCluster(cl)
