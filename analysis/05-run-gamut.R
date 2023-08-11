####################################
## GAMuT applied to UK Biobank #####
####################################
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

run_gamut <- function(file) {
  library(tidyverse)
  library(lit)
  out <- gamut_plink(y = phenotype,
                     file = file,
                     pop_struct = PC,
                     verbose = FALSE)

  chr = str_split(file, pattern = "/")[[1]][4]
  save(out,
       file = paste0("./data/ukb_gamut/gamut_", chr, ".RData"))
  rm(out)
  return("completed")
}

# Run cluster on LD regions (estimated from PLINK only at lead SNPs)
chr_LD <- c(1, 2, 3, 5, 6, 7, 11, 12, 16, 18)
cl <- makeCluster(10, type = "PSOCK")
clusterExport(cl, varlist =  c("run_gamut", "chr_LD", "phenotype", "PC", "fnames"))
out <- parLapply(cl, 1:length(fnames), function(i) {
  return(run_gamut(file = paste0("./data/LD/LD_", chr_LD[i])))
})
stopCluster(cl)

# combine all files
dat <- NULL
for (i in chr_LD) {
  load(paste0("./data/gamut_LD_", i, ".RData"))
  dat <- rbind(out, dat)
}
saveRDS(dat, "./data/gamut.rds")
