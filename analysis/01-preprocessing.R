##################################################
### File preprocesses UK Biobank data  ###
## Interested in obesity-related traits for
## individuals (unrelated) with British ancestry
##################################################
library(tidyverse)

# Sample information
sample_info <- read_delim("./data/ukb_sqc_v2.txt",
                          delim = " ",
                          col_names = c("x1",
                                        "x2",
                                        "genotyping.array",
                                        "Batch",
                                        "Plate.Name",
                                        "Well",
                                        "Cluster.CR",
                                        "dQC",
                                        "Internal.Pico..ng.uL.",
                                        "Submitted.Gender",
                                        "Inferred.Gender",
                                        "X.intensity",
                                        "Y.intensity",
                                        "Submitted.Plate.Name",
                                        "Submitted.Well",
                                        "sample.qc.missing.rate",
                                        "heterozygosity",
                                        "heterozygosity.pc.corrected",
                                        "het.missing.outliers",
                                        "putative.sex.chromosome.aneuploidy",
                                        "in.kinship.table",
                                        "excluded.from.kinship.inference",
                                        "excess.relatives",
                                        "in.white.British.ancestry.subset",
                                        "used.in.pca.calculation",
                                        paste("pc", 1:40, sep = ""),
                                        "in.Phasing.Input.chr1_22",
                                        "in.Phasing.Input.chrX",
                                        "in.Phasing.Input.chrXY"
                          ))

# Focus on those used in PC calculation as these represent unrelated individuals
sample_info_onlyPCA <- sample_info %>%
  filter(used.in.pca.calculation == 1,
         in.white.British.ancestry.subset == 1,
         putative.sex.chromosome.aneuploidy == 0,
         excess.relatives==0,
         in.Phasing.Input.chr1_22==1)

## PCs
PC <- sample_info_onlyPCA %>%
  select(Batch, Submitted.Gender, pc1:pc40)
write_csv(PC, file = "./data/PCs.txt", col_names = TRUE)

## Sample IDs
fam_info <- read_delim("./data/ukb22418_c2_b0_v2.fam",
                       delim = " ", col_names = FALSE)
id <- fam_info[sample_info$used.in.pca.calculation == 1 &
               sample_info$in.white.British.ancestry.subset == 1 &
               sample_info$putative.sex.chromosome.aneuploidy == 0 &
               sample_info$excess.relatives==0 &
               sample_info$in.Phasing.Input.chr1_22==1, 1]
write_csv(id, file = "./data/samples_in_PCA.txt", col_names = FALSE)

## Load UK Biobank phenotypes
df = read_delim("./data/UKB_phenotypes.tab", delim = "\t",
                col_select = c(1, 4266, 4254, 28, 32, 4436))
colnames(df) <- c("IID", "age", "BMI", "WC", "HC", "BFP")

# Remove individuals that do not have measured values across all traits
id <- apply(df, 1, anyNA)
df <- df[!id,]

# British ancestry + unrelated individuals
samples_PCA <- read_delim("./data/samples_in_PCA.txt", delim = "\t", col_names = "IID")
df <- df %>% filter(IID %in% samples_PCA$IID)

## Filter PCs with people that have phenotype measurements
pcs <- read_delim("./data/PCs.txt", delim = ",")
id_samples <- samples_PCA$IID %in% df$IID
pcs <- pcs[id_samples,]

## Organize for plink
sex <- model.matrix(~as.factor(pcs$Submitted.Gender))[,-1]
out <- data.frame(FID = samples_PCA$IID[id_samples],
                  IID = samples_PCA$IID[id_samples],
                  sex = sex,
                  pcs[, -(1:2)])
out <- df %>%
  select(IID, age) %>%
  left_join(out)
out <- out[, c("IID", "FID", "age", "sex", colnames(out)[-(1:4)])]
out_save <- out
df$sex <- out$sex

## Withdrawals
withdrawals <- read_tsv("./data/withdrawals.csv", col_names = FALSE)
withdrawals <- withdrawals %>% rename(IID = X1)
df <- df %>% filter(!(IID %in% withdrawals$IID))
out_save <- out_save %>% filter(IID %in% df$IID, !(IID %in% withdrawals$IID))

## Remove top 20 PCs and age
df$WC <- residuals(lm(df$WC ~ df$age + as.matrix(out_save[, 5:24])))
df$HC <- residuals(lm(df$HC ~ df$age + as.matrix(out_save[, 5:24])))
df$BMI <- residuals(lm(df$BMI ~ df$age + as.matrix(out_save[, 5:24])))
df$BFP <- residuals(lm(df$BFP ~ df$age + as.matrix(out_save[, 5:24])))

# Remove outliers
df <- df %>%
  group_by(sex) %>%
  filter(abs(scale(WC)) < 4 &
           abs(scale(HC)) < 4 &
           abs(scale(BMI)) < 4 &
           abs(scale(BFP) < 4))

# Scale by sex
df <- df %>% group_by(sex) %>%
  mutate(WC = as.numeric(scale(WC)),
         HC = as.numeric(scale(HC)),
         BMI = as.numeric(scale(BMI)),
         BFP = as.numeric(scale(BFP)))

out_save <- out_save %>% filter(IID %in% df$IID)

## Save files
write_tsv(out_save[, 1:2], file = "./data/obesity_individuals.tab")
write_tsv(out_save, file = "./data/obesity_covariates.tab")
write_tsv(df[, c("IID", "BMI", "WC", "HC", "BFP")],
          file = "./data/obesity_phenotypes.tab")
