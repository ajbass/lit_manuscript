#!/bin/bash
#$ -cwd
#$ -N plink
#$ -o 02-plink.o
#$ -e 02-plink.e
#$ -m e

module load plink/2.0a

for chr in {1..22}; do \
  plink --bgen /mnt/EpsteinFSS/data/UKBiobank/genotype_data/imputed_genotypes/ukb22828_c${chr}_b0_v3.bgen --sample /mnt/EpsteinFSS/data/UKBiobank/genotype_data/imputed_genotypes/ukb22828_c${chr}_b0_v3.sample --keep ./data/obesity_individuals.tab  --hwe 1e-5 --maf 0.05 --geno .05 --make-bed --out ./plink/chr_${chr}; \
done
