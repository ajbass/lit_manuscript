###
# 04-PC.R performs an association test between the genotype and eigenvectors of
# the SQ/CP matrix
###
source("./00-functions.R")
run_study <- function(N,
                      num_traits,
                      num_null,
                      var_part,
                      positive,
                      counteract,
                      error,
                      seed) {
  # Function run simulation study
  # N: sample size
  # num_traits: number of traits
  # num_null: number of null traits
  # var_part: variance partition of components
  # positive: positive pleiotropy
  # counteract (boolean): interactive env. effect size opposite of interaction
  # error: distribution of error term (gaussian)
  # seed: seed for simulations
  set.seed(seed)
  dat <- simulation_study(N = N,
                          num_traits = num_traits,
                          num_null = num_null,
                          var_part = var_part,
                          counteract = counteract,
                          positive = positive,
                          error = error)
  dat$X <- dat$X[1,]

  tmp <- residuals(lm(dat$Y ~ dat$X))
  oo <- cbind(tmp ^ 2, lit:::.pairwise_prod(tmp))
  svd.out <- RSpectra::svds(scale(oo), k = 44)
  p <- sapply(1:44, FUN= function(x) broom::tidy(anova(lm(svd.out$u[,x] ~ dat$X)))$p.value[1])

  df <- data.frame(PC = 1:44, p = p)
  rm(dat)
  return(df)
}

library(tidyverse)
library(digest)
library(parallel)
library(lit)
# Alternative simulations for 5 traits
design <- expand.grid(N = c(300000),
                      id = 1:500,
                      heritability = seq(0, 0.95, 0.05),
                      env = c(0.0),
                      counteract = c(FALSE),
                      num_traits = 10,
                      prop_null =  0:9 / 10,
                      positive = c(TRUE),
                      error = "gaussian")

design <- design %>%
  group_by(N, id, heritability, env, counteract, num_traits, prop_null, positive, error) %>%
  mutate(seed = readBin(digest(c(N, id, heritability, env, counteract, num_traits, prop_null, positive, error), raw=TRUE), "integer"))

cl <- makeCluster(25, type = "PSOCK")
clusterExport(cl, varlist = c("design", "lit", "run_study", "simulation_study", "phenotype", "generate_parameters", "generate_predictors"))
out <- parLapply(cl, 1:nrow(design), function(ii) {
  library(tidyverse)
  return(run_study(design[ii,]$N,
                   var_part = c(design[ii,]$heritability, design[ii,]$env),
                   num_traits = design[ii,]$num_traits,
                   num_null = round(design[ii,]$num_traits * design[ii,]$prop_null),
                   positive = design[ii,]$positive,
                   counteract = design[ii,]$counteract,
                   error = design[ii,]$error,
                   seed = design[ii,]$seed))
})
stopCluster(cl)

out <- dplyr::bind_rows(out, .id = "column_label")
design$column_label <- as.character(1:nrow(design))
out <- design %>%
  right_join(out)

saveRDS(out, file = "./data/04-PC.rds")
