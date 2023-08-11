###
# 02-lit-alt.R compares the different LIT implementations
# discussed in the manuscript.
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
  # Function runs simulation study
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
  p_lit <- lit::lit(y = dat$Y,
                    x = as.matrix(dat$X))

  p_marginal <- lit::marginal(y = dat$Y,
                              x = as.matrix(dat$X))

  # estimate the number of independent tests for SQ + CP
  tmp <- scale(dat$Y)
  oo <- cbind(tmp ^ 2, lit:::.pairwise_prod(tmp))
  oo <- cor(oo)
  d <- eigen(oo, symmetric = TRUE, only.values = TRUE)$values
  tot_tests <- sum(cumsum(d/sum(d)) < 0.95)
  marginal_SQCP <- min(min(as.numeric(p_marginal)) * tot_tests, 1)

  # same for SQ
  oo <- tmp ^ 2
  oo <- cor(oo)
  d <- eigen(oo, symmetric = TRUE, only.values = TRUE)$values
  tot_tests2 <- sum(cumsum(d/sum(d)) < 0.95)
  tt = choose(num_traits, 2)
  marginal_SQ <- min(min(as.numeric(p_marginal[(tt+1):(tt+num_traits)])) * tot_tests2, 1)

  df <- data.frame(method =  c("wLIT", "uLIT", "aLIT", "Marginal", "Marginal_SQ"),
                   p.value = c(as.numeric(p_lit), marginal_SQCP, marginal_SQ),
                   Neff = c(0, 0, 0, tot_tests, tot_tests2))

  rm(dat)
  as_tibble(df)
}

library(tidyverse)
library(digest)
library(parallel)
library(lit)
# Alternative simulations for 5 traits
design <- expand.grid(N = c(300000),
                      id = 1:500,
                      heritability = c(0.1, 0.35, 0.6),
                      env = c(0.15),
                      counteract = c(TRUE, FALSE),
                      num_traits = 5,
                      prop_null =  0:4 / 5,
                      positive = c(TRUE, FALSE),
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

saveRDS(out, file = "./data/02-lit-alt-5-traits.rds")

# Alternative simulations for 10 traits
design <- expand.grid(N = c(300000),
                      id = 1:500,
                      heritability = c(0.1, 0.35, 0.6),
                      env = c(0.15),
                      counteract = c(TRUE, FALSE),
                      num_traits = 10,
                      prop_null =  0:9 / 10,
                      positive = c(TRUE, FALSE),
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

saveRDS(out, file = "./data/02-lit-alt-10-traits.rds")
