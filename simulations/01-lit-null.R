###
# 01-lit-null.R compares the different LIT implementations
# under the null hypothesis of no latent interactions.
###
source("./00-functions.R")
run_study <- function(N,
                      num_traits,
                      num_null,
                      var_part,
                      positive,
                      B = 10000,
                      counteract,
                      error,
                      seed) {
  # Function runs simulation study
  # N: sample size
  # num_traits: number of traits
  # num_null: number of null traits
  # var_part: variance partition of components
  # positive: positive pleiotropy
  # B: number of null iterations
  # counteract: effect sign of interactive environment opposite of interaction
  # error: distribution of error r.v.: t, Chi-squared, and Gaussian
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

  p <- matrix(nrow = 3, ncol = B)
  for (i in 1:B) {
    p0 <- runif(n = 1, 0.1, 0.4)
    X <- matrix(rbinom(N,
                       size = 2,
                       prob = p0),
                nrow = N,
                ncol = 1)
    p[, i] <- as.numeric(lit(y = as.matrix(dat$Y), x = X))
  }

  df <- data.frame(p.value = as.numeric(p),
                   method =  rep(c("lit_ev", "lit_ev_eq", "lit"), B),
                   N = N,
                   error = error,
                   positive = positive,
                   correlation = sum(var_part[1:2]),
                   seed = seed,
                   num_traits = num_traits,
                   B = B)

  saveRDS(df, file = paste0("./data/null/01-lit-null-", seed, ".rds"))
  return(NULL)
}

library(tidyverse)
library(digest)
library(parallel)
library(lit)
design <- expand.grid(N = c(300000),
                      id = 1:100,
                      heritability = c(0.10, 0.35, 0.6),
                      env = 0.15,
                      counteract = TRUE,
                      num_traits = c(5, 10),
                      prop_null =  1,
                      positive = TRUE,
                      error = c("gaussian", "chisq", "t"))

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
