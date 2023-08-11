###
# 03-lit-time.R calculates computational time as a function of sample size
# and number of traits in the simulation study
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
  dat <- simulation_study(N = as.integer(N),
                          num_traits = num_traits,
                          num_null = num_null,
                          var_part = var_part,
                          counteract = counteract,
                          positive = positive,
                          error = error)

  dat$X <- dat$X[1,]
  X <- as.matrix(dat$X)
  Y <- as.matrix(dat$Y)
  t1 <- proc.time()[3]
  p <- lit::lit(y = Y, x = X)
  t2 <- proc.time()[3] - t1

  t3 <- proc.time()[3]
  p_marginal <- lit::marginal(y = dat$Y,
                              x = as.matrix(dat$X))
  t4 <- proc.time()[3] - t3

  df <- data.frame(method = c("aLIT", "Marginal"),
                   time = c(as.numeric(t2), as.numeric(t4)))
  rm(dat)
  as_tibble(df)
}

library(tidyverse)
library(digest)
library(parallel)
library(lit)
design <- expand.grid(N = seq(0.1, 5, .7) * 10 ^ 5,
                      id = 1:500,
                      heritability = 0.1,
                      env = 0.15,
                      counteract = FALSE,
                      num_traits = c(5, 10),
                      prop_null = 0,
                      positive = TRUE,
                      error = "gaussian")

design <- design %>%
  group_by(N, id, heritability, env, counteract, num_traits, prop_null, positive, error) %>%
  mutate(seed = readBin(digest(c(N, id, heritability, env, counteract, num_traits, prop_null, positive, error), raw=TRUE), "integer"))

t1 <- proc.time()
cl <- makeCluster(1, type = "PSOCK") #"PSOCK
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
t2 <- proc.time() - t1

out <- dplyr::bind_rows(out, .id = "column_label")
design$column_label <- as.character(1:nrow(design))
out <- design %>% right_join(out)

saveRDS(out, file = "./data/03-lit-time.rds")
