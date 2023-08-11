generate_predictors <- function(out, o_genetic = 99) {
  # Function generates components of the polygenic trait model
  # out: list of input parameters
  # o_genetic: the total number of SNPs that contribute to additive signal

  # Initialization and population structure
  N <- out$N[1]
  p <- 0.25

  # Generate genotypes for SNP involved in interaction
  X <- t(matrix(rbinom(N, size = 2, prob = p),
                nrow = N,
                ncol = 1))

  # Generate genotypes for rest of additive signal
  p1 <- runif(n = o_genetic, 0.1, 0.4)
  X_other <- t(matrix(rbinom(o_genetic * N, size = 2, prob = t(p1)),
                      nrow = N, ncol = o_genetic))

  # Environmental variable involved in interaction
  W <- matrix(rnorm(N), ncol=1)

  list(X = X,
       X_other = X_other,
       W = as.vector(W),
       W_other = rnorm(N), # environmental variable shared across traits
       XW = as.vector(as.numeric(X) * as.numeric(W)))
}

generate_parameters <- function(N, num_traits, positive, counteract, o_genetic = 99) {
  # Function generates effect sizes
  # N: sample size
  # num_traits: number of traits
  # positive: effect size direction across traits
  # counteract: effect size of environment opposite of interaction
  # o_genetic: total number of SNPs that contribute to additive signal

  # Generate effect sizes
  es_X <- rnorm(num_traits, sd = 0.01)
  es_oX <- matrix(rep(rnorm(1 * o_genetic, sd = 0.01), each = num_traits),
                  nrow = num_traits,
                  ncol = o_genetic)
  es_XW <- rnorm(num_traits, sd = 0.01)
  es_W <- rnorm(num_traits, sd = 0.01)

  # Assign directions of effect sizes
  if (positive) {
      es_XW <- abs(es_XW)
      if (counteract) {
        es_W <- -1 * abs(es_W)
      } else {
        es_W <- 1 * abs(es_W)
      }
  } else {
    sign <- sample(c(1, -1), size = num_traits, replace = T)
    es_XW <- abs(es_XW)
    if (counteract) {
      es_W <- -1 * abs(es_W)
    } else {
      es_W <- 1 * abs(es_W)
    }
    es_W <- sign * es_W
    es_XW <- sign * es_XW
    }

  es_X <- sign(es_XW) * abs(es_X) # same sign as interaction
  intercept <- rnorm(num_traits, sd = 5)
  list(N = N,
       intercept = intercept,
       es_X = es_X,
       es_W = es_W,
       es_XW = es_XW,
       es_oX = es_oX)
}

phenotype <- function(out, pred, var_part, id,  null, error = "gaussian") {
  # Function generates a trait
  # out: list of parameters
  # pred: list of predictos
  # var_part: variance partition of each component
  # id: indicator for the trait
  # null: indicator of whether interaction term exists for trait
  # error: distribution of error

  # Extract variables and coefficients for traits
  es_XW <- out$es_XW[id]; es_X <- out$es_X[id];es_W <- out$es_W[id]; es_oX <- out$es_oX[id,,drop = F];
  N <- out$N[1]; intercept <- out$intercept[id]; es_Z <- out$es_Z[id];
  X <- pred$X; W <- pred$W; W_other <- pred$W_other; X_other <- pred$X_other;
  XW <- pred$XW;
  X <- as.vector(X)

  # Generate signal components
  signal_prim_genetic <- X * es_X
  signal_other_genetic <- as.numeric(es_oX %*% X_other)
  signal_interaction <- XW * es_XW
  W <- es_W * W
  if (error == "gaussian") {
    error <- rnorm(N)
  } else if (error == "chisq") {
    error <- rchisq(N, df = 3)
  } else{
    error <- rt(N, df = 5)
  }

  # Standardize components
  signal_prim_genetic <- (signal_prim_genetic - mean(signal_prim_genetic)) / sd(signal_prim_genetic)
  signal_other_genetic <- (signal_other_genetic - mean(signal_other_genetic)) / sd(signal_other_genetic)
  signal_interaction <- (signal_interaction - mean(signal_interaction)) / sd(signal_interaction)
  env <- (W - mean(W)) / sd(W)
  env_other <- (W_other - mean(W_other)) / sd(W_other)
  error <- matrix((error - mean(error)) / sd(error), ncol = 1)

  # Polygenic trait model
  additive <- 0.002 # PVE of additive component
  env_pve <- runif(n = 1, 0.005, 0.02)#runif(n = 1, 0.01, 0.025) #1-2.5 0.005-0.02
  int_pve <- runif(n = 1, 0.001, 0.0015)#runif(n = 1, 0.0005, 0.0015) #0.005-0.0015 ... 0.001-0.002
  if (!null) {
    y <- intercept +
      sqrt(additive) * signal_prim_genetic +
      sqrt(var_part[1]) * signal_other_genetic +
      sqrt(env_pve) * env +
      sqrt(var_part[2]) * env_other +
      sqrt(int_pve) * signal_interaction +
      sqrt(1 - var_part[1] - var_part[2] - int_pve - env_pve - additive) * error
  } else {
    y <- intercept +
      sqrt(additive) * signal_prim_genetic +
      sqrt(var_part[1]) * signal_other_genetic +
      sqrt(env_pve) * env +
      sqrt(var_part[2]) * env_other +
      sqrt(1 - var_part[1] - var_part[2] - env_pve - additive) * error
  y
  }
}

simulation_study <- function(N,
                             num_traits,
                             num_null,
                             var_part,
                             counteract,
                             positive,
                             error) {
  # Function generates a set of traits
  # N: sample size
  # num_traits: number of traits
  # num_null: number of null traits
  # var_part: variance partition of the components
  # counteract: effect size direction of interactive env relative to interaction
  # positive: Indicator for positive pleiotropy
  # error: distribution of error term (gaussian, t, chisq)

  # Generate coefficients/predictors
  coef <- generate_parameters(N = N,
                              num_traits = num_traits,
                              positive = positive,
                              counteract = counteract)
  pred <- generate_predictors(coef)

  Y <- matrix(nrow = N, ncol = num_traits)
  Y <- matrix(rnorm(N * num_traits), ncol = num_traits)
  if (num_traits != num_null) {
    for (i in 1:(num_traits - num_null)) {
      Y[, i] <- phenotype(coef, pred, var_part = var_part,
                          null = FALSE, error = error, id = i)
    }
  }
  if (num_null != 0) {
    for (i in (num_traits-num_null+1):num_traits) {
      # Generate predictors
      Y[, i] <- phenotype(coef, pred, var_part = var_part,
                          null = TRUE, error = error, id = i)
    }
  }
  list(Y = Y,
       X = pred$X,
       W = pred$W)
}
