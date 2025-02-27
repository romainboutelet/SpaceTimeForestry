# Developing MCMC functions for spatio-temporal markov binary process

### Maybe I can start with just the betas

my_nghb <- function(x, X, D){
  if (is.null(dim(X))) X <- matrix(X, ncol = length(x))
  D[D == 0] <- Inf
  nghb <- which.min(D)
  x0 <- X[nghb,]
  mask <- (X[,1]*(x[1]-x0[1]) > (x[1]^2+x[2]^2-x0[1]^2-x0[2]^2)/2 -
             X[,2]*(x[2]-x0[2]))
  mask[D == Inf] <- F
  count <- 1
  while (any(mask) && count < 5){
    D[!mask] <- Inf
    nghb <- c(nghb, which.min(D))
    x0 <- X[tail(nghb,1),]
    mask <- mask & (X[,1]*(x[1]-x0[1]) > (x[1]^2+x[2]^2-x0[1]^2-x0[2]^2)/2 -
                      X[,2]*(x[2]-x0[2]))
    count <- count + 1
  }
  return(nghb)
}

vecchia_lik_betas <- function(y, x, z, a1, a2, a3, beta0, beta1){
  A <- matrix(c(a1, a2, a2, a3), ncol = 2)
  n <- length(y)
  x <- x[order(x[,1]),]
  y <- y[order(x[,1])]
  z <- z[order(x[,1])]
  dist <- unname(as.matrix(dist(x, diag = T, upper = T)))
  p <- exp(-(beta0+beta1*z))/(1+exp(-(beta0+beta1*z)))
  lik <-  p[1]
  for (i in 2:n){
    D <- dist[i,1:(i-1)]
    nghb <- my_nghb(x[i,], x[1:(i-1),], D)
    h <- sqrt(sum((x[i,]-x[j,])^2))
    p0_tmp <- 1-p[i]
    p1_tmp <- p[i]
    for (j in nghb){
      h <- sqrt(sum((x[i,]-x[j,])^2))
      if (y[j] == 0){
        p1_tmp <- p1_tmp*(1-p[j])*(1-exp(-h/A[2,1]))
        p0_tmp <- p0_tmp*(1-p[j]*(1-exp(-h/A[1,1])))
      }
      if (y[j] == 1){
        p1_tmp <- p1_tmp*(1-(1-p[j])*(1-exp(-h/A[2,2])))
        p0_tmp <- p0_tmp*p[j]*(1-exp(-h/A[1,2]))
      }
    }
    if (y[i] == 0) lik <- lik * p0_tmp / (p0_tmp + p1_tmp)
    else if (y[i] == 1) lik <- lik * p1_tmp / (p0_tmp + p1_tmp)
  }
  return(lik)
}

MCMC_betas <- function(x, y, z, a1, a2, a3, n_samples = 4000, n_burnin = 10,
                       thin = 10, batch = 200){
  n_keep <- n_samples %/% thin
  beta0_out <- rep(NA, n_keep)
  beta0_now <- rnorm(1)
  beta1_out <- rep(NA, n_keep)
  beta1_now <- rnorm(1)
  p_out <- rep(NA, n_samples %/% batch)
  pb1 <- txtProgressBar(min = 0, max = n_burnin, initial = 0, style = 3) 
  for (i in 1:n_burnin){
    beta0_prop <- rnorm(1,beta0_now,0.5)
    beta1_prop <- rnorm(1,beta1_now,0.5)
    r <- log(vecchia_lik_betas(y, x, z, a1, a2, a3, beta0_prop, beta1_prop)) + 
      log(dnorm(beta0_prop,0,10)) + log(dnorm(beta1_prop,0,10)) - (
        log(vecchia_lik_betas(y, x, z, a1, a2, a3, beta0_now, beta1_now)) +
          log(dnorm(beta0_now,0,10)) + log(dnorm(beta1_now,0,10)))
    p <- min(c(1,exp(r)))
    cond <- sample(c(0,1),1,prob = c(1-p, p)) == 1
    if (cond){
      beta1_now <- beta1_prop
      beta0_now <- beta0_prop
    }
    # range parameters sampler
    if ((i %% 100) == 0) setTxtProgressBar(pb1,i)
  }
  count <- 1
  count_batch <- 1
  p_count <- 0
  pb2 <- txtProgressBar(min = 0, max = n_samples, initial = 0, style = 3) 
  for (i in 1:n_samples){
    # beta0 sampler
    beta0_prop <- rnorm(1,beta0_now,0.5)
    beta1_prop <- rnorm(1,beta1_now,0.5)
    r <- log(vecchia_lik_betas(y, x, z, a1, a2, a3, beta0_prop, beta1_prop)) + 
      log(dnorm(beta0_prop,0,10)) + log(dnorm(beta1_prop,0,10)) - (
        log(vecchia_lik_betas(y, x, z, a1, a2, a3, beta0_now, beta1_now)) +
          log(dnorm(beta0_now,0,10)) + log(dnorm(beta1_now,0,10)))
    p <- min(c(1,exp(r)))
    cond <- sample(c(0,1),1,prob = c(1-p, p)) == 1
    if (cond){
      p_count <- p_count + 1
      beta1_now <- beta1_prop
      beta0_now <- beta0_prop
    }
    if ((i %% thin) == 0){
      beta0_out[count] <- beta0_now
      beta1_out[count] <- beta1_now
      count <- count + 1
      setTxtProgressBar(pb2,i)
    }
    if ((i %% batch) == 0){
      p_out[count_batch] <- p_count/batch
      print(paste0("Acceptance rate: ", p_count/batch))
      p_count <- 0
      count_batch <- count_batch + 1
    }
  }
  return(data.frame(beta0 = beta0_out, beta1 = beta1_out, p = p_out))
}

{
  Rprof()
  MCMC_out <- MCMC_betas(xobs, yobs, Z_obs, 3/5, 3/5, 3/5)
  Rprof(NULL)
  print(summaryRprof())
  {
    par(mfrow = c(2,2))
    plot(MCMC_out$beta0, type = "l")
    plot(MCMC_out$beta1, type = "l")
    plot(density(MCMC_out$beta0))
    plot(density(MCMC_out$beta1))
  }
}

MCMC_as <- function(x, y, z, beta0, beta1, n_samples = 20000, n_burnin = 10,
                       thin = 50, batch = 1000){
  n_keep <- n_samples %/% thin
  a_min <- 0.1
  a_max <- 10
  a1_out <- rep(NA, n_keep)
  a1_now <- exp(runif(1,log(a_min),log(a_max)))
  a2_out <- rep(NA, n_keep)
  a2_now <- exp(runif(1,log(a_min),log(a_max)))
  a3_out <- rep(NA, n_keep)
  a3_now <- exp(runif(1,log(a_min),log(a_max)))
  p_out <- rep(NA, n_samples %/% batch)
  pb1 <- txtProgressBar(min = 0, max = n_burnin, initial = 0, style = 3) 
  for (i in 1:n_burnin){
    a1_prop <- exp(rnorm(1,log(a1_now), 0.1))
    a2_prop <- exp(rnorm(1,log(a2_now), 0.1))
    a3_prop <- exp(rnorm(1,log(a3_now), 0.1))
    r <- log(vecchia_lik_betas(y, x, z, a1_prop, a2_prop, a3_prop,
                               beta0, beta1)) +
      log(dunif(a1_prop,a_min,a_max)) + log(dunif(a2_prop,a_min,a_max)) +
      log(dunif(a3_prop,a_min,a_max)) - (
        log(vecchia_lik_betas(y, x, z, a1_now, a2_now, a3_now,
                              beta0, beta1)) +
          log(dunif(a1_now,a_min,a_max)) + log(dunif(a2_now,a_min,a_max)) +
          log(dunif(a3_now,a_min,a_max)))
    p <- min(c(1,exp(r)))
    cond <- sample(c(0,1),1,prob = c(1-p, p)) == 1
    if (cond){
      a1_now <- a1_prop
      a2_now <- a2_prop
      a3_now <- a3_prop
    }
    # range parameters sampler
    if ((i %% 100) == 0) setTxtProgressBar(pb1,i)
  }
  count <- 1
  count_batch <- 1
  p_count <- 0
  pb2 <- txtProgressBar(min = 0, max = n_samples, initial = 0, style = 3) 
  for (i in 1:n_samples){
    # beta0 sampler
    a1_prop <- exp(rnorm(1,log(a1_now), 1))
    a2_prop <- exp(rnorm(1,log(a2_now), 1))
    a3_prop <- exp(rnorm(1,log(a3_now), 1))
    r <- log(vecchia_lik_betas(y, x, z, a1_prop, a2_prop, a3_prop,
                               beta0, beta1)) + 
      log(dunif(a1_prop,a_min,a_max)) + log(dunif(a2_prop,a_min,a_max)) +
      log(dunif(a3_prop,a_min,a_max)) - (
        log(vecchia_lik_betas(y, x, z, a1_now, a2_now, a3_now,
                              beta0, beta1)) + 
          log(dunif(a1_now,a_min,a_max)) + log(dunif(a2_now,a_min,a_max)) +
          log(dunif(a3_now,a_min,a_max)))
    p <- min(c(1,exp(r)))
    cond <- sample(c(0,1),1,prob = c(1-p, p)) == 1
    if (cond){
      p_count <- p_count + 1
      a1_now <- a1_prop
      a2_now <- a2_prop
      a3_now <- a3_prop
    }
    if ((i %% thin) == 0){
      a1_out[count] <- a1_now
      a2_out[count] <- a2_now
      a3_out[count] <- a3_now
      count <- count + 1
      setTxtProgressBar(pb2,i)
    }
    if ((i %% batch) == 0){
      p_out[count_batch] <- p_count/batch
      print(paste0("Acceptance rate: ", p_count/batch))
      p_count <- 0
      count_batch <- count_batch + 1
    }
  }
  return(data.frame(a1 = a1_out,a2 = a2_out, a3 = a3_out, p = p_out))
}

{
  Rprof()
  MCMC_out <- MCMC_as(xobs, yobs, Z_obs, beta0, beta1)
  Rprof(NULL)
  print(summaryRprof())
  {
    par(mfrow = c(2,3))
    plot(log(MCMC_out$a1), type = "l")
    plot(log(MCMC_out$a2), type = "l")
    plot(log(MCMC_out$a3), type = "l")
    plot(density(MCMC_out$a1, from = 0))
    plot(density(MCMC_out$a2, from = 0))
    plot(density(MCMC_out$a3, from = 0))
  }
}

MCMC_all <- function(x, y, z, n_samples = 5000, n_burnin = 1000,
                    thin = 10, batch = 200){
  n_keep <- n_samples %/% thin
  beta0_out <- rep(NA, n_keep)
  beta0_now <- rnorm(1)
  beta1_out <- rep(NA, n_keep)
  beta1_now <- rnorm(1)
  a_min <- 0.1
  a_max <- 10
  a1_out <- rep(NA, n_keep)
  a1_now <- exp(runif(1,log(a_min),log(a_max)))
  a2_out <- rep(NA, n_keep)
  a2_now <- exp(runif(1,log(a_min),log(a_max)))
  a3_out <- rep(NA, n_keep)
  a3_now <- exp(runif(1,log(a_min),log(a_max)))
  p_out <- rep(NA, n_samples %/% batch)
  print("Burn-in phase")
  pb1 <- txtProgressBar(min = 0, max = n_burnin, initial = 0, style = 3) 
  for (i in 1:n_burnin){
    beta0_prop <- rnorm(1,beta0_now,1)
    beta1_prop <- rnorm(1,beta1_now,1)
    a1_prop <- exp(rnorm(1,log(a1_now),1))
    a2_prop <- exp(rnorm(1,log(a2_now),1))
    a3_prop <- exp(rnorm(1,log(a3_now),1))
    r <- log(vecchia_lik_betas(y, x, z, a1_prop, a2_prop, a3_prop,
                               beta0_prop, beta1_prop)) +
      log(dnorm(beta0_prop,0,10)) + log(dnorm(beta1_prop,0,10)) +
      log(dunif(a1_prop,a_min,a_max)) + log(dunif(a2_prop,a_min,a_max)) +
      log(dunif(a3_prop,a_min,a_max)) - (
        log(vecchia_lik_betas(y, x, z, a1_now, a2_now, a3_now,
                              beta0_now, beta1_now)) +
          log(dnorm(beta0_now,0,10)) + log(dnorm(beta1_now,0,10)) +
          log(dunif(a1_now,a_min,a_max)) + log(dunif(a2_now,a_min,a_max)) +
          log(dunif(a3_now,a_min,a_max)))
    p <- min(c(1,exp(r)))
    cond <- sample(c(0,1),1,prob = c(1-p, p)) == 1
    if (cond){
      beta1_now <- beta1_prop
      beta0_now <- beta0_prop
      a1_now <- a1_prop
      a2_now <- a2_prop
      a3_now <- a3_prop
    }
    # range parameters sampler
    if ((i %% 100) == 0) setTxtProgressBar(pb1,i)
  }
  print("Sampled phase")
  count <- 1
  count_batch <- 1
  p_count <- 0
  pb2 <- txtProgressBar(min = 0, max = n_samples, initial = 0, style = 3) 
  for (i in 1:n_samples){
    beta0_prop <- rnorm(1,beta0_now,0.25)
    beta1_prop <- rnorm(1,beta1_now,0.25)
    a1_prop <- exp(rnorm(1,log(a1_now),0.5))
    a2_prop <- exp(rnorm(1,log(a2_now),0.5))
    a3_prop <- exp(rnorm(1,log(a3_now),0.5))
    r <- log(vecchia_lik_betas(y, x, z, a1_prop, a2_prop, a3_prop,
                               beta0_prop, beta1_prop)) + 
      log(dnorm(beta0_prop,0,10)) + log(dnorm(beta1_prop,0,10)) +
      log(dunif(a1_prop,a_min,a_max)) + log(dunif(a2_prop,a_min,a_max)) +
      log(dunif(a3_prop,a_min,a_max)) - (
        log(vecchia_lik_betas(y, x, z, a1_now, a2_now, a3_now,
                              beta0_now, beta1_now)) + 
          log(dnorm(beta0_now,0,10)) + log(dnorm(beta1_now,0,10)) +
          log(dunif(a1_now,a_min,a_max)) + log(dunif(a2_now,a_min,a_max)) +
          log(dunif(a3_now,a_min,a_max)))
    p <- min(c(1,exp(r)))
    cond <- sample(c(0,1),1,prob = c(1-p, p)) == 1
    if (cond){
      p_count <- p_count + 1
      beta1_now <- beta1_prop
      beta0_now <- beta0_prop
      a1_now <- a1_prop
      a2_now <- a2_prop
      a3_now <- a3_prop
    }
    if ((i %% thin) == 0){
      beta0_out[count] <- beta0_now
      beta1_out[count] <- beta1_now
      a1_out[count] <- a1_now
      a2_out[count] <- a2_now
      a3_out[count] <- a3_now
      count <- count + 1
      setTxtProgressBar(pb2,i)
    }
    if ((i %% batch) == 0){
      p_out[count_batch] <- p_count/batch
      print(paste0("Acceptance rate: ", p_count/batch))
      p_count <- 0
      count_batch <- count_batch + 1
    }
  }
  return(data.frame(beta0 = beta0_out, beta1 = beta1_out, a1 = a1_out,
                    a2 = a2_out, a3 = a3_out, p = p_out))
}

{
  Rprof()
  MCMC_out <- MCMC_all(xobs, yobs, Z_obs)
  Rprof(NULL)
  print(summaryRprof())
}
{
  par(mfrow = c(2,2))
  plot(MCMC_out$beta0, type = "l")
  plot(MCMC_out$beta1, type = "l")
  plot(density(MCMC_out$beta0))
  plot(density(MCMC_out$beta1))
}
{
  par(mfrow = c(2,3))
  plot(log(MCMC_out$a1), type = "l")
  plot(log(MCMC_out$a2), type = "l")
  plot(log(MCMC_out$a3), type = "l")
  plot(density(MCMC_out$a1, from = 0))
  plot(density(MCMC_out$a2, from = 0))
  plot(density(MCMC_out$a3, from = 0))
}

MCMC_pred <- function(MCMC_samples, X, x, y, z, dist){
  n_samples <- nrow(MCMC_samples)
  p_out <- z*0
  pb <- txtProgressBar(min = 0, max = n_samples, initial = 0, style = 3) 
  for (s in 1:n_samples){
    a1 <- MCMC_samples$a1[s]
    a2 <- MCMC_samples$a2[s]
    a3 <- MCMC_samples$a3[s]
    A <- matrix(c(a1, a2, a2, a3), nrow = 2, ncol = 2)
    b0 <- MCMC_samples$beta0[s]
    b1 <- MCMC_samples$beta1[s]
    p <- exp(-(b0+b1*z))/(1+exp(-(b0+b1*z)))
    Pcond <- p
    for (i in 1:npred){
      D <- dist[i,]
      nghb <- my_nghb(X[i,], x, D)
      tmp <- 1-Pcond[i]
      j <- nghb[1]
      h <- sqrt(sum((X[i,]-x[j,])^2))
      for (j in nghb){
        h <- dist[i,j]
        if (y[j] == 0){
          Pcond[i] <- Pcond[i]*(1-p[j])*(1-exp(-h/A[2,1]))
          tmp <- tmp*(1-p[j]*(1-exp(-h/A[1,1])))
        }
        if (y[j] == 1){
          Pcond[i] <- Pcond[i]*(1-(1-p[j])*(1-exp(-h/A[2,2])))
          tmp <- tmp*p[j]*(1-exp(-h/A[1,2]))
        }
      }
      Pcond[i] <- Pcond[i]/(Pcond[i]+tmp)
    }
    p_out <- p_out + Pcond/n_samples
    if ((s %% 10) == 0) setTxtProgressBar(pb, s)
  }
  return(p_out)
}

