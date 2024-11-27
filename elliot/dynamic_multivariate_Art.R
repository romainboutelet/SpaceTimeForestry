library(sf)
library(dplyr)
library(spdep)
library(ggplot2)
library(MCMCpack)
library(coda)
library(patchwork)
library(wesanderson)
library(pals)

set.seed(123)

rmvn <- function(n, mu=0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p))))
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
}

rinvgamma <- function(n, shape, scale){
  1/rgamma(n, shape = shape, rate = scale)
}

logit <- function(theta, a, b){log((theta-a)/(b-theta))}

logit.inv <- function(z, a, b){b-(b-a)/(1+exp(z))}

inv <- function(X) {
  chol2inv(chol(X))
}

##########################
## Data.
##########################

dat <- readRDS("../../data_prep/plot_level/output/plot_level_carbon_pools_conus.rds")

shp <- dat$shp
n_jt <- dat$n_jt

y <- dat$y
X <- dat$X
X_tilde <- dat$X_tilde

years <- dat$years

T <- dat$T
J <- dat$J
M <- dat$M
P <- dat$P
Q <- dat$Q
R <- dat$R

W <- dat$W
D <- dat$D
d_j <- diag(D)
N_j <- poly2nb(shp)
N_j[[2494]] <- which(W[2494,] == 1) # Nantucket
N_j[[3095]] <- which(W[3095,] == 1) # San Juan

# Fast CAR stuff
D.sqrt.inv <- diag(1/sqrt(diag(D)))
eig <- eigen(D.sqrt.inv %*% W %*% D.sqrt.inv)
DP <- sqrt(D)%*%eig$vectors

S.1 <- matrix(0, nrow = J, ncol = J)
S.2 <- matrix(0, nrow = J, ncol = J)
for(j in 1:J){
  S.1 <- S.1 + (DP[,j] %*% t(DP[,j]))
  S.2 <- S.2 + eig$values[j] * (DP[,j] %*% t(DP[,j]))
}

regions <- dat$regions
r_j <- match(st_drop_geometry(shp$eco), regions)
J_r <- split(seq_along(r_j), r_j)

##########################
## Priors.
##########################

mu_beta <- rep(0, times = M*P)
Sigma_beta <- diag(1000, nrow = M*P)

a_eta <- 0
b_eta <- 1
c_eta <- 2
d_eta <- 100

a_v <- 0
b_v <- 1

H_delta <- diag(1, nrow = M, ncol = M)
nu_delta <- M+1

H_u <- diag(1, nrow = M, ncol = M)
nu_u <- M+1

H_epsilon <- diag(1, nrow = M, ncol = M)
nu_epsilon <- M+1

##########################
## Starting.
##########################

y.s <- y
mu.s <- array(0, dim = c(M, J, T))
beta.s <- rep(0, times = M*P)
u.s <- array(0, dim = c(J, M, T))
eta.s <- array(0, dim = c(M, Q, J))
Sigma_epsilon.s <- array(NA, dim = c(M, M, R))
for (r in 1:R) {
  Sigma_epsilon.s[,,r] <- diag(100, nrow = M, ncol = M)
}
Sigma_u.s <- array(NA, dim = c(M, M, R, T))
for (t in 1:T) {
  for (r in 1:R) {
    Sigma_u.s[,,r,t] <- diag(100, nrow = M, ncol = M)
  }
}
Sigma_delta.s <- array(NA, dim = c(M, M, J))
for (j in 1:J) {
  Sigma_delta.s[,,j] <- diag(100, nrow = M, ncol = M)
}
tau_sq_eta.s <- matrix(100, nrow = M, ncol = Q)
rho_v.s <- rep(0.5, times = M)
rho_eta.s <- matrix(0.5, nrow = M, ncol = Q)

##########################
## Sampler Prep.
##########################

n.iter <- 500
keep <- seq(250, n.iter, by=2)
n.keep <- length(keep)

rho_v_tuning <- 0.01
rho_eta_tuning <- 0.01

##########################
## Storage.
##########################

mu.samples <- array(NA, dim = c(M, J, T, n.keep))
u.samples <- array(NA, dim = c(J, M, T, n.keep))
beta.samples <- matrix(NA, nrow = n.keep, ncol = M*P)
eta.samples <- array(NA, dim = c(M, Q, J, n.keep))
Sigma_epsilon.samples <- array(NA, dim = c(M, M, R, n.keep))
Sigma_u.samples <- array(NA, dim = c(M, M, R, T, n.keep))
Sigma_delta.samples <- array(NA, dim = c(M, M, J, n.keep))
tau_sq_eta.samples <- array(NA, dim = c(M, Q, n.keep))
rho_v.samples <- matrix(NA, nrow = n.keep, ncol = M)
rho_eta.samples <- array(NA, dim = c(M, Q, n.keep))

##########################
## Sampler.
##########################

s.keep <- 1

batch.iter <- 0
batch.length <- 10
batch.accept.v <- rep(0, times = M)
batch.accept.eta <- matrix(0, nrow = M, ncol = Q)

time.start <- Sys.time()

for (s in 1:n.iter) {
  # update missing elements of y
  y.s <- y
  for (t in 1:T) {
    for (j in 1:J) {
      if (n_jt[j,paste0(years[t])] > 0) {
        for (i in 1:n_jt[j,paste0(years[t])]){
          missing_y <- which(is.na(y[[t]][[j]][,i]))
          if (length(missing_y) > 0) {
            obs_y <- which(!is.na(y[[t]][[j]][,i]))
            mean_missing_y <- mu.s[missing_y, j, t] + Sigma_delta.s[missing_y, obs_y, j] %*% inv(Sigma_delta.s[obs_y, obs_y, j]) %*% (y[[t]][[j]][obs_y, i] - mu.s[obs_y, j, t])
            Sigma_missing_y <- Sigma_delta.s[missing_y, missing_y, j] - Sigma_delta.s[missing_y, obs_y, j] %*% inv(Sigma_delta.s[obs_y, obs_y, j]) %*% Sigma_delta.s[obs_y, missing_y, j]
            y.s[[t]][[j]][missing_y, i] <- rmvn(n = 1, mu = mean_missing_y, V = Sigma_missing_y)
          }
        }
      }
    }
  }
  
  # update mu
  
  for (t in 1:T) {
    for (j in 1:J) {
      mu.V <- inv(n_jt[j, paste0(years[t])] * inv(Sigma_delta.s[,,j]) + inv(Sigma_epsilon.s[,,r_j[j]]))
      mu.v <- inv(Sigma_epsilon.s[,,r_j[j]]) %*% (u.s[j,,t] + X[,,j,t] %*% beta.s + X_tilde[,,j,t] %*% c(t(eta.s[,,j])))
      if (n_jt[j, paste0(years[t])] > 0) {
        for (i in 1:n_jt[j, paste0(years[t])]) {
          mu.v <- mu.v + inv(Sigma_delta.s[,,j]) %*% y.s[[t]][[j]][,i]
        }
      }
      mu.s[,j,t] <- rmvn(n = 1, mu = mu.V %*% mu.v, V = mu.V)
    }
  }
  # mu.s <- dat$mu_true
  
  # update beta
  
  beta.V <- inv(Sigma_beta)
  beta.v <- inv(Sigma_beta) %*% mu_beta
  for (t in 1:T) {
    for (j in 1:J) {
      beta.V <- beta.V + t(X[,,j,t]) %*% inv(Sigma_epsilon.s[,,r_j[j]]) %*% X[,,j,t]
      beta.v <- beta.v + t(X[,,j,t]) %*% inv(Sigma_epsilon.s[,,r_j[j]]) %*% mu.s[,j,t] -
        t(X[,,j,t]) %*% inv(Sigma_epsilon.s[,,r_j[j]]) %*% u.s[j,,t] -
        t(X[,,j,t]) %*% inv(Sigma_epsilon.s[,,r_j[j]]) %*% X_tilde[,,j,t] %*% c(t(eta.s[,,j]))
    }
  }
  beta.V <- inv(beta.V)
  
  beta.s <- rmvn(n = 1, mu = beta.V %*% beta.v, V = beta.V)
  # beta.s <- dat$beta_true
  
  # update u
  t <- 1
  for (j in 1:J) {
    u.V <- inv(inv(Sigma_epsilon.s[,,r_j[j]]) + inv(1/d_j[j] * Sigma_u.s[,,r_j[j], t]) +
                 inv(1/d_j[j] * Sigma_u.s[,,r_j[j], t+1]))
    u.v <- inv(Sigma_epsilon.s[,,r_j[j]]) %*% (mu.s[,j,t] - X[,,j,t] %*% beta.s - X_tilde[,,j,t] %*% c(t(eta.s[,,j])))
    mu.u.t <- rep(0, times = M)
    mu.u.tp1 <- rep(0, times = M)
    for (k in N_j[[j]]) {
      mu.u.t <- mu.u.t + (1/d_j[j]) * diag(rho_v.s) %*% (u.s[k,,t])
      mu.u.tp1 <- mu.u.tp1 + (1/d_j[j]) * diag(rho_v.s) %*% (u.s[k,,t+1] - u.s[k,,t])
    }
    u.v <- u.v + inv(1/d_j[j] * Sigma_u.s[,,r_j[j], t]) %*% mu.u.t +
      inv(1/d_j[j] * Sigma_u.s[,,r_j[j], t + 1]) %*% (u.s[j,,t+1] - mu.u.tp1)
    u.s[j,,t] <- rmvn(n = 1, mu = u.V %*% u.v, V = u.V)
  }
  
  for (t in 2:(T-1)) {
    for (j in 1:J) {
      u.V <- inv(inv(Sigma_epsilon.s[,,r_j[j]]) + inv(1/d_j[j] * Sigma_u.s[,,r_j[j], t]) +
                   inv(1/d_j[j] * Sigma_u.s[,,r_j[j], t+1]))
      u.v <- inv(Sigma_epsilon.s[,,r_j[j]]) %*% (mu.s[,j,t] - X[,,j,t] %*% beta.s - X_tilde[,,j,t] %*% c(t(eta.s[,,j])))
      mu.u.t <- rep(0, times = M)
      mu.u.tp1 <- rep(0, times = M)
      for (k in N_j[[j]]) {
        mu.u.t <- mu.u.t + (1/d_j[j]) * diag(rho_v.s) %*% (u.s[k,,t] - u.s[k,,t-1])
        mu.u.tp1 <- mu.u.tp1 + (1/d_j[j]) * diag(rho_v.s) %*% (u.s[k,,t+1] - u.s[k,,t])
      }
      u.v <- u.v + inv(1/d_j[j] * Sigma_u.s[,,r_j[j], t]) %*% (u.s[j,,t-1] + mu.u.t) + inv(1/d_j[j] * Sigma_u.s[,,r_j[j], t+1]) %*% (u.s[j,,t+1] - mu.u.tp1)
      u.s[j,,t] <- rmvn(n = 1, mu = u.V %*% u.v, V = u.V)
    }
  }
  
  t <- T
  for (j in 1:J) {
    u.V <- inv(inv(Sigma_epsilon.s[,,r_j[j]]) + inv(1/d_j[j] * Sigma_u.s[,,r_j[j], t]))
    u.v <- inv(Sigma_epsilon.s[,,r_j[j]]) %*% (mu.s[,j,t] - X[,,j,t] %*% beta.s - X_tilde[,,j,t] %*% c(t(eta.s[,,j])))
    mu.u.t <- rep(0, times = M)
    for (k in N_j[[j]]) {
      mu.u.t <- mu.u.t + (1/d_j[j]) * diag(rho_v.s) %*% (u.s[k,,t] - u.s[k,,t-1])
    }
    u.v <- u.v + inv(1/d_j[j] * Sigma_u.s[,,r_j[j], t]) %*% (u.s[j,,t-1] + mu.u.t)
    u.s[j,,t] <- rmvn(n = 1, mu = u.V %*% u.v, V = u.V)
  }
  # u.s <- dat$u_true
  
  # update eta
  for (j in 1:J) {
    Sigma_eta_j <- 1/d_j[j] * diag(c(t(tau_sq_eta.s)))
    mu_eta_j <- rep(0, times = M*Q)
    for (k in N_j[[j]]) {
      mu_eta_j <- mu_eta_j + (1/d_j[j]) * diag(c(t(rho_eta.s))) %*% c(t(eta.s[,,k]))
    }
    eta.V <- inv(Sigma_eta_j)
    eta.v <- inv(Sigma_eta_j) %*% mu_eta_j
    for (t in 1:T) {
      eta.V <- eta.V + t(X_tilde[,,j,t]) %*% inv(Sigma_epsilon.s[,,r_j[j]]) %*% X_tilde[,,j,t]
      eta.v <- eta.v + t(X_tilde[,,j,t]) %*% inv(Sigma_epsilon.s[,,r_j[j]]) %*% mu.s[,j,t] -
        t(X_tilde[,,j,t]) %*% inv(Sigma_epsilon.s[,,r_j[j]]) %*% u.s[j,,t] -
        t(X_tilde[,,j,t]) %*% inv(Sigma_epsilon.s[,,r_j[j]]) %*% X[,,j,t] %*% beta.s
    }
    eta.V <- inv(eta.V)
    
    eta.s[,,j] <- matrix(rmvn(n = 1, mu = eta.V %*% eta.v, V = eta.V),
                         nrow = M, ncol = Q, byrow = TRUE)
  }
  # eta.s <- dat$eta_true
  
  # update Sigma_u
  t <- 1
  for (r in 1:R) {
    nu_tilde <- nu_u + (length(J_r[[r]]))
    H_tilde <- H_u
    for (j in J_r[[r]]) {
      mu.u.t <- rep(0, times = M)
      for (k in N_j[[j]]) {
        mu.u.t <- mu.u.t + (1/d_j[j]) * diag(rho_v.s) %*% u.s[k,,t]
      }
      H_tilde <- H_tilde + d_j[j] * (u.s[j,,t] - mu.u.t) %*% t(u.s[j,,t] - mu.u.t)
    }
    Sigma_u.s[,,r, t] <- riwish(nu_tilde, H_tilde)
  } 
  
  for (t in 2:T) {
    for (r in 1:R) {
      nu_tilde <- nu_u + (length(J_r[[r]]))
      H_tilde <- H_u
      for (j in J_r[[r]]) {
        mu.u.t <- rep(0, times = M)
        for (k in N_j[[j]]) {
          mu.u.t <- mu.u.t + (1/d_j[j]) * diag(rho_v.s) %*% (u.s[k,,t] - u.s[k,,t-1])
        }
        H_tilde <- H_tilde + d_j[j] * (u.s[j,,t] - u.s[j,,t-1] - mu.u.t) %*% t(u.s[j,,t] - u.s[j,,t-1] - mu.u.t)
      }
      Sigma_u.s[,,r,t] <- riwish(nu_tilde, H_tilde)
    } 
  }
  
  # Sigma_u.s <- array(apply(dat$A_true, 3, function(A) A %*% t(A)), dim = c(M, M, R))
  
  # update Sigma_delta
  
  for (j in 1:J) {
    nu_tilde <- nu_delta + sum(n_jt[j, paste0(years)])
    H_tilde <- H_delta
    for (t in 1:T) {
      if (n_jt[j, paste0(years[t])] > 0) {
        for (i in 1:n_jt[j, paste0(years[t])]) {
          H_tilde <- H_tilde + (y.s[[t]][[j]][,i] - mu.s[,j,t]) %*% t(y.s[[t]][[j]][,i] - mu.s[,j,t])
        }
      }
    }
    Sigma_delta.s[,,j] <- riwish(nu_tilde, H_tilde)
  }
  # Sigma_delta.s <- dat$Sigma_delta_true
  
  # update Sigma_epsilon
  
  for (r in 1:R) {
    nu_tilde <- nu_epsilon + (T * length(J_r[[r]]))
    H_tilde <- H_epsilon
    for (t in 1:T) {
      for (j in J_r[[r]]) {
        H_tilde <- H_tilde + (mu.s[,j,t] - u.s[j,,t] - X[,,j,t] %*% beta.s - X_tilde[,,j,t] %*% c(t(eta.s[,,j]))) %*%
          t(mu.s[,j,t] - u.s[j,,t] - X[,,j,t] %*% beta.s - X_tilde[,,j,t] %*% c(t(eta.s[,,j])))
      }
    }
    Sigma_epsilon.s[,,r] <- riwish(nu_tilde, H_tilde)
  }
  # Sigma_epsilon.s <- dat$Sigma_epsilon_true
  
  # update tau_sq_eta
  
  for (m in 1:M) {
    for (q in 1:Q) {
      a_tilde <- a_eta + ((J - 1)/2) # change to J instead of J - 1 when rho neq 1
      Q_rho_m_q <- S.1 - rho_eta.s[m,q] * S.2
      b_tilde <- b_eta + 0.5 * t(eta.s[m, q,]) %*% Q_rho_m_q %*% eta.s[m, q,]
      tau_sq_eta.s[m, q] <- rinvgamma(n = 1, shape = a_tilde, scale = b_tilde)
    }
  }
  
  # tau_sq_eta.s <- dat$tau_sq_true
  
  # update rho_v
  
  # for (m in 1:M) {
  #   rho_v.cand <- rho_v.s
  #   rho_v.cand[m] <- logit.inv(rnorm(1, mean = logit(rho_v.s[m], a_v, b_v), sd = sqrt(rho_v_tuning)), a_v, b_v)
  #   current.ltd <- log(rho_v.s[m] - a_v) + log(b_v - rho_v.s[m])
  #   cand.ltd <- log(rho_v.cand[m] - a_v) + log(b_v - rho_v.cand[m])
  #   t <- 1
  #   for (j in 1:J) {
  #     mu.u.j.t.current <- rep(0, times = M)
  #     mu.u.j.t.cand <- rep(0, times = M)
  #     for (k in N_j[[j]]) {
  #       mu.u.j.t.current <- mu.u.j.t.current + (1/d_j[j]) * diag(c(rho_v.s)) %*% (u.s[k,,t])
  #       mu.u.j.t.cand <- mu.u.j.t.cand + (1/d_j[j]) * diag(c(rho_v.cand)) %*% (u.s[k,,t])
  #     }
  #     current.ltd <- current.ltd - 0.5 * d_j[j] * t(u.s[j,,t] - mu.u.j.t.current) %*%
  #       inv(Sigma_u.s[,,r_j[j]]) %*% (u.s[j,,t] - mu.u.j.t.current)
  #     cand.ltd <- cand.ltd - 0.5 * d_j[j] * t(u.s[j,,t] - mu.u.j.t.cand) %*%
  #       inv(Sigma_u.s[,,r_j[j]]) %*% (u.s[j,,t] - mu.u.j.t.cand)
  #   }
  # 
  #   for (t in 2:T) {
  #     for (j in 1:J) {
  #       mu.u.j.t.current <- rep(0, times = M)
  #       mu.u.j.t.cand <- rep(0, times = M)
  #       for (k in N_j[[j]]) {
  #         mu.u.j.t.current <- mu.u.j.t.current + (1/d_j[j]) * diag(c(rho_v.s)) %*% (u.s[k,,t] - u.s[k,,t-1])
  #         mu.u.j.t.cand <- mu.u.j.t.cand + (1/d_j[j]) * diag(c(rho_v.cand)) %*% (u.s[k,,t] - u.s[k,,t-1])
  #       }
  #       current.ltd <- current.ltd - 0.5 * d_j[j] * t(u.s[j,,t] - u.s[j,,t-1] - mu.u.j.t.current) %*%
  #         inv(Sigma_u.s[,,r_j[j]]) %*% (u.s[j,,t] - u.s[j,,t-1] - mu.u.j.t.current)
  #       cand.ltd <- cand.ltd - 0.5 * d_j[j] * t(u.s[j,,t] - u.s[j,,t-1] - mu.u.j.t.cand) %*%
  #         inv(Sigma_u.s[,,r_j[j]]) %*% (u.s[j,,t] - u.s[j,,t-1] - mu.u.j.t.cand)
  #     }
  #   }
  # 
  #   if(runif(1,0,1) < exp(cand.ltd-current.ltd)){
  #     rho_v.s[m] <- rho_v.cand[m]
  #     batch.accept.v[m] <- batch.accept.v[m]+1
  #   }
  # }
  
  rho_v.s <- rep(1, times = M)
  
  # update rho_eta
  
  # for (m in 1:M) {
  #   for (q in 1:Q) {
  #     
  #   }
  # }
  
  rho_eta.s <- matrix(1, nrow = M, ncol = Q)
  
  # store samples
  
  if(s %in% keep){
    mu.samples[,,,s.keep] <- mu.s
    u.samples[,,,s.keep] <- u.s
    beta.samples[s.keep,] <- beta.s
    eta.samples[,,,s.keep] <- eta.s
    Sigma_epsilon.samples[,,,s.keep] <- Sigma_epsilon.s
    Sigma_u.samples[,,,,s.keep] <- Sigma_u.s
    Sigma_delta.samples[,,,s.keep] <- Sigma_delta.s
    tau_sq_eta.samples[,,s.keep] <- tau_sq_eta.s
    rho_v.samples[s.keep,] <- rho_v.s
    rho_eta.samples[,,s.keep] <- rho_eta.s
    
    s.keep <- s.keep+1
  }
  
  ## Progress and reporting.
  batch.iter <- batch.iter + 1
  
  if(batch.iter == batch.length){
    print(paste("Complete:",round(100*s/n.iter), "%"))
    print(paste("Metrop acceptance v:", 100*batch.accept.v/batch.length))
    print(paste("Metrop acceptance eta:", 100*batch.accept.eta/batch.length))
    print(paste("time:", round(as.numeric(difftime(Sys.time(), time.start, units = "hours")), digits = 3), "hours"))
    print("--------------------------")
    batch.iter <- 0
    batch.accept.v <- rep(0, times = M)
    batch.accept.eta <- matrix(0, nrow = M, ncol = Q)
    # time.start <- Sys.time()
  }
}

apply(tau_sq_eta.samples, c(1,2), FUN = mean)

colMeans(rho_v.samples)
traceplot(as.mcmc(rho_v.samples[,1]), ylim = c(0, 1))

plot(200:250, rho_v.samples[1,5,200:250], col = "red", type = "l", ylim = c(0,1))
points(200:250, rho_v.samples[2,5,200:250], col = "blue", type = "l", ylim = c(0, 1))
points(200:250, rho_v.samples[3,5,200:250], col = "forestgreen", type = "l", ylim = c(0, 1))


r <- 3
apply(Sigma_u.samples, c(1,2,3), FUN = mean)[,,r]

j <- 51
apply(Sigma_delta.samples, c(1,2,3), FUN = mean)[,,j]

r <- 3
apply(Sigma_epsilon.samples, c(1,2,3), FUN = mean)[,,r]

m <- 3
t <- 1
ggplot() + geom_sf(data = shp, aes(fill = apply(mu.samples, c(1,2,3), FUN = mean)[m,,t])) + 
  labs(fill = "fit") +
  scale_fill_viridis_c(limits = range(dat$mu_true[m,,t], apply(mu.samples, c(1,2,3), FUN = mean)[m,,t]))

traceplot(as.mcmc(Sigma_u.samples[2,2,3,]))

m <- 1
t <- 1
ggplot() + geom_sf(data = shp, aes(fill = apply(u.samples, c(1,2,3), FUN = mean)[,m,t])) + 
  labs(fill = "fit") +
  scale_fill_viridis_c(limits = range(dat$u_true[,m,t], apply(u.samples, c(1,2,3), FUN = mean)[,m,t]))

u_quants <- apply(u.samples, c(1,2,3), FUN = quantile, probs = c(0.025, 0.5, 0.975))

for (t in 1:T) {
  p.u.t <- ggplot() + geom_sf(data = shp, aes(fill = u_quants[2,,1,t]), lwd = 0) +
    scale_fill_gradientn(colors = brewer.spectral(n = 100),
                         limits = c(-max(abs(u_quants[2,,1,t])), max(abs(u_quants[2,,1,t])))) +
    labs(fill = paste0("Live u")) +
    theme_void() +
    theme(text = element_text(family = "serif")) +
    ggplot() + geom_sf(data = shp, aes(fill = u_quants[2,,2,t]), lwd = 0) +
    scale_fill_gradientn(colors = brewer.spectral(n = 100),
                         limits = c(-max(abs(u_quants[2,,2,t])), max(abs(u_quants[2,,2,t])))) +
    labs(fill = paste0("Dead u")) +
    theme_void() +
    theme(text = element_text(family = "serif")) +
    ggplot() + geom_sf(data = shp, aes(fill = u_quants[2,,3,t]), lwd = 0) +
    scale_fill_gradientn(colors = brewer.spectral(n = 100),
                         limits = c(-max(abs(u_quants[2,,3,t])), max(abs(u_quants[2,,3,t])))) +
    labs(fill = paste0("CWD u")) +
    theme_void() +
    theme(text = element_text(family = "serif")) +
    plot_annotation(
      title = paste(years[t]),
      theme = theme(
        plot.title = element_text(size = 20, family = "serif", hjust = 0.5)))
  
  print(p.u.t)
  readline(prompt = "Pause. Press <Enter> to continue...")
}

J.colors <- jet(n = J)
m <- 3
plot(years, u_quants[2,1,m,], type = "l", col = J.colors[1], ylim = range(u_quants[2,,m,]))

for (j in 2:J) {
  points(years, u_quants[2,j,m,], type = "l", col = J.colors[j], ylim = range(u_quants[2,,m,]))
}

eta_quants <- apply(eta.samples, c(1,2,3), FUN = quantile, probs = c(0.025, 0.5, 0.975))

p.eta <- ggplot() + geom_sf(data = shp, aes(fill = eta_quants[2,1,1,]), lwd = 0) +
  scale_fill_gradientn(colors = brewer.spectral(n = 100),
                       limits = c(-max(abs(eta_quants[2,1,1,])), max(abs(eta_quants[2,1,1,])))) +
  labs(fill = paste0("Live eta")) +
  theme_void() +
  theme(text = element_text(family = "serif")) +
  ggplot() + geom_sf(data = shp, aes(fill = eta_quants[2,2,1,]), lwd = 0) +
  scale_fill_gradientn(colors = brewer.spectral(n = 100),
                       limits = c(-max(abs(eta_quants[2,2,1,])), max(abs(eta_quants[2,2,1,])))) +
  labs(fill = paste0("Dead eta")) +
  theme_void() +
  theme(text = element_text(family = "serif")) +
  ggplot() + geom_sf(data = shp, aes(fill = eta_quants[2,3,1,]), lwd = 0) +
  scale_fill_gradientn(colors = brewer.spectral(n = 100),
                       limits = c(-max(abs(eta_quants[2,3,1,])), max(abs(eta_quants[2,3,1,])))) +
  labs(fill = paste0("CWD eta")) +
  theme_void() +
  theme(text = element_text(family = "serif")) +
  plot_annotation(
    theme = theme(
      plot.title = element_text(size = 20, family = "serif", hjust = 0.5)))

print(p.eta)

m <- 3
q <- 1
ggplot() + geom_sf(data = shp, aes(fill = apply(eta.samples, c(1,2,3), FUN = mean)[m,q,])) + 
  labs(fill = "fit") +
  scale_fill_viridis_c(limits = range(dat$eta_true[m,q,], apply(eta.samples, c(1,2,3), FUN = mean)[m,q,]))

direct_mu <- array(NA, dim = c(J, M, T))
direct_se <- array(NA, dim = c(J, M, T))

for (t in 1:T) {
  for (j in 1:J) {
    for (m in 1:M) {
      direct_mu[j, m, t] <- mean(y[[t]][[j]][m,], na.rm = TRUE)
      direct_se[j, m, t] <- sd(y[[t]][[j]][m,], na.rm = TRUE) / sum(!is.na(y[[t]][[j]][m,]))
    }
  }
}

mu_quants <- apply(mu.samples, c(1,2,3), FUN = quantile, probs = c(0.025, 0.5, 0.975))
my.pal <- wes_palette("AsteroidCity1", n = 3, type = "discrete")
my.pal.2 <- brewer.paired(n = 8)
mu_sds <- apply(mu.samples, c(1,2,3), FUN = sd)
for (j in 1:J) {
  fip <- st_drop_geometry(shp$fips)[j]
  
  p.separate <- ggplot() + 
    geom_point(aes(x = years - 0.1, y = direct_mu[j,1,], col = "Live Trees Direct"))  +
    geom_linerange(aes(x = years - 0.1, 
                       ymax = direct_mu[j,1,] + 1.96*direct_se[j,1,],
                       ymin = direct_mu[j,1,] - 1.96*direct_se[j,1,],
                       col = "Live Trees Direct")) +
    geom_point(aes(x = years + 0.1, y = mu_quants[2,1,j,], col = "Live Trees Fit"))  +
    geom_linerange(aes(x = years + 0.1, 
                       ymax = mu_quants[3,1,j,],
                       ymin = mu_quants[1,1,j,],
                       col = "Live Trees Fit")) +
    ylim(range(0, direct_mu[j,1,] + 1.96*direct_se[j,1,],
               direct_mu[j,1,] - 1.96*direct_se[j,1,], 
               mu_quants[3,1,j,], mu_quants[1,1,j,], na.rm = TRUE)) +
    theme_minimal() + 
    xlab(element_blank()) +
    ylab("Carbon Mg/ha") +
    theme(plot.title = element_text(hjust = 0.5),
          text = element_text(size = 14, family = "serif"),
          legend.direction = 'horizontal') +
    scale_color_manual(name=element_blank(),
                       breaks=c('Live Trees Direct', 'Live Trees Fit', 'Dead Trees Direct', 'Dead Trees Fit', 'CWD Direct', 'CWD Fit'),
                       values=c('Live Trees Direct'=my.pal.2[3], 'Live Trees Fit' = my.pal.2[4], 'Dead Trees Direct'=my.pal.2[5], 'Dead Trees Fit' = my.pal.2[6], 'CWD Direct'=my.pal.2[7], 'CWD Fit' = my.pal.2[8])) +
    ggplot() + 
    geom_point(aes(x = years - 0.2, y = direct_mu[j,2,], col = "Dead Trees Direct"))  +
    geom_linerange(aes(x = years - 0.2, 
                       ymax = direct_mu[j,2,] + 1.96*direct_se[j,2,],
                       ymin = direct_mu[j,2,] - 1.96*direct_se[j,2,],
                       col = "Dead Trees Direct")) +
    geom_point(aes(x = years - 0.1, y = mu_quants[2,2,j,], col = "Dead Trees Fit"))  +
    geom_linerange(aes(x = years - 0.1, 
                       ymax = mu_quants[3,2,j,],
                       ymin = mu_quants[1,2,j,],
                       col = "Dead Trees Fit")) +
    geom_point(aes(x = years + 0.1, y = direct_mu[j,3,], col = "CWD Direct"))  +
    geom_linerange(aes(x = years + 0.1, 
                       ymax = direct_mu[j,3,] + 1.96*direct_se[j,3,],
                       ymin = direct_mu[j,3,] - 1.96*direct_se[j,3,],
                       col = "CWD Direct")) +
    geom_point(aes(x = years + 0.2, y = mu_quants[2,3,j,], col = "CWD Fit"))  +
    geom_linerange(aes(x = years + 0.2, 
                       ymax = mu_quants[3,3,j,],
                       ymin = mu_quants[1,3,j,],
                       col = "CWD Fit")) +
    theme_minimal() + 
    xlab(element_blank()) +
    ylab(element_blank()) +
    theme(plot.title = element_text(hjust = 0.5),
          text = element_text(size = 14, family = "serif"),
          legend.direction = 'horizontal') +
    scale_color_manual(name=element_blank(),
                       breaks=c('Live Trees Direct', 'Live Trees Fit', 'Dead Trees Direct', 'Dead Trees Fit', 'CWD Direct', 'CWD Fit'),
                       values=c('Live Trees Direct'=my.pal.2[3], 'Live Trees Fit' = my.pal.2[4], 'Dead Trees Direct'=my.pal.2[5], 'Dead Trees Fit' = my.pal.2[6], 'CWD Direct'=my.pal.2[7], 'CWD Fit' = my.pal.2[8])) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = paste0(st_drop_geometry(shp$county)[j], ", ", st_drop_geometry(shp$state)[j]),
      theme = theme(
        legend.position = 'bottom',
        plot.title = element_text(size = 20, family = "serif", hjust = 0.5)))
  
  p.both <- ggplot() + 
    geom_point(aes(x = years - 0.1, y = direct_mu[j,1,], col = "Live Trees Direct"))  +
    geom_linerange(aes(x = years - 0.1, 
                       ymax = direct_mu[j,1,] + 1.96*direct_se[j,1,],
                       ymin = direct_mu[j,1,] - 1.96*direct_se[j,1,],
                       col = "Live Trees Direct")) +
    geom_point(aes(x = years, y = direct_mu[j,2,], col = "Dead Trees Direct"))  +
    geom_linerange(aes(x = years, 
                       ymax = direct_mu[j,2,] + 1.96*direct_se[j,2,],
                       ymin = direct_mu[j,2,] - 1.96*direct_se[j,2,],
                       col = "Dead Trees Direct")) +
    geom_point(aes(x = years + 0.1, y = direct_mu[j,3,], col = "CWD Direct"))  +
    geom_linerange(aes(x = years + 0.1, 
                       ymax = direct_mu[j,3,] + 1.96*direct_se[j,3,],
                       ymin = direct_mu[j,3,] - 1.96*direct_se[j,3,],
                       col = "CWD Direct")) +
    geom_point(aes(x = years - 0.2, y = mu_quants[2,1,j,], col = "Live Trees Fit"))  +
    geom_linerange(aes(x = years - 0.2, 
                       ymax = mu_quants[3,1,j,],
                       ymin = mu_quants[1,1,j,],
                       col = "Live Trees Fit")) +
    geom_point(aes(x = years, y = mu_quants[2,2,j,], col = "Dead Trees Fit"))  +
    geom_linerange(aes(x = years, 
                       ymax = mu_quants[3,2,j,],
                       ymin = mu_quants[1,2,j,],
                       col = "Dead Trees Fit")) +
    geom_point(aes(x = years + 0.2, y = mu_quants[2,3,j,], col = "CWD Fit"))  +
    geom_linerange(aes(x = years + 0.2, 
                       ymax = mu_quants[3,3,j,],
                       ymin = mu_quants[1,3,j,],
                       col = "CWD Fit")) +
    theme_minimal() + 
    ggtitle(paste0(st_drop_geometry(shp$county)[j], ", ", st_drop_geometry(shp$state)[j], " Direct")) +
    xlab(element_blank()) +
    ylab("Carbon Mg/ha") +
    theme(plot.title = element_text(hjust = 0.5),
          text = element_text(size = 14, family = "serif")) +
    scale_color_manual(name=element_blank(),
                       breaks=c('Live Trees Direct', 'Live Trees Fit', 'Dead Trees Direct', 'Dead Trees Fit', 'CWD Direct', 'CWD Fit'),
                       values=c('Live Trees Direct'=my.pal.2[1], 'Live Trees Fit' = my.pal.2[2], 'Dead Trees Direct'=my.pal.2[3], 'Dead Trees Fit' = my.pal.2[4], 'CWD Direct'=my.pal.2[5], 'CWD Fit' = my.pal.2[6]))
  
  print(p.separate)
  readline(prompt = "Pause. Press <Enter> to continue...")
}

m <- 1
plot(years, mu_quants[2,m,1,], type = "l", col = J.colors[1], ylim = range(mu_quants[2,m,,]))

for (j in 2:J) {
  points(years, mu_quants[2,m,j,], type = "l", col = J.colors[j], ylim = range(mu_quants[2,m,,]))
}

for (t in 1:T) {
  direct_t <- as.data.frame(direct_mu[,,t])
  direct_se_t <- as.data.frame(direct_se[,,t])
  colnames(direct_t) <- c("m_live", "m_dead", "m_cwd")
  colnames(direct_se_t) <- c("se_live", "se_dead", "se_cwd")
  means_t <- cbind(shp, direct_t, direct_se_t)
  for (m in 1:M) {
    pools <- c("Live", "Dead", "CWD")
    fill_text <- paste0(pools[m], " C\nMg/ha") 
    maps <- ggplot() + 
      geom_sf(data = means_t, aes(fill = direct_t[,m]), lwd = 0) +
      scale_fill_gradientn(colors = rev(ocean.algae(n = 100)), limits = range(direct_t[,m], mu_quants[2,m,,t], na.rm = TRUE)) +
      theme_void() +
      labs(fill = fill_text) +
      ggtitle("Direct") +
      theme(text = element_text(family = "serif"),
            plot.title = element_text(hjust = 0.5)) +
      ggplot() +
      geom_sf(data = means_t, aes(fill = mu_quants[2,m,,t]), lwd = 0) +
      scale_fill_gradientn(colors = rev(ocean.algae(n = 100)), limits = range(direct_t[,m], mu_quants[2,m,,t], na.rm = TRUE)) +
      theme_void() +
      labs(fill = fill_text) +
      ggtitle("Fit") +
      theme(text = element_text(family = "serif"),
            plot.title = element_text(hjust = 0.5)) +
      ggplot() + 
      geom_sf(data = means_t, aes(fill = direct_se_t[,m]), lwd = 0) +
      scale_fill_gradientn(colors = brewer.ylorrd(n = 100), limits = range(mu_sds[m,,t], mu_quants[3,m,,t] - mu_quants[1,m,,t], na.rm = TRUE)) +
      theme_void() +
      labs(fill = fill_text) +
      ggtitle("Direct SE") +
      theme(text = element_text(family = "serif"),
            plot.title = element_text(hjust = 0.5)) +
      ggplot() +
      geom_sf(data = means_t, aes(fill = mu_sds[m,,t]), lwd = 0) +
      scale_fill_gradientn(colors = brewer.ylorrd(n = 100), limits = range(mu_sds[m,,t], mu_quants[3,m,,t] - mu_quants[1,m,,t], na.rm = TRUE)) +
      theme_void() +
      labs(fill = fill_text) +
      ggtitle("Fit SD") +
      theme(text = element_text(family = "serif"),
            plot.title = element_text(hjust = 0.5)) +
      plot_layout(guides = "collect") +
      plot_annotation(
        title = paste(pools[m], years[t]),
        theme = theme(
          plot.title = element_text(size = 20, family = "serif", hjust = 0.5)))
    
    print(maps)
    readline(prompt = "Pause. Press <Enter> to continue...")
  }
}

# linear trends
trends <- array(NA, dim = c(M, J, n.keep))

for (i in 1:n.keep) {
  for (j in 1:J) {
    for (m in 1:M) {
      tmp <- as.data.frame(mu.samples[m,j,,i])
      colnames(tmp) <- "mu"
      trends[m, j, i] <- unname(coefficients(lm(data = tmp, mu ~ years + 1))[2])
    }
  }
}

trend_quants <-  apply(trends, c(1,2), FUN = quantile, probs = c(0.025, 0.5, 0.975))
which_sig <- matrix(FALSE, nrow = J, ncol = M)
for (j in 1:J) {
  for (m in 1:M) {
    if (trend_quants[1,m,j] * trend_quants[3,m,j] > 0) {
      which_sig[j, m] <- TRUE
    }
  }
}

sig_trend_quants <- matrix(0, nrow = J, ncol = M)

for (m in 1:M) {
  sig_trend_quants[which(which_sig[,m]), m] <- trend_quants[2,m, which(which_sig[,m])]
}

(ggplot() + 
    geom_sf(data = shp, aes(fill = sig_trend_quants[,1]), lwd = 0) +
    scale_fill_gradientn(colors = brewer.spectral(n = 100),
                         limits = c(-max(abs(sig_trend_quants[,1])), max(abs(sig_trend_quants[,1])))) +
    theme_void() +
    labs(fill = "Live C\nMg/ha/year") +
    theme(text = element_text(family = "serif"))) +
  (ggplot() + 
     geom_sf(data = shp, aes(fill = sig_trend_quants[,2]), lwd = 0) +
     scale_fill_gradientn(colors = brewer.spectral(n = 100),
                          limits = c(-max(abs(sig_trend_quants[,2])), max(abs(sig_trend_quants[,2])))) +  theme_void() +
     labs(fill = "Dead C\nMg/ha/year") +
     theme(text = element_text(family = "serif"))) +
  (ggplot() + 
     geom_sf(data = shp, aes(fill = sig_trend_quants[,3]), lwd = 0) +
     scale_fill_gradientn(colors = brewer.spectral(n = 100),
                          limits = c(-max(abs(sig_trend_quants[,3])), max(abs(sig_trend_quants[,3])))) +  theme_void() +
     labs(fill = "CWD C\nMg/ha/year") +
     theme(text = element_text(family = "serif"))) +
  plot_annotation(
    title = "Significant Carbon Trends",
    theme = theme(
      plot.title = element_text(size = 20, family = "serif", hjust = 0.5)))

ggplot() + 
  geom_sf(data = shp, aes(fill = trend_quants[2,1,]), lwd = 0) +
  scale_fill_gradientn(colors = brewer.spectral(n = 100),
                       limits = c(-max(abs(trend_quants[2,1,])), max(abs(trend_quants[2,1,])))) +
  theme_void() +
  labs(fill = "Live C\nMg/ha/year") +
  theme(text = element_text(family = "serif")) +
  ggplot() + 
  geom_sf(data = shp, aes(fill = trend_quants[2,2,]), lwd = 0) +
  scale_fill_gradientn(colors = brewer.spectral(n = 100),
                       limits = c(-max(abs(trend_quants[2,2,])), max(abs(trend_quants[2,2,])))) +  theme_void() +
  labs(fill = "Dead C\nMg/ha/year") +
  theme(text = element_text(family = "serif")) +
  ggplot() + 
  geom_sf(data = shp, aes(fill = trend_quants[2,3,]), lwd = 0) +
  scale_fill_gradientn(colors = brewer.spectral(n = 100),
                       limits = c(-max(abs(trend_quants[2,3,])), max(abs(trend_quants[2,3,])))) +  theme_void() +
  labs(fill = "CWD C\nMg/ha/year") +
  theme(text = element_text(family = "serif")) +
  plot_annotation(
    title = "Carbon Trends",
    theme = theme(
      plot.title = element_text(size = 20, family = "serif", hjust = 0.5)))
