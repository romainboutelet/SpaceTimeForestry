library(terra)
library(sf)
library(stars)
library(FIESTA)
library(dplyr)
library(ggplot2)
library(patchwork)
library(exactextractr)
library(spNNGP)
library(parallel)

# Create shape file

str_folder <- '/home/romain/Documents/Research/fia_tcc_CONUS_all'

# Load in CONUS shapefile from FIESTA package
shp <- stunitco[(stunitco$STATENM == "Delaware"), c("STATENM", "COUNTYNM", "COUNTYFIPS")] %>%
  st_union()

plot(shp)

st_write(shp, paste0(str_folder, "/delaware.shp"), driver = "ESRI Shapefile")

## Need to do this in terminal:
#gdalwarp -cutline delaware.shp -crop_to_cutline -dstalpha science_tcc_CONUS_2019_v2021-4/science_tcc_conus_2019_v2021-4.tif delaware_tcc.tif

# Exctracting TCC data

#str_name <- "/science_tcc_CONUS_2019_v2021-4/science_tcc_conus_2019_v2021-4.tif"
str_name <- "/delaware_tcc.tif"

myraster <- rast(paste0(str_folder,str_name))

# plot(myraster)

tcc_del <- myraster[myraster$delaware_tcc_2 == 255][,1]
mask_del <- c(as.matrix(myraster$delaware_tcc_2 == 255))
coords_del <- crds(myraster)[mask_del,]

mydf_del <- data.frame(tcc = tcc_del, x = coords_del[,1],
                       y = coords_del[,2])
head(mydf_del)

# saveRDS(mydf_del, file = "tcc_data/tcc_delaware.rds"))
mydf_del <- readRDS("tcc_data/tcc_delaware.rds")

mygrid <- st_make_grid(st_bbox(shp), cellsize = 5295, square = F)
mygrid_fia <- mygrid[shp]
sample_fia_init <- st_sample(mygrid_fia, size = c(1,1))
gap_fia <- 37
sample_fia_up <- sample_fia_init
for (i in 1:length(sample_fia_init)){
  sample_fia_up[[i]] <- sample_fia_init[[i]] + c(0,gap_fia)
}
sample_fia_left <- sample_fia_init
for (i in 1:length(sample_fia_init)){
  sample_fia_left[[i]] <- sample_fia_init[[i]]+c(-gap_fia*cos(pi/6),
                                                 -gap_fia*sin(pi/6))
}
sample_fia_right <- sample_fia_init
for (i in 1:length(sample_fia_init)){
  sample_fia_right[[i]] <- sample_fia_init[[i]]+c(gap_fia*cos(pi/6),
                                                  -gap_fia*sin(pi/6))
}
sample_fia <- c(sample_fia_init, sample_fia_up, sample_fia_left,
                sample_fia_right)[shp]

par(mfrow = c(1,1))
plot(mygrid)
plot(st_geometry(shp), add = T)
plot(mygrid_fia, col = "#ff000088", add = T)
plot(sample_fia, cex = 0.05, add = T)

samples_coords <- unname(st_coordinates(sample_fia))
samples_obs <- exact_extract(myraster, st_buffer(sample_fia, 30),
                             "mean")$mean.delaware_tcc_1
samples_tcc <- as.numeric(extract(myraster, vect(samples_coords))$delaware_tcc_1)

mydf_del_samples <- data.frame(x = samples_coords[,1], y = samples_coords[,2],
                               obs = samples_obs, z = samples_tcc)

head(mydf_del_samples)

# saveRDS(mydf_del_samples, file = "tcc_data/tcc_del_samples.rds")
mydf_del_samples <- readRDS("tcc_data/tcc_del_samples.rds")


### Trying MCMC on TCC data

X_del <- cbind(mydf_del$x, mydf_del$y)
Z_del <- mydf_del$tcc
xobs_del <- cbind(mydf_del_samples$x, mydf_del_samples$y)
p_del <- matrix((mydf_del_samples$obs/100)^(log(0.5)/log(0.1)), ncol = 1)

yobs_del <- apply(p_del, 1, 
                  function(x) return(sample(c(0,1), 1, prob = c(1-x, x))))
mydf_del_samples$y <- yobs_del
zobs_del <- mydf_del_samples$z
nn_del <- nn2(xobs_del, X_del, k = 20)
nn_obs_del <- nn2(xobs_del, xobs_del, k = 20)

MCMC_out_del <- MCMC_all(xobs_del, yobs_del, zobs_del, nn_obs_del,
                         beta0_var_prop = 0.15, beta1_var_prop = 0.005,
                         a_var_prop = 0.15, beta1_var_prior = 0.15)
# saveRDS(MCMC_out_del, file = "spMC_MCMCsamples.rds")
# MCMC_out_del <- readRDS("spMC_MCMCsamples.rds")
ggplot(MCMC_out_del$diagnostics, aes(x = x, y = chain)) +
  geom_line() + 
  facet_wrap(~factor(par_name), scales = "free")
p1 <- ggplot(MCMC_out_del$subsample, aes(beta0)) +
  geom_density() 
p2 <- ggplot(MCMC_out_del$subsample, aes(beta1)) +
  geom_density() 
p1 + p2
p3 <- ggplot(MCMC_out_del$subsample, aes(a1)) +
  geom_density() 
p4 <- ggplot(MCMC_out_del$subsample, aes(a2)) +
  geom_density() 
p5 <- ggplot(MCMC_out_del$subsample, aes(a3)) +
  geom_density() 
p3 + p4 + p5
b0_SMC <- MCMC_out_del$subsample$beta0
b1_SMC <- MCMC_out_del$subsample$beta1
myseq <- matrix(seq(0,100,0.1), ncol = 1)
mycurve_SMC <- apply(myseq,1,function(x) return(mean(1-1/(1+exp(-b0_SMC-b1_SMC*x)))))
plot(myseq, mycurve_SMC, "l")
lines(myseq, rep(0.5, 1001), col = "red")
lines(rep(10, 1001), myseq/100, col = "red")

fn <- function(i){
  if (i != 8) idx_seq <- (1:(nrow(X_del) %/% 8)) + (i-1)*nrow(X_del) %/% 8
  else idx_seq <- (1 + (i-1)*nrow(X_del) %/% 8):nrow(X_del)
  return(MCMC_pred(MCMC_out_del$subsample, X_del[idx_seq,,drop = F], xobs_del,
                   yobs_del, z = Z_del[idx_seq],
                   nn = list(nn.idx = nn_del$nn.idx[idx_seq,,drop = F],
                             nn.dists = nn_del$nn.dists[idx_seq,,drop = F]),
                   ProgressFile = paste0("/mcmc_",i)))
}

p_del_MCMC <- unlist(mclapply(1:8, fn, mc.cores = 8, mc.preschedule = F))
# saveRDS(p_del_MCMC, file = "spMC_MCMCpredictions.rds")



par(mfrow = c(1,2))
plot(1*(myraster$delaware_tcc_1>=10))
myraster_MCMC <- myraster
# myraster_MCMC$delaware_tcc_1[myraster_MCMC$delaware_tcc_2 == 255] <-
#   160*p_del_MCMC^2 - 60*p_del_MCMC
myraster_MCMC$delaware_tcc_1[myraster_MCMC$delaware_tcc_2 == 255] <- p_del_MCMC
plot(1*(myraster_MCMC$delaware_tcc_1 >= 0.5))

## Now need to run spNNGP  

### Fitting NNGP

sigma.sq <- 5
tau.sq <- 1
phi <- 3/100
##Fit a Response and Latent NNGP model
starting <- list("phi"=phi, "sigma.sq"=5, "tau.sq"=1)
tuning <- list("phi"=1.5)
priors <- list(phi.unif = c(3/5000,3/37), sigma.sq.IG = c(2,5),
               tau.sq.IG = c(2,5))

my_logit_del <- spNNGP(y ~ z, data = mydf_del_samples,
                   coords = xobs_del, method = "latent",
                   family = "binomial", cov.model = "exponential", 
                   priors = priors, starting = starting, tuning = tuning,
                   n.samples = 10000)

plot(my_logit_del$p.beta.samples)
plot((my_logit_del$p.theta.samples))

idx_NNGP <- seq(5000, 10000, 50)
b0_NNGP <- my_logit_del$p.beta.samples[idx_NNGP,1]
b1_NNGP <- my_logit_del$p.beta.samples[idx_NNGP,2]
myseq <- matrix(seq(0,100,0.1), ncol = 1)
mycurve_NNGP <- apply(myseq,1,function(x) return(mean(1-1/(1+exp(b0_NNGP+b1_NNGP*x)))))

plot(myseq, mycurve_NNGP, "l")
lines(myseq, rep(0.5, 1001), col = "red")
lines(rep(10, 1001), myseq/100, col = "red")

fn_spNNGP <- function(i){
  if (i != 8) idx_seq <- (1:(nrow(X_del) %/% 8)) + (i-1)*nrow(X_del) %/% 8
  else idx_seq <- (1 + (i-1)*nrow(X_del) %/% 8):nrow(X_del)
  myX <- X_del[idx_seq,,drop = F]
  myZ <- Z_del[idx_seq]
  return(predict(my_logit_del,
                  matrix(c(rep(1, nrow(myX)),myZ), ncol = 2),
                  myX, sub.sample = list(start = 5000, end = 10000,
                                           thin = 50), verbose = T))
}

# logit_out_del <- unlist(mclapply(1:8, fn_spNNGP, mc.cores = 8,
#                                  mc.preschedule = F))
logit_out_del <- predict(my_logit_del,
                         matrix(c(rep(1, nrow(X_del)),Z_del), ncol = 2),
                         X_del, sub.sample = list(start = 5000, end = 10000,
                                                thin = 50), verbose = T,
                         n.omp.threads = 8)
 
