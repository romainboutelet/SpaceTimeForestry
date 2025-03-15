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

tcc_del <- myraster[myraster$delaware_tcc_2 == 255][,1]
mask_del <- c(as.matrix(myraster$delaware_tcc_2 == 255))
coords_del <- crds(myraster)[mask_del,]

mydf_del <- data.frame(x1 = coords_del[,1], x2 = coords_del[,2], z = tcc_del)
head(mydf_del)

# saveRDS(mydf_del, file = "tcc_data/tcc_delaware.rds")
set.seed(1)
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
p_del <- matrix((samples_obs/100)^(log(0.5)/log(0.1)), ncol = 1)
yobs_del <- apply(p_del, 1, 
                  function(x) return(sample(c(0,1), 1, prob = c(1-x, x))))
mydf_del_samples <- data.frame(x1 = samples_coords[,1], x2 = samples_coords[,2],
                               z = samples_tcc, y = yobs_del)

head(mydf_del_samples)

# saveRDS(mydf_del_samples, file = "tcc_data/tcc_del_samples.rds")


### Trying MCMC on TCC data

mydf_del <- readRDS("tcc_data/tcc_delaware.rds")
mydf_del_samples <- readRDS("tcc_data/tcc_del_samples.rds")
X_del <- cbind(mydf_del$x1, mydf_del$x2)
Z_del <- mydf_del$z
xobs_del <- cbind(mydf_del_samples$x1, mydf_del_samples$x2)
yobs_del <- mydf_del_samples$y
zobs_del <- mydf_del_samples$z
nn_del <- nn2(xobs_del, X_del, k = 50)
nn_obs_del <- nn2(xobs_del, xobs_del, k = 50)
dist_obs_del <- unname(as.matrix(dist(xobs_del, diag = T, upper = T)))

MCMC_out_del <- MCMC_all(xobs_del, yobs_del, zobs_del, n_burnin = 5000,
                         n_samples = 10000, thin = 50, batch = 500,
                         beta0_var_prop = 0.1, beta0_start = 0,
                         a_var_prop = 0.1, a_max = 10000,
                         beta1_var_prop = 0.002, 
                         beta1_var_prior = 0.15, beta1_start = 0, dist.used = F)
# saveRDS(MCMC_out_del, file = "spMC_MCMCsamples.rds")
MCMC_out_del <- readRDS("spMC_MCMCsamples.rds")
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



fn <- function(i){
  if (i != 6) idx_seq <- (1:(nrow(X_del) %/% 6)) + (i-1)*nrow(X_del) %/% 6
  else idx_seq <- (1 + (i-1)*nrow(X_del) %/% 6):nrow(X_del)
  return(MCMC_pred(MCMC_out_del$subsample, X_del[idx_seq,,drop = F], xobs_del,
                   yobs_del, z = Z_del[idx_seq],
                   nn = list(nn.idx = nn_del$nn.idx[idx_seq,,drop = F],
                             nn.dists = nn_del$nn.dists[idx_seq,,drop = F]),
                   ProgressFile = paste0("/mcmc_",i)))
}

# p_del_MC <- unlist(mclapply(1:6, fn, mc.cores = 6, mc.preschedule = F))
# saveRDS(p_del_MC, file = "spMC_del_MCMCpredictions.rds")
p_del_MC <- readRDS("spMC_del_MCMCpredictions.rds")

## Now need to run spNNGP  

### Fitting NNGP

sigma.sq <- 5
tau.sq <- 1
phi <- 1/100
##Fit a Response and Latent NNGP model
starting <- list("phi"=phi, "sigma.sq"=5, "tau.sq"=1)
tuning <- list("phi"=1.5)
priors <- list(phi.unif = c(1/1000,1/11), sigma.sq.IG = c(2,5),
               tau.sq.IG = c(2,5))

spNNGP_del <- spNNGP(y ~ z, data = mydf_del_samples,
                   coords = xobs_del, method = "latent",
                   family = "binomial", cov.model = "exponential",
                   priors = priors, starting = starting, tuning = tuning,
                   n.samples = 10000)
# saveRDS(spNNGP_del, file = "spNNGP_del_MCMCsamples.rds")
spNNGP_del <- readRDS("spNNGP_del_MCMCsamples.rds")

plot(spNNGP_del$p.beta.samples)
plot(spNNGP_del$p.theta.samples)


# {
#   n <- 10
#   p_del_NNGP <- rep(0, nrow(X_del))
#   for (i in (1:n)){
#     thin <- 50
#     start <- 5000 + (i-1)*500
#     end <- 5000 + i*500 - 1
#     if (i == n) end <- 5000 + i*500
#     logit_out_del <- predict(spNNGP_del,
#                              matrix(c(rep(1, length(Z_del)),Z_del), ncol = 2),
#                              X_del, sub.sample = list(start = start, end = end,
#                                                     thin = thin), verbose = T,
#                              n.omp.threads = 1, n.report = 5000)
#     pred_del_NNGP <- logit_out_del$p.y.0
#     w_p <- 10
#     if (i == n) w_p <- 11
#     p_del_NNGP <- p_del_NNGP + 
#       apply(1-1/(1+exp(pred_del_NNGP)), 1, mean)/w_p
#     rm(logit_out_del)
#     rm(pred_del_NNGP)
#     gc()
#   }
# }
# saveRDS(p_del_NNGP, file = "spNNGP_del_predictions.rds")
p_del_NNGP <- readRDS("spNNGP_del_predictions.rds")

### Now trying pure GLM
 
model_nospat <- glm(y~z, family = "binomial", data = mydf_del_samples)
p_del_GLM <- as.numeric(1-1/(1+exp(predict(model_nospat, mydf_del))))


### Checking maps

{
  str_name <- "/delaware_tcc.tif"
  myraster <- rast(paste0(str_folder,str_name))
  par(mfrow = c(1,4))
  plot(1*(myraster$delaware_tcc_1>=10))
  myraster_MC <- myraster
  myraster_MC$delaware_tcc_1[myraster_MC$delaware_tcc_2 == 255] <- p_del_MC
  plot(1*(myraster_MC$delaware_tcc_1 >= 0.5))
  myraster_GLM <- myraster
  myraster_GLM$delaware_tcc_1[myraster_GLM$delaware_tcc_2 == 255] <- p_del_GLM
  plot(1*(myraster_GLM$delaware_tcc_1 >= 0.5))
  myraster_NNGP <- myraster
  myraster_NNGP$delaware_tcc_1[myraster_NNGP$delaware_tcc_2 == 255] <- p_del_NNGP
  plot(1*(myraster_NNGP$delaware_tcc_1 >= 0.5))
  rm(myraster,myraster_MC,myraster_GLM,myraster_NNGP)
  gc()
}

{
  str_name <- "/delaware_tcc.tif"
  myraster <- rast(paste0(str_folder,str_name))
  par(mfrow = c(1,3))
  # plot(100*(myraster$delaware_tcc_1/100)^(log(0.5)/log(0.1)))
  myraster_MC <- myraster
  myraster_MC$delaware_tcc_1[myraster_MC$delaware_tcc_2 == 255] <- p_del_MC
  plot(100*(myraster_MC$delaware_tcc_1) + 50 - 
         100*(myraster$delaware_tcc_1/100)^(log(0.5)/log(0.1)))
  myraster_GLM <- myraster
  myraster_GLM$delaware_tcc_1[myraster_GLM$delaware_tcc_2 == 255] <- p_del_GLM
  plot(100*(myraster_GLM$delaware_tcc_1) + 50 -
         100*(myraster$delaware_tcc_1/100)^(log(0.5)/log(0.1)))
  myraster_NNGP <- myraster
  myraster_NNGP$delaware_tcc_1[myraster_NNGP$delaware_tcc_2 == 255] <- p_del_NNGP
  plot(100*(myraster_NNGP$delaware_tcc_1) + 50 -
         100*(myraster$delaware_tcc_1/100)^(log(0.5)/log(0.1)))
  rm(myraster,myraster_MC,myraster_GLM,myraster_NNGP)
  gc()
}



### Diagnostics

mean((mydf_del$z >= 10) != (p_del_GLM >= 0.5))
mean((mydf_del$z >= 10) != (p_del_MC >= 0.5))
mean((mydf_del$z >= 10) != (p_del_NNGP >= 0.5))

##



b0_SMC <- MCMC_out_del$subsample$beta0
b1_SMC <- MCMC_out_del$subsample$beta1
myseq <- matrix(seq(0,100,0.1), ncol = 1)
mycurve_SMC <- apply(myseq,1,function(x) return(mean(1-1/(1+exp(-b0_SMC-b1_SMC*x)))))

idx_NNGP <- seq(5000, 10000, 50)
b0_NNGP <- spNNGP_del$p.beta.samples[idx_NNGP,1]
b1_NNGP <- spNNGP_del$p.beta.samples[idx_NNGP,2]
mycurve_NNGP <- apply(myseq,1,function(x) return(mean(1-1/(1+exp(b0_NNGP+b1_NNGP*x)))))

b0_glm <- model_nospat$coefficients[1]
b1_glm <- model_nospat$coefficients[2]

par(mfrow=c(1,1))
plot(myseq, mycurve_SMC, "l", col = "blue", ylim = c(0,1))
lines(myseq, 1-1/(1+exp(b0_glm+b1_glm*myseq)), col = "green")
lines(myseq, mycurve_NNGP, col = "red")
lines(myseq, (myseq/100)^(log(0.5)/log(0.1)))
lines(myseq, rep(0.5, 1001))
lines(rep(10, 1001), myseq/100)
