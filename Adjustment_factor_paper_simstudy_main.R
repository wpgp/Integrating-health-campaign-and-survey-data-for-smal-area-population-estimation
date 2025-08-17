#-----------------------------------------------------------------------------
# Adjustment factor-based small area population modelling - SIMULATION STUDY
# Developer: Dr Chris Nnanatu, WorldPop, University of Southampton
# February, 2025
#----------------------------------------------------------------------------
rm(list=ls()) #----Clear the workspace
packages <- c("raster", "haven", "sf","sp", "tmap","tmaptools","tidyverse","terra",
              "lattice", "gridExtra", "devtools", "rlang", "viridis", "spdep",
              "car", "MASS", "maps", "spData", "spDataLarge", "ceramic", "basemaps", "ggmap")
if(length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages()))) }


#Installing INLA!!
if(length(setdiff("INLA", rownames(installed.packages()))) > 0){
  install.packages("INLA", type="binary", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
}
library(INLA)
lapply(packages, library, character.only = TRUE) ##--access the libraries


##----------------------------------------------------------------------------------------------------------------------------------
# Here, we will fit the scale factor model first and then use it to predict population count
# In provinces where there are no Mcensus vs PNLP overlap, scale factors are first predicted for the province
# And then used to predict the gridded population estimates for that province
#-----------------------------------------------------------------------------------------------------------------------------------
# Set working directory for input and output files
path <- "//worldpop.files.soton.ac.uk/Worldpop/Projects/WP517763_GRID3_Phase2/Working/DRC/CHRIS_N/province"
data_path <- paste0(path, "/data/")#---paths for survey data:
out_path <- paste0(path , "/Haut-Lomami/paper/simstudy/output_main")



#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##-------------------STEPS:
#--1) Simulate Grid cells pop counts with a set of specified parameters and covariates
#--2) Aggregate to clusters within provinces
#--3) Aggregate to provinces
#--4) Calculate scale factor where there are Microcensus overlap with the pop counts
#--5) Fit scale factor model at cluster level and predict across the entire grid cells where no overlap exists
#--6) Fit pop model at cluster level and predict across the entire grid cells for both scaled and unscaled
#--7) Carry out model checks and cross validations, aggregate and compare
#--8) Withhold data for a couple of provinces and refit
#--9) Predict pop count for all grid cells including the unsampled provinces for bath scaled and unscaled
#--10) Repeat 1 - 10 for a p% biased population with each cluster biased with b% for p and b in 0.1 to 1.



#---Select Admin unit Size
##--Province
pr=5 #- 25 provinces
prov <- as(raster(nrow=pr, ncol=pr, xmn=0, xmx=1, ymn=0, ymx=1), "SpatialPolygons")
proj4string(prov) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
plot(prov, lwd=3)


##--Cluster
cl = 30#--900 clusters with 36 clusters per province
clust <- as(raster(nrow=cl, ncol=cl, xmn=0, xmx=1, ymn=0, ymx=1), "SpatialPolygons")
proj4string(clust) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
plot(clust, lwd=1, col="green")


#-- grid cells
pg = 120 #--14400 prediction grid cells with 16 grid cells per cluster
#-- 16 times 36 (576) grid cells in each of the 25 provinces
grid = raster(nrow=pg, ncol=pg, xmn=0, xmx=1, ymn=0, ymx=1) #--prediction grid
res(grid) = c(0.0083, 0.0083) #--approx 1km by 1km resolution
ncell(grid)

##--Vectorize the prediction raster
coords.grid = xyFromCell(grid, 1:ncell(grid), spatial=FALSE)


##---Visualise
plot(coords.grid, col="green",
     pch=16, cex.axis=1.5)
#plot(clust, lwd=1.5)
plot(clust, add=T, lwd=1.5)
#plot(prov, lwd=3)
plot(prov, lwd=3, add=T)


##---Convert all polygons and grid cells to points
coords.prov <- coordinates(prov)
coords.clust <- coordinates(clust)
coords.grid <- coordinates(grid)


###--Various combined plots
plot(prov, lwd=2)
points(coords.clust, col="green", pch=15)
points(coords.grid, col="cyan")


####----Generate id of grid cells within the provinces
c.g=SpatialPoints(coords.grid)
proj4string(c.g) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
gp <- rep(NA, nrow(coords.grid))
for(i in 1:length(prov)){
  gp[as.vector(which(!is.na(over(c.g, prov[i]))))] <- i
}
grid4prov <- gp
table(grid4prov)
unique(grid4prov)


####----Generate id of grid cells within the Clusters
c.g=SpatialPoints(coords.grid)
proj4string(c.g) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
gc <- rep(NA, nrow(coords.grid))
for(i in 1:length(clust)){
  gc[as.vector(which(!is.na(over(c.g, clust[i]))))] <- i
}
grid4clust <- gc
table(grid4clust)
unique(grid4clust)


####----Generate id of cluster within provinces
c.c=SpatialPoints(coords.clust)
proj4string(c.c) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
cp <- rep(NA, nrow(coords.clust))
for(i in 1:length(prov)){
  cp[as.vector(which(!is.na(over(c.c, prov[i]))))] <- i
}
clust4prov <- cp
table(clust4prov)
unique(clust4prov)

##--rename
grid_id <- 1:ncell(grid)
prov_id <- grid4prov
clust_id <- grid4clust


##----Check data for consistency
df1 <- data.frame(grd_id = grid_id,
                  clst_id = clust_id,
                  prv_id = prov_id,
                  lon = coords.grid[,1],
                  lat = coords.grid[,2])
df1[order(df1$clst_id),] #--passed
plot(df1$lon, df1$lat)
par(mfrow=c(1,1))

dat_all <- df1
dat.grid <- df1

####------Build the mesh
bnd_grd <- inla.nonconvex.hull(coords.grid, -0.035, -0.05, resolution = c(100, 100))
mesh.grd<- inla.mesh.2d(boundary = bnd_grd,
                        offset=c(0.01, 0.15), max.edge=c(0.01,1),
                        cutoff = 0.0303)
plot(mesh.grd)
plot(clust,add=T, col="green", lwd=2)
points(dat.grid$lon, dat.grid$lat, col="red", pch=15, cex=0.5)
plot(prov, add=T, lwd=3)
plot(mesh.grd, add=T)
mesh.grd$n  ##999 nodes

#plot(dat.grid$lon, dat.grid$lat, col="red", pch=15, cex=0.5)
#plot(prov, lwd=3)
#--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##---Specify SPDE parameters
r0 <- 0.3 #--range
nu <- 1  #--smooth parameter
sigma0 <- 1 #--marginal variance
kappa0 <- sqrt(8*nu)/r0 #--scale parameters
tau0 <- 1/(sqrt(4*pi)*kappa0*sigma0) #--precision

#---SPDE
spde.grd <- inla.spde2.matern(mesh.grd,
                              B.tau = matrix(c(log(tau0), -1, +1), nrow=1, ncol=3),
                              B.kappa = matrix(c(log(kappa0), 0, -1), nrow=1, ncol=3),
                              theta.prior.mean = c(0,0),
                              theta.prior.prec = c(0.1, 0.1))



#--Precision Matrix
Q <- inla.spde2.precision(spde=spde.grd, theta=c(0,0))


#---Simulate the GMRF
sam <- as.vector(inla.qsample(
  n = 1, Q = Q, seed=100))
#length(sam)


###---Build projector matrix A
A <- inla.spde.make.A(mesh=mesh.grd, loc=coords.grid);dim(A)
S.pred <- as.vector(A %*% sam)
#hist(S.pred)


#####---Spatial random effects
library(MASS)
tau2 =1
dim(Q); class(Q)
Q <- as.matrix(Q)
sig <- solve(Q)
mm <- mesh.grd$n
phi.a = mvrnorm(1,rep(0,mm), sig)
phi.a = as.vector(phi.a)  #

#Construct phi for each grid cell
phi.grid <- numeric()
for(i in 1:length(grid4prov)) phi.grid[i] <- phi.a[grid4prov[i]]


##----specify the observation indices for estimation
iset.grd <- inla.spde.make.index(name = "s", spde.grd$n.spde)
#length(iset)


##---Nugget effect/iid term for population count
sigma_e <- 0.05
ee.pred <- rnorm(nrow(dat_all), 0, sigma_e)


#-- generate the Covariates
###-----Simulate covariates ----------------
nn <- nrow(dat_all)
covs = cbind(1,runif(nn, 1, 5),abs(rnorm(nn, 4, 9)),rpois(nn,6),
             abs(rnorm(nn, 7, 9)), runif(nn, 2, 6))#--GEOSPATIAL COVARIATES
dim(zpred <- covs)
head(zpred)
dim(dat_all)

##--Add to dataset
dim(ddat <- cbind(dat_all, covs))
ddat <- data.frame(ddat)
dim(ddat)
head(ddat)


names(ddat)[6:11] <- paste0("x", 1:6)
head(ddat)
#names(ddat)

###############
#Parameters
##----Nugget effect
sigma_e <- 0.344 ##--
eps <- rnorm(nrow(ddat), 0, sigma_e)

###---Simulate Province level effects
n.prov <- length(unique(ddat$prv_id))
prov_sigma <- diag(runif(n.prov, 0.0001, 0.0005))
prov_re <- mvrnorm(1, rep(0, n.prov), prov_sigma)


##--extract province ranef for each grid cell
rand_eff <- function(dat, ranef, d2)
{
  effect <- rep(1, nrow(dat))
  uniq <- unique(d2)

  for(i in  1:nrow(dat))
  {
    for(j in 1:length(uniq))
    {
      if(d2[i]==uniq[j]) effect[i] = ranef[j]
    }
  }
  effect
}
prov_re_all <- rand_eff(ddat, prov_re,
                        ddat$prv_id)

#zpred[,2:6] <- apply(zpred[,2:6], 2, stdize)

###---Simulate building counts
set.seed(505)
sigma_eb <- 0.32
epsb <- rnorm(nrow(ddat), 0, sigma_eb)
sd(epsb)

betaB <- c(4.020, 0.250, 0.0211, 0.012, 0.110, 0.105) #---betas- fixed effects

bld <- lambdaB <- numeric(nn) #

for (i in 1:nn)
{
  lambdaB[i] <- exp(zpred[i,1]*betaB[1] + zpred[i,2]*betaB[2] +
                      zpred[i,3]*betaB[3] + zpred[i,4]*betaB[4] +
                      zpred[i,5]*betaB[5] + zpred[i,6]*betaB[6] + prov_re_all[i]  + epsb[i])
  bld[i] <- rpois(1, lambdaB[i])
}
bld
min(bld)
hist(bld); hist(log(bld))
mean(bld); sd(bld); #--approximately ave of 50 bldgs per cell
length(bld[bld==0])
#bld[bld==0]=1



# Simulate Population Counts
sigma_ep <-0.001
epsp <- rnorm(nrow(ddat), 0, sigma_ep)
sd(epsp)

betaP <- c(3.711, 0.433, 0.104, 0.062, 0.033, 0.082) #--betas - fixed effects coeffs
##-------------
pop <- lambdaP <- numeric(nn)
for (i in 1:nn)
{
  lambdaP[i] <- exp(zpred[i,1]*betaP[1] + zpred[i,2]*betaP[2] +
                      zpred[i,3]*betaP[3] + zpred[i,4]*betaP[4] +
                      zpred[i,5]*betaP[5] + zpred[i,6]*betaP[6] + prov_re_all[i] + S.pred[i] + epsp[i])
  pop[i] <- rpois(1, lambdaP[i])
}
#pop
hist(pop);mean(pop); sd(pop)
sum(pop) 
boxplot(pop)


###--------Add to dataset
min(bld); min(pop); max(bld); max(pop)
ddat$bld <- bld
ddat$pop <- pop
ddat$resp <- pop
ddat$dens <- ddat$pop/ddat$bld #----population density
hist(ddat$dens); hist(log(ddat$dens))


###
class(ddat)
dim(ddat)
#----Scale covariates for stack
stdize2 <- function(x)
{
  #stdz <- (x - mean(x, na.rm=T))/sd(x, na.rm=T)
  stdz <- (x - min(x, na.rm=T))/(max(x, na.rm=T)-min(x, na.rm=T))
  return(stdz)
}


#=-----------------------------------------------------------------------------------
#    Run simulation to address overcounts
#-----------------------------------------------------------------------------------
coverp <- c(0.1, 0.2, 0.3, 0.4, 0.5,0.6, 0.7,0.8,0.9,1) # proportion of biased population
coverb <- c(0.1, 0.2, 0.3, 0.4, 0.5,0.6, 0.7,0.8,0.9, 1) # magnitude of bias

for(i in 1:length(coverp))
{

  # i =1
  result_path1 <- paste0(out_path,"/outputs_for_", coverp[i]*100,"%","_biased_sample")
  if (file.exists(result_path1)){
    setwd(file.path(result_path1))
  } else {
    dir.create(file.path(result_path1))
    setwd(file.path(result_path1))
  }

  pm <- coverp[i]
  print(pm*100)#--

  for(j in 1:length(coverb))
  {

    # j=4
    bm = coverb[j]
    result_path2 <- paste0(result_path1,"/", coverb[j]*100, "%","_bias")
    if (file.exists(result_path2)){
      setwd(file.path(result_path2))
    } else {
      dir.create(file.path(result_path2))
      setwd(file.path(result_path2))
    }
    print(c(pm*100,bm*100))#---

    ##--Aggregate Data to cluster level
    # Calculate/extract the cluster level variables
    dat.grid <- data.frame(ddat)
    dim(dat.clust <- data.frame(dat.grid %>%
                                  group_by(clst_id,prv_id) %>%
                                  dplyr::summarise(bld = sum(bld),
                                            pop = sum(pop),
                                            x1 = mean(x1),
                                            x2 = mean(x2),
                                            x3 = mean(x3),
                                            x4 = mean(x4),
                                            x5 = mean(x5),
                                            x6 = mean(x6))))


    # Apply missing values and inflated values to the cluster data
    dat.clust$popm <- dat.clust$pop
    ind.obs <- sample(nrow(dat.clust), pm*nrow(dat.clust))

    dat.clust$popm[ind.obs] = dat.clust$pop[ind.obs] + round(dat.clust$popm[ind.obs]*bm)

    # Add scale variable to the cluster data
    dat.clust$scale <- dat.clust$pop/dat.clust$popm


    ##---Add lon - lat to cluster data
    dat.clust$lon <- coords.clust[,1]
    dat.clust$lat <- coords.clust[,2]
    with(dat.clust, plot(lon, lat, pch=16, col="red"))
    coords.clust <- cbind(dat.clust$lon, dat.clust$lat)

    ###--Cluster level mesh, spde and projection matrix
    bnd_clust <- inla.nonconvex.hull(coords.clust, -0.035, -0.05, resolution = c(100, 100))
    mesh.clst <- inla.mesh.2d(boundary =bnd_clust, max.edge=c(0.01,1),
                              offset = c(0.01, 0.15),
                              cutoff = 0.0309)

    par(mfrow=c(1,1))
    plot(mesh.clst)
    with(dat.clust, points(lon, lat, pch=16, col="red"))
    plot(prov, lwd=2, add=T)
    mesh.clst$n
    plot(prov, lwd=2)
    ###---Build projector matrix A for clusters
    A.clst<-inla.spde.make.A(mesh=mesh.clst,loc=as.matrix(coords.clust));dim(A.clst)

    ##---Create the SPDE for clusters
    spde.clst <- inla.spde2.matern(mesh.clst, alpha=2)

    ##----specify the observation indices for estimation  for clusters
    iset.clst <- inla.spde.make.index(name = "s", spde.clst$n.spde)

    #@@@@@@@@@@------------------------------------------------------------------------------
    library(car)
    dat.fit.clst <- dat.clust

    ## Covariates scaling and data stacking for both grid and cluster
    dat.fit.clst[, c("x2", "x3", "x4", "x5", "x6")] <- apply(dat.fit.clst[, c("x2",
                                                                              "x3", "x4", "x5", "x6")], 2, stdize2)

    #---
    #names(dat.fit.grd);
    names(dat.fit.clst)
    dat.fit.clst$prv_id2 <- factor(dat.fit.clst$prv_id)
    dat.fit.clst$ID <- 1:nrow(dat.fit.clst)
    covars.fit.clst <- dat.fit.clst[,c("ID","x2", "x3", "x4", "x5", "x6", "bld", "clst_id",
                                       "prv_id", "prv_id2")]; dim(covars.fit.clst) ##---Population density


    #---Build the stacks-----------------------------------
    ##-- for scale  at cluster leve;
    stk.clst.sc <- inla.stack(data=list(y=dat.fit.clst$scale), #the response

                              A=list(A.clst,1),  #the A matrix; the 1 is included to make the list(covariates)

                              effects=list(c(list(Intercept=1), #the Intercept
                                             iset.clst),  #the spatial index
                                           #the covariates
                                           list(covars.fit.clst)
                              ),
                              #this is a quick name so you can call upon easily
                              tag='est')


    ##-- for density at cluster level
    dat.fit.clst$pop2 <- dat.fit.clst$popm
    dat.fit.clst$pop2[dat.fit.clst$pop2==0]= NA
    stk.clst.dens <- inla.stack(data=list(y=dat.fit.clst$pop2/dat.fit.clst$bld), #the response

                                A=list(A.clst,1),  #the A matrix; the 1 is included to make the list(covariates)

                                effects=list(c(list(Intercept=1), #the Intercept
                                               iset.clst),  #the spatial index
                                             #the covariates
                                             list(covars.fit.clst)
                                ),
                                #this is a quick name so you can call upon easily
                                tag='est')

# specify and fit INLA models
    
    # scale factor
    fsc <- y ~ -1 + Intercept + x2 + x3 + x4 + x5 + x6 + f(s, model=spde.clst) +
      f(prv_id, model='iid') +
      f(clst_id, model='iid')
    msc <-inla(fsc, #the formula
               data=inla.stack.data(stk.clst.sc,spde=spde.clst),  #the data stack
               family= 'gamma',   #which family the data comes from
               control.predictor=list(A=inla.stack.A(stk.clst.sc),compute=TRUE),  #compute gives you the marginals of the linear predictor
               control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
               verbose = FALSE)
    summary(msc)

    ind.clst.sc <-inla.stack.index(stk.clst.sc, "est")$data
    fitsc <- exp(msc$summary.linear.predictor[ind.clst.sc,"mean"])


    # population density
    fc3 <- y ~ -1 + Intercept + x2 + x3 + x4 + x5 + x6 + f(s, model=spde.clst) +
      f(prv_id, model='iid') +
      f(clst_id, model='iid')
    mc3 <-inla(fc3, #the formula
               data=inla.stack.data(stk.clst.dens,spde=spde.clst),  #the data stack
               family= 'gamma',   #which family the data comes from
               control.predictor=list(A=inla.stack.A(stk.clst.dens),compute=TRUE),  #compute gives you the marginals of the linear predictor
               control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
               verbose = FALSE)
    summary(mc3)



## Simulation and predictions at grid cells based on the cluster level parameters
    dat_pred <- dat.grid
    dat_pred[, c("x2", "x3", "x4", "x5", "x6")] <- apply(dat_pred[, c("x2",
                                                                      "x3", "x4", "x5", "x6")], 2, stdize2)

    # Grid predictions for scale factor
    A_pred_sc <- inla.spde.make.A(mesh = mesh.clst, loc = cbind(dat_pred$lon, dat_pred$lat))
    dim(A_pred_sc)

    sim_scale <- function(model, dat, Aprediction, run)
    {
      fixedeff  <- scale_hat <- matrix(0, nrow=nrow(dat), ncol = run)
      inla.seed = as.integer(runif(1)*.Machine$integer.max)
      #inla.seed =  1657687559
      set.seed(inla.seed)
      print(inla.seed)
      m1.samp <- inla.posterior.sample(run, model,
                                       seed = inla.seed,
                                       selection=list(x2=1,
                                                      x3=1,
                                                      x4=1,
                                                      x5=1,
                                                      x6=1),num.threads="1:1")

      sfield_nodes_mean <- model$summary.random$s['mean']
      field_mean <- (Aprediction%*% as.data.frame(sfield_nodes_mean)[, 1])
      for(i in 1:run)
      {
        #
        fixedeff[,i] <-
          model$summary.fixed['Intercept', 'mean'] +
          m1.samp[[i]]$latent[1,] * dat[,'x2'] +
          m1.samp[[i]]$latent[2,] * dat[,'x3'] +
          m1.samp[[i]]$latent[3,] * dat[,'x4'] +
          m1.samp[[i]]$latent[4,] * dat[,'x5'] +
          m1.samp[[i]]$latent[5,] * dat[,'x6'] +
          dat$clust_ranef +
          dat$prov_ranef +
          field_mean[,1]

        scale_hat[,i]<- exp(fixedeff[,i])
      }

      dat$mean_scale_hat <- apply(scale_hat, 1, mean, na.rm=T) # mean
      dat$lower_scale_hat <- apply(scale_hat, 1, quantile, probs=c(0.025), na.rm=T) # lower
      dat$upper_scale_hat <- apply(scale_hat, 1, quantile, probs=c(0.975), na.rm=T) # upper
      dat$uncert_scale_hat <- (dat$upper_scale_hat - dat$lower_scale_hat)/dat$mean_scale_hat# uncertainty


      output <- list(scale_hat = scale_hat,
                     est_data = dat)

    }


    # Grid Predictions for cluster populations
    A_pred <- inla.spde.make.A(mesh = mesh.clst, loc = cbind(dat_pred$lon, dat_pred$lat))
    dim(A_pred)

    # Population density
    post_dens <- function(model, dat, Aprediction, scale, run)
    {
      fixedeff <- dens_hat <- pop_hat <- pop_hat_sc1 <- matrix(0, nrow=nrow(dat), ncol = run)
      inla.seed = as.integer(runif(1)*.Machine$integer.max)
      #inla.seed =

      set.seed(inla.seed)
      print(inla.seed)
      m1.samp <- inla.posterior.sample(run, model,
                                       seed = inla.seed,
                                       selection=list(x2=1,
                                                      x3=1,
                                                      x4=1,
                                                      x5=1,
                                                      x6=1),num.threads="1:1")

      sfield_nodes_mean <- model$summary.random$s['mean']
      field_mean <- (Aprediction%*% as.data.frame(sfield_nodes_mean)[, 1])
      for(i in 1:run)
      {
        fixedeff[,i] <-
          model$summary.fixed['Intercept', 'mean'] +
          m1.samp[[i]]$latent[1,] * dat[,'x2'] +
          m1.samp[[i]]$latent[2,] * dat[,'x3'] +
          m1.samp[[i]]$latent[3,] * dat[,'x4'] +
          m1.samp[[i]]$latent[4,] * dat[,'x5'] +
          m1.samp[[i]]$latent[5,] * dat[,'x6'] +
          dat$clust_ranef +
          dat$prov_ranef +
          field_mean[,1]

        ##
        dens_hat[,i] <- exp(fixedeff[,i]) #--density
        pop_hat[,i] <- dens_hat[,i]*dat$bld #--population count
      }

      #--Unscaled
      dat$mean_dens_hat <- apply(dens_hat, 1, mean, na.rm=T) #
      dat$mean_pop_hat <- apply(pop_hat, 1, mean, na.rm=T) #
      dat$lower_pop_hat <- apply(pop_hat, 1, quantile, probs=c(0.025), na.rm=T) #
      dat$upper_pop_hat <- apply(pop_hat, 1, quantile, probs=c(0.975), na.rm=T) #
      dat$uncert_pop_hat <- (dat$upper_pop_hat - dat$lower_pop_hat)/dat$mean_pop_hat#


      #--Scale 1 #--Scaled
      pop_hat_sc1 <- pop_hat*scale # allows for scale-pop uncertainty propagation
      dat$mean_pop_hat_sc1 <- apply(pop_hat_sc1, 1, mean, na.rm=T) #
      dat$lower_pop_hat_sc1 <- apply(pop_hat_sc1, 1, quantile, probs=c(0.025), na.rm=T) #
      dat$upper_pop_hat_sc1 <- apply(pop_hat_sc1, 1, quantile, probs=c(0.975), na.rm=T) #
      dat$uncert_pop_hat_sc1 <- (dat$upper_pop_hat_sc1 - dat$lower_pop_hat_sc1)/dat$mean_pop_hat_sc1


      output <- list(pop_hat = pop_hat,
                     pop_hat_sc1 = pop_hat_sc1,
                     est_data = dat)

    }


    run=100#--Note that run=100 or above would suffice but the larger the better

    # Scale factor gridding
    #--add province random effects for scale predictions
    dat_pred_sc <- dat_pred
    dat_pred_sc$prov_ranef <- rand_eff(dat_pred_sc, msc$summary.random$prv_id$mean,
                                       dat_pred_sc$prv_id)
    dat_pred_sc$clust_ranef <- rand_eff(dat_pred_sc, msc$summary.random$clst_id$mean,
                                        dat_pred_sc$clst_id)
    sim.scale <- sim_scale(msc, dat_pred_sc, A_pred_sc, run)
    scale1 <- sim.scale$scale_hat


    #dat_pred$ID_ranef <- mc3$summary.random$ID$mean
    dat_pred$prov_ranef <- rand_eff(dat_pred, mc3$summary.random$prv_id$mean,
                                    dat_pred$prv_id)
    dat_pred$clust_ranef <- rand_eff(dat_pred, mc3$summary.random$clst_id$mean,
                                     dat_pred$clst_id)
    ##---add scale factors to the pop grid data
    dat_pred$scale1 <- sim.scale$est_data$mean_scale_hat #-- add the predicted scale factor to the grid data

    #------------------------------------------------------------------------------
    system.time(str(sim.sam <-  post_dens(mc3,dat_pred, A_pred, scale1, run)))
    sum(sim.sam$est_data$mean_pop_hat, na.rm=T)
    sum(sim.sam$est_data$mean_pop_hat_sc1, na.rm=T)

    ##---------------------------------------------------------------------------------
    dat3 <- sim.sam$est_data


    model_metrics <- function(obs, pred, upper, lower)
    {
      residual = pred - obs
      MAE = mean(abs(residual), na.rm=T)
      MSE = mean(residual^2, na.rm=T)
      MAPE = 100*sum(abs((obs-pred)/obs))/length(obs)
      RMSE = sqrt(MSE)
     # In_IC = mean(obs<upper & obs> lower, na.rm=T)*100
      corr = cor(obs[!is.na(obs)],pred[!is.na(obs)])

      output <- list(MAE  = MAE,
                     RMSE = RMSE,
                     CC = corr,
                     MAPE = MAPE)
      return(output)
    }

    ###----Model cross validation and fit metrics




    #  Unscaled
    met <- model_metrics(dat_pred$pop, sim.sam$est_data$mean_pop_hat, sim.sam$est_data$upper_pop_hat,
                         sim.sam$est_data$lower_pop_hat)
    ## - scaled
    metsc <- model_metrics(dat_pred$pop, sim.sam$est_data$mean_pop_hat_sc1, sim.sam$est_data$upper_pop_hat_sc1,
                           sim.sam$est_data$lower_pop_hat_sc1)

    ##
    met_dat <- data.frame(rbind(unlist(met), unlist(metsc)))
    rownames(met_dat) <- c("Base", "Proposed")
    print(met_dat)
    write.csv(met_dat, file=paste0(result_path2,"/fit_metrics.csv"))
#

  }
}

path <- "//worldpop.files.soton.ac.uk/Worldpop/Projects/WP517763_GRID3_Phase2/Working/DRC/CHRIS_N/province"
out_path <- paste0(path , "/Haut-Lomami/paper/simstudy/output_main")



pop.bias <- seq(0.1,1, 0.1) # sample of the population with biased data
mag.bias <- seq(0.1,1, 0.1) # magnitude of bias in the biased population


metric <- c("mae", "rmse", "mape")
method <- c("unscaled", "scaled")

n.pop.bias <- length(pop.bias)
n.mag.bias <- length(mag.bias)
n.metric <- length(metric)
n.method <- length(method)

##---build the dataframe for the metrics
n.size <- n.pop.bias*n.mag.bias*n.metric*n.method
dim(dat.met <- data.frame(expand.grid(method=method,
                                      mag.bias=mag.bias,
                                      pop.bias=pop.bias,
                                      metric=metric)))

###

#dat_met <- matrix(0, nrow=n.size, ncol=5)
dat_met1 <- list()
dat_met2 <- list()
for(j in 1:n.pop.bias)
{
  pathp <- paste0(out_path,"/outputs_for_", pop.bias[j]*100,"%","_biased_sample")
  for(k in 1:n.mag.bias)
  {
    pathb <- paste0(pathp,"/", mag.bias[k]*100, "%","_bias")
    met0 <- read.csv(paste0(pathb, "/fit_metrics.csv"))
    met0$mag_bias <- rep(mag.bias[k]*100, 2) # biased population proportion
    met0$pop_bias <- rep(pop.bias[j]*100, 2) # magnitude of bias
    #met0 <- c(pop.bias[j], met0)
    dat_met1[[k]] = met0
  }
  dat_met2[[j]] = dat_met1
}
#dat_met2
unnest_list <- unlist(dat_met2, recursive = FALSE)  #--unnest the list
metrics <- do.call(rbind, unnest_list)
metrics <- metrics %>% rename(method="X") # rename the variable 'X' to model

# pop_bias is the proportion of the entire population with biased values
# mag_bias is the magnitude of bias in each biased sample

setwd(out_path)
write.csv(metrics, "combined_fit_metrics.csv", row.names=FALSE)


# Convert to long format for plotting
# install.packages("reshape2")
library(reshape2)
dim(met_long <- melt(metrics, id.vars=c("method","pop_bias","mag_bias"),
                     value.name="estimate", variable.name = "metric"))

met_long$method = factor(met_long$method)
met_long$pop_bias2 = factor(met_long$pop_bias)
write.csv(met_long, "combined_fit_metrics_long.csv", row.names=FALSE)

library(ggpubr)
##----- MAE
dim(mae <- met_long[met_long$metric=="MAE",])
#plot_mae <- ggplot(mae, aes(x=bld_cover, y=estimate, color = pop_cover))+
#geom_point()+
#geom_line()+
#facet_wrap(~method)

plot_mae <- ggline(mae, x = "mag_bias", y = "estimate",
                   error.plot = "estimate",
                   facet.by = "method",
                   panel.labs= list(method=c("Base", "Proposed")),
                   color = "pop_bias2",
                   palette = "aaas",
                   point.size=1.2,
                   ggtheme = theme_bw(),
                   #linetype = "pop_bias",
                   size=1)
rmae <- ggpar(plot_mae,
              main = "",
              #font.main = c(14,"bold.italic", "red"),
              legend = "top",
              legend.title="Proportion of \n biased population (%)",
              font.x = c(18),
              font.y = c(18),
              font.legend = c(18),
              font.label=list(size=18),
              font.xtickslab =c(18),
              font.ytickslab =c(18),
              xlab="Magnitude of bias(%)",
              ggtheme = theme(strip.text.x = element_text(size = 20)),

              ylab = "MAE",
              xtickslab.rt = 45, ytickslab.rt = 45)
rmae






##----- RMSE
dim(rmse <- met_long[met_long$metric=="RMSE",])
#plot_rmse <- ggplot(rmse, aes(x=bld_cover, y=estimate, color = pop_cover))+
#geom_point()+
#geom_line()+
#facet_wrap(~method)

plot_rmse <- ggline(rmse, x = "mag_bias", y = "estimate",
                    error.plot = "estimate",
                    facet.by = "method",
                    panel.labs= list(method=c("Base", "Proposed")),
                    color = "pop_bias2",
                    palette = "aaas",
                    #palette = c("#00AFBB", "#E7B800"))
                    point.size=1.2,
                    #linetype = "pop_bias",
                    ggtheme = theme_bw(),
                    #ggtheme = theme(strip.text.x = element_text(size = 20)),
                    size=1)

rrmse <- ggpar(plot_rmse,
               main = "",
               legend = "top", legend.title=element_text("Proportion of \n Biased sample(%)"),
               font.legend=c(18),
               # palette = c("aaa"),
               font.label = list(size = 18, face = "bold", color ="red"),
               font.x = c(18),
               font.y = c(18),
               font.main=c(18),
               font.xtickslab =c(16),
               font.ytickslab =c(18),
               xtickslab.rt = 45, ytickslab.rt = 45,

               ggtheme = theme(strip.text.x = element_text(size = 20)),
               xlab = "Magnitude of bias(%)",
               ylab = "RMSE")
rrmse






##----- MAPE
dim(mape <- met_long[met_long$metric=="MAPE",])
mape$mag_bias2 <- factor(mape$mag_bias)
mape$method2 <- factor(mape$method)
plot_mape <- ggdotplot(mape, x = "mag_bias2", y = "estimate",
                   # error.plot = "estimate",
                    facet.by = "method",
                    panel.labs= list(method=c("Base", "Proposed")),
                    fill = "pop_bias2",
                    palette = "aaas",
                    #palette = c("#00AFBB", "#E7B800"))
                    point.size=0.5,
                    #linetype = "pop_bias",
                    ggtheme = theme_bw(),
                    #ggtheme = theme(strip.text.x = element_text(size = 20)),
                    size=0.8)

rmape <- ggpar(plot_mape,
               main = "",
               legend = "top", legend.title=element_text("Proportion of \n Biased sample(%)"),
               font.legend=c(18),
                palette = c("aaa"),
               font.label = list(size = 18, face = "bold", color ="red"),
               font.x = c(18),
               font.y = c(18),
               font.main=c(18),
               font.xtickslab =c(16),
               font.ytickslab =c(18),
               xtickslab.rt = 45, ytickslab.rt = 45,
               
               ggtheme = theme(strip.text.x = element_text(size = 20)),
               xlab = "Magnitude of bias(%)",
               ylab = "MAPE(%)")
rmape


ggdotchart(cc, x = "mag_bias2", y = "estimate",
           color = "pop_bias2", size = 3,
           add = "segment",
           facet.by = "method",
           add.params = list(color = "lightgray", size = 1.5),
           position = position_dodge(0.3),
           palette = "jco",
           ggtheme = theme_pubclean()
)



##----- CC
dim(cc <- met_long[met_long$metric=="CC",])
cc$mag_bias2 <- factor(mape$mag_bias)
cc$method2 <- factor(mape$method)
cc$pop_bias2 <- factor(mape$pop_bias)
plot_cc <- ggdotchart(cc, x = "mag_bias", y = "estimate",
                 #error.plot = "estimate",
                  facet.by = "method",
                  panel.labs= list(method=c("Base", "Proposed")),
                  color = "pop_bias2",
                  palette = "aaas",
                  #palette = c("#00AFBB", "#E7B800"))
                  point.size=1.8,
                  #linetype = "pop_bias",
                  ggtheme = theme_bw(),
                  size=2)
rcc <- ggpar(plot_cc,
               main = "",
               legend = "top", legend.title=element_text("Proportion of \n Biased sample(%)"),
               font.legend=c(18),
               palette = c("aaa"),
               font.label = list(size = 18, face = "bold", color ="red"),
               font.x = c(18),
               font.y = c(18),
               font.main=c(18),
               font.xtickslab =c(16),
               font.ytickslab =c(18),
               xtickslab.rt = 45, ytickslab.rt = 45,
               
               ggtheme = theme(strip.text.x = element_text(size = 20)),
               xlab = "Magnitude of bias(%)",
               ylab = "Correlation Coefficient (CC)")
rcc


ggarrange(rrmse, rmape,nrow = 1, ncol=2)

##   Plot all together
##---Plots
variable_names <- list(
  "mae" = "MAE" ,
  "rmse" = "RMSE"
)
variable_names["mae"]


met_long$metric2 <- met_long$metric2
levels(met_long$method) <- c("Scaled","Unscaled")
grp <- levels(met_long$method)

variable_labeller2 <- function(variable,value){
  if (variable=='metric') {
    return(variable_names[value])
  } else {
    return(grp)
  }
}



plot_metrics <- ggplot(met_long, aes(x=mag_bias, y=estimate))+
  geom_point(aes(colour=pop_bias2), size=2)+
  geom_line(aes(colour=pop_bias2), size=1)+
  theme_bw()+
  theme(strip.text = element_text(size = 15),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=14))+
  facet_grid(metric~method, scales="free",
             space="free_x",  labeller= variable_labeller2)
plot_metrics

plot_met <- ggpar(plot_metrics, xlab="Magnitude of bias (%)", ylab="Estimate",
                  legend = "right", legend.title = "Proportion of biased \n population (%)",
                  palette = "jco",
                  font.label = list(size = 15, face = "bold", color ="red"),
                  font.x = c(16),
                  font.y = c(16),
                  xtickslab.rt = 45, ytickslab.rt = 45)

plot_met
##





### Reduction in relative bias
rrmae_dat <- as.data.frame(read.csv(paste0(out_path, "/reduction_in_rrmae.csv")))
rrmae_dat <- rrmae_dat %>% drop_na(rrmae)%>% dplyr::select(c(mag_bias,
                                                                pop_bias,
                                                                rrmae))


rrmae_dat$pop_bias <- factor(rrmae_dat$pop_bias)
rrmae_dat$mag_bias <- factor(rrmae_dat$mag_bias,
                              labels=c("Bias:10%", "Bias:20%",
                                       "Bias:30%", "Bias:40%",
                                       "Bias:50%", "Bias:60%",
                                       "Bias:70%", "Bias:80%",
                                       "Bias:90%", "Bias:100%"))

levels(rrmae_dat$mag_bias)
barp <- rrmae_dat %>%
  ggplot( aes(x=pop_bias, y=rrmae, fill=mag_bias)) +
  geom_bar(stat = "identity",
           position = position_dodge()) +
  #scale_x_continuous(breaks = seq(20,90, by=10))+
  #scale_y_continuous(breaks=seq(500,7000,1000))+
  coord_flip()+
  theme_bw()+
  theme(strip.text = element_text(size = 16),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=14))+
  facet_wrap(~mag_bias)
#


rcval5 <-  ggpar(barp, xlab="Proportion of biased population(%)",
                 ylab="Reduction in relative MAE(%)",
                 legend = "none", legend.title=element_text("Method:"),
                 font.legend=c(20),
                 palette = c("aaas"),
                 font.label = list(size = 16, face = "bold", color ="red"),
                 font.x = c(16),
                 font.y = c(16),
                 font.main=c(16),
                 font.xtickslab =c(16),
                 font.ytickslab =c(12),
                 xtickslab.rt = 45, ytickslab.rt = 45)
rcval5


min(rrmae_dat$rrmae) # 0.17%
max(rrmae_dat$rrmae) # 51.66%





#-------------------------------------------------------------
#     undercount
#----------------------------------------------------------
path <- "//worldpop.files.soton.ac.uk/Worldpop/Projects/WP517763_GRID3_Phase2/Working/DRC/CHRIS_N/province"
data_path <- paste0(path, "/data/")#---paths for survey data:
out_path <- paste0(path , "/Haut-Lomami/paper/simstudy/output_main2")
coverp <- c(0.1, 0.2, 0.3, 0.4, 0.5,0.6, 0.7,0.8,0.9,1)
coverb <- c(0.1, 0.2, 0.3, 0.4, 0.5,0.6, 0.7,0.8,0.9)

for(i in 1:length(coverp))
{
  
  # i =1
  result_path1 <- paste0(out_path,"/outputs_for_", coverp[i]*100,"%","_biased_sample")
  if (file.exists(result_path1)){
   # setwd(file.path(result_path1))
  } else {
    dir.create(file.path(result_path1))
    #setwd(file.path(result_path1))
  }
  
  pm <- coverp[i]
  print(pm*100)#--
  
  for(j in 1:length(coverb))
  {
    
    # j=5
    bm = coverb[j]
    result_path2 <- paste0(result_path1,"/", coverb[j]*100, "%","_bias")
    if (file.exists(result_path2)){
      #setwd(file.path(result_path2))
    } else {
      dir.create(file.path(result_path2))
      #setwd(file.path(result_path2))
    }
    print(c(pm*100,bm*100))#---
    
    ##--Aggregate Data to cluster level
    # Calculate/extract the cluster level variables
    dat.grid <- ddat
    dim(dat.clust <- data.frame(dat.grid %>%
                                  group_by(clst_id,prv_id) %>%
                                  dplyr::summarise(bld = sum(bld),
                                            pop = sum(pop),
                                            x1 = mean(x1),
                                            x2 = mean(x2),
                                            x3 = mean(x3),
                                            x4 = mean(x4),
                                            x5 = mean(x5),
                                            x6 = mean(x6))))
    
    
    # Apply missing values and inflated values to the cluster data
    dat.clust$popm <- dat.clust$pop
    ind.obs <- sample(nrow(dat.clust), pm*nrow(dat.clust))
    
    dat.clust$popm[ind.obs] = dat.clust$pop[ind.obs] - round(dat.clust$pop[ind.obs]*bm)
    
    # Add scale variable to the cluster data
    dat.clust$scale <- dat.clust$pop/dat.clust$popm
    
    
    ##---Add lon - lat to cluster data
    dat.clust$lon <- coords.clust[,1]
    dat.clust$lat <- coords.clust[,2]
    with(dat.clust, plot(lon, lat, pch=16, col="red"))
    coords.clust <- cbind(dat.clust$lon, dat.clust$lat)
    
    ###--Cluster level mesh, spde and projection matrix
    bnd_clust <- inla.nonconvex.hull(coords.clust, -0.035, -0.05, resolution = c(100, 100))
    mesh.clst <- inla.mesh.2d(boundary =bnd_clust, max.edge=c(0.01,1),
                              offset = c(0.01, 0.15),
                              cutoff = 0.0309)
    
    par(mfrow=c(1,1))
    plot(mesh.clst)
    with(dat.clust, points(lon, lat, pch=16, col="red"))
    plot(prov, lwd=2, add=T)
    mesh.clst$n
    plot(prov, lwd=2)
    ###---Build projector matrix A for clusters
    A.clst<-inla.spde.make.A(mesh=mesh.clst,loc=as.matrix(coords.clust));dim(A.clst)
    
    ##---Create the SPDE for clusters
    spde.clst <- inla.spde2.matern(mesh.clst, alpha=2)
    
    ##----specify the observation indices for estimation  for clusters
    iset.clst <- inla.spde.make.index(name = "s", spde.clst$n.spde)
    
    #@@@@@@@@@@------------------------------------------------------------------------------
    library(car) ##--For calculating variance inflation factor (vif)
    #dat.fit.grd <- dat.grid
    dat.fit.clst <- dat.clust
    
    ## Covariates and data stacking for both grid and cluster
    #dat.fit.grd[, c("x2", "x3", "x4", "x5", "x6")] <- apply(dat.grid[, c("x2", "x3", "x4", "x5", "x6")], 2, stdize2)
    dat.fit.clst[, c("x2", "x3", "x4", "x5", "x6")] <- apply(dat.fit.clst[, c("x2",
                                                                              "x3", "x4", "x5", "x6")], 2, stdize2)
    
    #---
    #names(dat.fit.grd);
    names(dat.fit.clst)
    dat.fit.clst$prv_id2 <- factor(dat.fit.clst$prv_id)
    dat.fit.clst$ID <- 1:nrow(dat.fit.clst)
    covars.fit.clst <- dat.fit.clst[,c("ID","x2", "x3", "x4", "x5", "x6", "bld", "clst_id",
                                       "prv_id", "prv_id2")]; dim(covars.fit.clst) ##---Population density
    
    
    #---Build the stacks-----------------------------------
    ##-- for scale  at cluster leve;
    stk.clst.sc <- inla.stack(data=list(y=dat.fit.clst$scale), #the response
                              
                              A=list(A.clst,1),  #the A matrix; the 1 is included to make the list(covariates)
                              
                              effects=list(c(list(Intercept=1), #the Intercept
                                             iset.clst),  #the spatial index
                                           #the covariates
                                           list(covars.fit.clst)
                              ),
                              #this is a quick name so you can call upon easily
                              tag='est')
    
    
    ##-- for density at cluster level
    dat.fit.clst$pop2 <- dat.fit.clst$popm
    dat.fit.clst$pop2[dat.fit.clst$pop2==0]= NA
    stk.clst.dens <- inla.stack(data=list(y=dat.fit.clst$pop2/dat.fit.clst$bld), #the response
                                
                                A=list(A.clst,1),  #the A matrix; the 1 is included to make the list(covariates)
                                
                                effects=list(c(list(Intercept=1), #the Intercept
                                               iset.clst),  #the spatial index
                                             #the covariates
                                             list(covars.fit.clst)
                                ),
                                #this is a quick name so you can call upon easily
                                tag='est')
    
    
    ##---Quasi Poisson Models
    ######-Scale---------------------------------------------
    
    #--Gamma
    fsc <- y ~ -1 + Intercept + x2 + x3 + x4 + x5 + x6 + f(s, model=spde.clst) +
      f(prv_id, model='iid') +
      f(clst_id, model='iid')
    msc <-inla(fsc, #the formula
               data=inla.stack.data(stk.clst.sc,spde=spde.clst),  #the data stack
               family= 'gamma',   #which family the data comes from
               control.predictor=list(A=inla.stack.A(stk.clst.sc),compute=TRUE),  #compute gives you the marginals of the linear predictor
               control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
               verbose = FALSE)
    summary(msc)
    
    ind.clst.sc <-inla.stack.index(stk.clst.sc, "est")$data
    fitsc <- exp(msc$summary.linear.predictor[ind.clst.sc,"mean"])
    
    
    ##
    fc3 <- y ~ -1 + Intercept + x2 + x3 + x4 + x5 + x6 + f(s, model=spde.clst) +
      f(prv_id, model='iid') +
      f(clst_id, model='iid')
    mc3 <-inla(fc3, #the formula
               data=inla.stack.data(stk.clst.dens,spde=spde.clst),  #the data stack
               family= 'gamma',   #which family the data comes from
               control.predictor=list(A=inla.stack.A(stk.clst.dens),compute=TRUE),  #compute gives you the marginals of the linear predictor
               control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
               verbose = FALSE)
    summary(mc3)
    
    
    
    #############################################################
    ##-----Posterior Explorations/Simulations-----------------
    ##-----Extract province level random effects----------------
    ##---Province level posterior random effects
    # Scale grid covs
    dat_pred <- dat.grid
    dat_pred[, c("x2", "x3", "x4", "x5", "x6")] <- apply(dat_pred[, c("x2",
                                                                      "x3", "x4", "x5", "x6")], 2, stdize2)
    
    # Grid predictions for scale factor
    A_pred_sc <- inla.spde.make.A(mesh = mesh.clst, loc = cbind(dat_pred$lon, dat_pred$lat))
    dim(A_pred_sc)
    
    sim_scale <- function(model, dat, Aprediction, run)
    {
      fixedeff  <- scale_hat <- matrix(0, nrow=nrow(dat), ncol = run)
      inla.seed = as.integer(runif(1)*.Machine$integer.max)
      #inla.seed =  1657687559
      set.seed(inla.seed)
      print(inla.seed)
      m1.samp <- inla.posterior.sample(run, model,
                                       seed = inla.seed,
                                       selection=list(x2=1,
                                                      x3=1,
                                                      x4=1,
                                                      x5=1,
                                                      x6=1),num.threads="1:1")
      
      sfield_nodes_mean <- model$summary.random$s['mean']
      field_mean <- (Aprediction%*% as.data.frame(sfield_nodes_mean)[, 1])
      for(i in 1:run)
      {
        #
        fixedeff[,i] <-
          model$summary.fixed['Intercept', 'mean'] +
          m1.samp[[i]]$latent[1,] * dat[,'x2'] +
          m1.samp[[i]]$latent[2,] * dat[,'x3'] +
          m1.samp[[i]]$latent[3,] * dat[,'x4'] +
          m1.samp[[i]]$latent[4,] * dat[,'x5'] +
          m1.samp[[i]]$latent[5,] * dat[,'x6'] +
          dat$clust_ranef +
          dat$prov_ranef +
          field_mean[,1]
        
        scale_hat[,i]<- exp(fixedeff[,i])
      }
      
      dat$mean_scale_hat <- apply(scale_hat, 1, mean, na.rm=T) #
      dat$lower_scale_hat <- apply(scale_hat, 1, quantile, probs=c(0.025), na.rm=T) #
      dat$upper_scale_hat <- apply(scale_hat, 1, quantile, probs=c(0.975), na.rm=T) #
      dat$uncert_scale_hat <- (dat$upper_scale_hat - dat$lower_scale_hat)/dat$mean_scale_hat#
      
      
      output <- list(scale_hat = scale_hat,
                     est_data = dat)
      
    }
    
    
    # Grid Predictions for cluster populations
    
    A_pred <- inla.spde.make.A(mesh = mesh.clst, loc = cbind(dat_pred$lon, dat_pred$lat))
    dim(A_pred)
    
    ##
    #-for density
    post_dens <- function(model, dat, Aprediction, scale, run)
    {
      fixedeff <- dens_hat <- pop_hat <- pop_hat_sc1 <- matrix(0, nrow=nrow(dat), ncol = run)
      inla.seed = as.integer(runif(1)*.Machine$integer.max)
      #inla.seed =
      
      set.seed(inla.seed)
      print(inla.seed)
      m1.samp <- inla.posterior.sample(run, model,
                                       seed = inla.seed,
                                       selection=list(x2=1,
                                                      x3=1,
                                                      x4=1,
                                                      x5=1,
                                                      x6=1),num.threads="1:1")
      
      sfield_nodes_mean <- model$summary.random$s['mean']
      field_mean <- (Aprediction%*% as.data.frame(sfield_nodes_mean)[, 1])
      for(i in 1:run)
      {
        fixedeff[,i] <-
          model$summary.fixed['Intercept', 'mean'] +
          m1.samp[[i]]$latent[1,] * dat[,'x2'] +
          m1.samp[[i]]$latent[2,] * dat[,'x3'] +
          m1.samp[[i]]$latent[3,] * dat[,'x4'] +
          m1.samp[[i]]$latent[4,] * dat[,'x5'] +
          m1.samp[[i]]$latent[5,] * dat[,'x6'] +
          dat$clust_ranef +
          dat$prov_ranef +
          field_mean[,1]
        
        ##
        dens_hat[,i] <- exp(fixedeff[,i]) #--density
        pop_hat[,i] <- dens_hat[,i]*dat$bld #--population count
        #pop_hat_sc1[,i] <- pop_hat[,i]*dat$scale1 #--national scale factor
      }
      
      #--Unscaled
      dat$mean_dens_hat <- apply(dens_hat, 1, mean, na.rm=T) #
      dat$mean_pop_hat <- apply(pop_hat, 1, mean, na.rm=T) #
      dat$lower_pop_hat <- apply(pop_hat, 1, quantile, probs=c(0.025), na.rm=T) #
      dat$upper_pop_hat <- apply(pop_hat, 1, quantile, probs=c(0.975), na.rm=T) #
      dat$uncert_pop_hat <- (dat$upper_pop_hat - dat$lower_pop_hat)/dat$mean_pop_hat#
      
      
      #--Scale 1 #--Scaled
      pop_hat_sc1 <- pop_hat*scale # allows for scale-pop uncertainty propagation
      dat$mean_pop_hat_sc1 <- apply(pop_hat_sc1, 1, mean, na.rm=T) #
      dat$lower_pop_hat_sc1 <- apply(pop_hat_sc1, 1, quantile, probs=c(0.025), na.rm=T) #
      dat$upper_pop_hat_sc1 <- apply(pop_hat_sc1, 1, quantile, probs=c(0.975), na.rm=T) #
      dat$uncert_pop_hat_sc1 <- (dat$upper_pop_hat_sc1 - dat$lower_pop_hat_sc1)/dat$mean_pop_hat_sc1
      
      
      output <- list(pop_hat = pop_hat,
                     pop_hat_sc1 = pop_hat_sc1,
                     est_data = dat)
      
    }
    
    
    run=100#--Note that run=100 or above would suffice but the larger the better
    
    # Scale factor gridding
    #--add province random effects for scale predictions
    dat_pred_sc <- dat_pred
    dat_pred_sc$prov_ranef <- rand_eff(dat_pred_sc, msc$summary.random$prv_id$mean,
                                       dat_pred_sc$prv_id)
    dat_pred_sc$clust_ranef <- rand_eff(dat_pred_sc, msc$summary.random$clst_id$mean,
                                        dat_pred_sc$clst_id)
    sim.scale <- sim_scale(msc, dat_pred_sc, A_pred_sc, run)
    scale1 <- sim.scale$scale_hat
    #sim.scale$scale_hat
    
    #dat_pred$ID_ranef <- mc3$summary.random$ID$mean
    dat_pred$prov_ranef <- rand_eff(dat_pred, mc3$summary.random$prv_id$mean,
                                    dat_pred$prv_id)
    dat_pred$clust_ranef <- rand_eff(dat_pred, mc3$summary.random$clst_id$mean,
                                     dat_pred$clst_id)
    ##---add scale factors to the pop grid data
    dat_pred$scale1 <- sim.scale$est_data$mean_scale_hat #-- add the predicted scale factor to the grid data
    
    #------------------------------------------------------------------------------
    system.time(str(sim.sam <-  post_dens(mc3,dat_pred, A_pred, scale1, run)))
    sum(sim.sam$est_data$mean_pop_hat, na.rm=T)
    sum(sim.sam$est_data$mean_pop_hat_sc1, na.rm=T)
    
    ##---------------------------------------------------------------------------------
    dat3 <- sim.sam$est_data
    # dim(data.sim1<- data.frame(cbind(dat3[,c("grd_id","prv_id","lon","lat","pop","bld")], sim.sam$pop_hat)))
    #dim(data.sim2<- data.frame(cbind(dat3[,c("grd_id","prv_id","lon","lat","pop","bld")], sim.sam$pop_hat_sc1)))
    
    
    model_metrics <- function(obs, pred, upper, lower)
    {
      residual = pred - obs
      MAE = mean(abs(residual), na.rm=T)
      MSE = mean(residual^2, na.rm=T)
      MAPE = 100*sum(abs((obs-pred)/obs))/length(obs)
      RMSE = sqrt(MSE)
      # In_IC = mean(obs<upper & obs> lower, na.rm=T)*100
      #corr = cor(obs[!is.na(obs)],pred[!is.na(obs)])
      
      output <- list(MAE  = MAE,
                     RMSE = RMSE,
                     MAPE = MAPE)
      return(output)
    }
    
    ###----Model cross validation and fit metrics
    
    
    
    
    #  Unscaled
    met <- model_metrics(dat_pred$pop, sim.sam$est_data$mean_pop_hat, sim.sam$est_data$upper_pop_hat,
                         sim.sam$est_data$lower_pop_hat)
    ## - scaled
    metsc <- model_metrics(dat_pred$pop, sim.sam$est_data$mean_pop_hat_sc1, sim.sam$est_data$upper_pop_hat_sc1,
                           sim.sam$est_data$lower_pop_hat_sc1)
    
    ##
    met_dat <- data.frame(rbind(unlist(met), unlist(metsc)))
    rownames(met_dat) <- c("Base", "Proposed")
    print(met_dat)
    write.csv(met_dat, file=paste0(result_path2,"/fit_metrics.csv"))
    #
    
  }
}

path <- "//worldpop.files.soton.ac.uk/Worldpop/Projects/WP517763_GRID3_Phase2/Working/DRC/CHRIS_N/province"
out_path <- paste0(path , "/Haut-Lomami/paper/simstudy/output_main2")



pop.bias <- seq(0.1,1, 0.1) # sample of the population with biased data
mag.bias <- seq(0.1,0.9, 0.1) # magnitude of bias in the biased population


metric <- c("mae", "rmse", "mape")
method <- c("Base", "Proposed")

n.pop.bias <- length(pop.bias)
n.mag.bias <- length(mag.bias)
n.metric <- length(metric)
n.method <- length(method)

##---build the dataframe for the metrics
n.size <- n.pop.bias*n.mag.bias*n.metric*n.method
dim(dat.met <- data.frame(expand.grid(method=method,
                                      mag.bias=mag.bias,
                                      pop.bias=pop.bias,
                                      metric=metric)))

###

#dat_met <- matrix(0, nrow=n.size, ncol=5)
dat_met1 <- list()
dat_met2 <- list()
for(j in 1:n.pop.bias)
{
  pathp <- paste0(out_path,"/outputs_for_", pop.bias[j]*100,"%","_biased_sample")
  for(k in 1:n.mag.bias)
  {
    pathb <- paste0(pathp,"/", mag.bias[k]*100, "%","_bias")
    met0 <- read.csv(paste0(pathb, "/fit_metrics.csv"))
    met0$mag_bias <- rep(mag.bias[k]*100, 2) # biased population proportion
    met0$pop_bias <- rep(pop.bias[j]*100, 2) # magnitude of bias
    #met0 <- c(pop.bias[j], met0)
    dat_met1[[k]] = met0
  }
  dat_met2[[j]] = dat_met1
}
#dat_met2
unnest_list <- unlist(dat_met2, recursive = FALSE)  #--unnest the list
metrics <- do.call(rbind, unnest_list)
metrics <- metrics %>% rename(method="X") # rename the variable 'X' to model

# pop_bias is the proportion of the entire population with biased values
# mag_bias is the magnitude of bias in each biased sample

setwd(out_path)
write.csv(metrics, "combined_fit_metrics.csv", row.names=FALSE)


# Convert to long format for plotting
# install.packages("reshape2")
metrics <- metrics[metrics$mag_bias < 70,] # view results for up to 60% bias
library(reshape2)
dim(met_long <- melt(metrics, id.vars=c("method","pop_bias","mag_bias"),
                     value.name="estimate", variable.name = "metric"))

met_long$method = factor(met_long$method)
met_long$pop_bias2 = factor(met_long$pop_bias)
write.csv(met_long, "combined_fit_metrics_long.csv", row.names=FALSE)

library(ggpubr)
##----- MAE
dim(mae <- met_long[met_long$metric=="MAE",])
#plot_mae <- ggplot(mae, aes(x=bld_cover, y=estimate, color = pop_cover))+
#geom_point()+
#geom_line()+
#facet_wrap(~method)

plot_mae <- ggline(mae, x = "mag_bias", y = "estimate",
                   error.plot = "estimate",
                   facet.by = "method",
                   panel.labs= list(method=c("Base", "Proposed")),
                   color = "pop_bias2",
                   palette = "aaas",
                   point.size=1.2,
                   ggtheme = theme_bw(),
                   #linetype = "pop_bias",
                   size=1)
rmae <- ggpar(plot_mae,
              main = "",
              #font.main = c(14,"bold.italic", "red"),
              legend = "top",
              legend.title="Proportion of \n biased population (%)",
              font.x = c(18),
              font.y = c(18),
              font.legend = c(18),
              font.label=list(size=18),
              font.xtickslab =c(18),
              font.ytickslab =c(18),
              xlab="Magnitude of bias(%)",
              ggtheme = theme(strip.text.x = element_text(size = 20)),
              
              ylab = "MAE",
              xtickslab.rt = 45, ytickslab.rt = 45)
rmae




##----- RMSE
dim(rmse <- met_long[met_long$metric=="RMSE",])
#plot_rmse <- ggplot(rmse, aes(x=bld_cover, y=estimate, color = pop_cover))+
#geom_point()+
#geom_line()+
#facet_wrap(~method)

plot_rmse <- ggline(rmse, x = "mag_bias", y = "estimate",
                    error.plot = "estimate",
                    facet.by = "method",
                    panel.labs= list(method=c("Base", "Proposed")),
                    color = "pop_bias2",
                    palette = "aaas",
                    #palette = c("#00AFBB", "#E7B800"))
                    point.size=1.2,
                    #linetype = "pop_bias",
                    ggtheme = theme_bw(),
                    #ggtheme = theme(strip.text.x = element_text(size = 20)),
                    size=1)

rrmse <- ggpar(plot_rmse,
               main = "",
               legend = "top", legend.title=element_text("Proportion of \n Biased sample(%)"),
               font.legend=c(18),
               # palette = c("aaa"),
               font.label = list(size = 18, face = "bold", color ="red"),
               font.x = c(18),
               font.y = c(18),
               font.main=c(18),
               font.xtickslab =c(16),
               font.ytickslab =c(18),
               xtickslab.rt = 45, ytickslab.rt = 45,
               
               ggtheme = theme(strip.text.x = element_text(size = 20)),
               xlab = "Magnitude of bias(%)",
               ylab = "RMSE")
rrmse
ggarrange(rmae, rrmse,nrow = 2, ncol=1)

save.image(paste0(out_path,"/HL_sim_study.Rdata"))

