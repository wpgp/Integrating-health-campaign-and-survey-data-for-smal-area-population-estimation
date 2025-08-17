#-----------------------------------------------------------------------------
#  Ratio-based Bayesian hierarchical population modelling - DRC Application
#  Developer: Dr Chris Nnanatu, WOrldPop, University of Southampton
# February, 2025
#----------------------------------------------------------------------------
rm(list=ls()) #----Clear the workspace

# install some key libraries
packages <- c("raster", "haven", "sf","sp", "tmap","tmaptools","tidyverse","terra",
              "lattice", "gridExtra", "devtools", "rlang", "viridis", "spdep",
              "car", "MASS", "maps", "spData", "spDataLarge", "ceramic", "basemaps", "ggmap")
if(length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages()))) }


#Install INLA
if(length(setdiff("INLA", rownames(installed.packages()))) > 0){
  install.packages("INLA", type="binary", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
}
library(INLA)
lapply(packages, library, character.only = TRUE) ##--access the libraries


# Set working directory for input and output files
path <- "//worldpop.files.soton.ac.uk/Worldpop/Projects/WP517763_GRID3_Phase2/Working/DRC/CHRIS_N/province/Haut-Lomami/paper/application"
data_path <- paste0(path, "/data/")#---paths for survey data:
results_path <- paste0(path , "/output")


## Key steps: 
# 1. Data preparation/cleaning
# 2. Standardize continuous covariates 
# 3. Fit the Bayesian models with the best fit covariates
# 4. Choose the best Bayesian models based on the DIC values
# 5. Predict both the adjustment factor and population density at grid cells
# 6. Extract the posterior results and write the raster files
# 7. Run model cross-validation based on the best Bayesian model
# 8. Extract and map other desired aggregated estimates 

#-------------------------------------------------
#### 1. Data preparation/cleaning
#-------------------------------------------------
# Load the datasets
# Data with overlapping clusters across 6 provinces in the DRC
dat_sc1 <- read.csv(paste0(data_path,"Scaling_Factor.csv"))# Cluster level data across 6 DRC provinces
dim(dat_sc1)
shp_sc1 <- st_read(paste0(data_path, "Scaling_Factor.shp")) # shapefile across the 6 overlapping provinces

# DRC province-level shapefile
drc_shp <- st_read(paste0(data_path, "DRC_Province_from_HZ_20220921.shp"))
dim(drc_shp)
plot(drc_shp["Shape_Area"]) # quick view

# Read in health area level shapefile
harea_hl <- st_read(paste0(data_path, "/Haut_Lomami_Health_Area.shp")) # need updated HA boundaries
dim(harea_hl)

# Define the adjustment factor (scale) as the ratio of Microcensus:PNLP observations
dat_sc1$scale <- dat_sc1$Microcensu/dat_sc1$PNLP_Pop
dat_sc1$scale[is.infinite(dat_sc1$scale)] = NA # set infinite values t NA
dat_sc1$scale[is.nan(dat_sc1$scale)] = NA # Set unusual values to NA
dat_sc1$scale[dat_sc1$scale==0] = NA # set 0 values to NA
length(dat_sc1$scale[!is.na(dat_sc1$scale)])#1207 clusters

# Explore the values of the scale factor
boxplot(dat_sc1$scale, cex.axis=1.5)
min(dat_sc1$scale, na.rm=T) # ~0.039
max(dat_sc1$scale, na.rm=T) # 435

# Visualise the boxplots of the scale factor
# include only values of the scale factor (k) between 0.1 and 10 as threshholds
boxplot(dat_sc1$scale[dat_sc1$scale >= 0.1 & dat_sc1$scale  <=10])
dim(dat_sc1[dat_sc1$scale >= 0.1 & dat_sc1$scale  <=10,])

hist(dat_sc1$scale[dat_sc1$scale >= 0.1 & dat_sc1$scale  <=10])
hist(log(dat_sc1$scale[dat_sc1$scale >= 0.1 & dat_sc1$scale  <=10]))
dat_sc1$scale[dat_sc1$scale < 0.1] = NA 
dat_sc1$scale[dat_sc1$scale > 10] = NA
min(dat_sc1$scale, na.rm=T)
boxplot(dat_sc1$scale)
hist(dat_sc1$scale, cex.axis=1.5)
hist(log(dat_sc1$scale))


# Explore the shapefile of the adjustment factor data
# extract the coordinates
shp_sc <- as(st_geometry(shp_sc1), "Spatial")
coords_sc <- as.data.frame(coordinates(shp_sc)) # extract the coordinates
head(coords_sc)

#  Visualise centroids overlaid over shapefiles
plot(shp_sc, col="grey", main="EA shapefile and centroids overlay")
points(coords_sc, col=2, pch="*", cex=1.5)

#  Add the lon-lat to the survey dat_sc1a data frame
dat_sc1$lon <- coords_sc$V1  #---Longitude
dat_sc1$lat <- coords_sc$V2  #--Latitude

# More exploratory checks
dim(dat_sc1)
names(dat_sc1)
barplot(table(dat_sc1$GHSL_SMOD))
zerop <- which(dat_sc1$PNLP_Pop==0)#----Only 1 cluster  with zero pop
zerom <- which(dat_sc1$Microcensus_Pop==0)#----Only 3cluster  with zero pop


dat_sc1$scale[is.infinite(dat_sc1$scale)] = NA #--set infinite values to NA
dat_sc1$scale[is.nan(dat_sc1$scale)] = NA #--set NAN values to NA

#  Visualise the values of the scale factor
hist(dat_sc1$scale, main="")
hist(log(dat_sc1$scale), main="")
boxplot(dat_sc1$scale, main="")

#-------------------------------------------------
### 2. Standardize continous covariates 
#-----------------------------------------------

names(dat_sc1)
dim(dat_sc1)
vars <- paste0("x", 1:85) # there are 85 covariates to be tested
#-------------------------------------

# Write the covariates scaling function for the continuous covariates
stdize <- function(x)
{
  stdz <- (x - mean(x, na.rm=T))/sd(x, na.rm=T)
  return(stdz)
}

# apply scaling to model covariates
dat_sc1[, vars] <- apply(dat_sc1[,vars], 2, stdize)

# Microcensus data for Kongo Central, Kinshasa, Kwilu were collected in 2017/2018 - dropped
# Microcensus data for Haut-Lomami, Sud-Kivu and Ituri were collected in 2021 (M4H)
# subset scale factor for Ituri and SUd-Kivu provinces (M4H 2021) for global scaling
dat_sc1a <- dat_sc1 %>% filter(province %in% c("Ituri", "Sud-Kivu"))
shp_sc1a <- shp_sc1 %>% filter(province %in% c("Ituri", "Sud-Kivu")) # 481 clusters

#--Best fit covariated based on GLM stepwise selection
cov_sc <- c("x61", "x77") # Best covariates selected for the global scaling

## For Population density modelling using Haut-Lomami data 
# subset the adjustment factor data for Haut-Lomami only
dim(dat_sc_hl <- dat_sc1[dat_sc1$province=="Haut-Lomami",]) # 233 clusters
dim(shp_sc_hl <- shp_sc1[shp_sc1$province=="Haut-Lomami",]) # 233 clusters


#-------------------------------------------------
### 3. Fit Bayesian models with the best fit covariates
#-------------------------------------------------
## Note: Global scaling refers to adjustment factor model using Ituri and Sud-Kivu
  #      Local scaling refers to adjustment factor model using Haut-Lomami data

# Extract the coordinates
coords_sc <- cbind(dat_sc1a$lon, dat_sc1a$lat) # Global

coords_sc_hl <- cbind(dat_sc_hl$lon, dat_sc_hl$lat) # Local
out <- which(coords_sc_hl[,2]> -6) # exclude 5 points outside HL
dat_sc_hl <- dat_sc_hl[-out,]
shp_sc_hl <- shp_sc_hl[-out,]
coords_sc_hl <- cbind(dat_sc_hl$lon, dat_sc_hl$lat) # Local
plot(coords_sc_hl)
# Build the mesh

 # Global 
bound_sc <- inla.nonconvex.hull(as.matrix(coords_sc), -0.035, -0.05, resolution = c(100, 100))
mesh_sc <- inla.mesh.2d(boundary = bound_sc, max.edge=c(0.3,5),
                        offset = c(0.5, 1),
                        cutoff = 0.1)

par(mfrow=c(1,1))
plot(mesh_sc)
plot(shp_sc, add=T)
points(coords_sc, col="red")
mesh_sc$n # 601 nodes


# Local 
bound_sc_hl <- inla.nonconvex.hull(as.matrix(coords_sc_hl), -0.0495, -0.06, resolution = c(100, 100))
mesh_sc_hl <- inla.mesh.2d(boundary = bound_sc_hl, max.edge=c(0.3,5),
                        offset = c(0.51, 0.52),
                        cutoff = 0.1)

par(mfrow=c(1,1))
plot(mesh_sc_hl)
plot(shp_sc_hl, add=T)
points(coords_sc_hl, col="red")
mesh_sc_hl$n # 379 nodes


# Global
# Build projector matrix A
A_sc<-inla.spde.make.A(mesh=mesh_sc,loc=as.matrix(coords_sc));dim(A_sc)
# Create the SPDE
spde_sc <- inla.spde2.matern(mesh_sc, alpha=2)
# Specify the observation indices for estimation
iset_sc <- inla.spde.make.index(name = "s", spde_sc$n.spde)

# Local
# Build projector matrix A
A_sc_hl<-inla.spde.make.A(mesh=mesh_sc_hl,loc=as.matrix(coords_sc_hl));dim(A_sc_hl)
# Create the SPDE
spde_sc_hl <- inla.spde2.matern(mesh_sc_hl, alpha=2)
# Specify the observation indices for estimation
iset_sc_hl <- inla.spde.make.index(name = "s", spde_sc_hl$n.spde)

# Recode random components 
  # Global
dat_sc1a$eps <- 1:nrow(dat_sc1a) # IID
dat_sc1a$set_typ <- factor(dat_sc1a$GHSL_SMOD)# settlement type
dat_sc1a$prov <- factor(dat_sc1a$province) # province
boxplot(dat_sc1a$scale)

# Local
dat_sc_hl$eps <- 1:nrow(dat_sc_hl) # IID
dat_sc_hl$set_typ <- factor(dat_sc_hl$GHSL_SMOD)# settlement type
dat_sc_hl$prov <- factor(dat_sc_hl$province) # province
boxplot(dat_sc_hl$scale)

# Model variables for the adjustment factor model
  # Global 
covars_sc <- dat_sc1a[,c(cov_sc,
                         "prov", "set_typ","eps")]; dim(covars_sc)

hist(log(dat_sc1a$scale))
stk_est_sc<- inla.stack(data=list(y=dat_sc1a$scale), #the response
                        
                        A=list(A_sc,1),  #the A matrix; the 1 is included to make the list(covariates)
                        
                        effects=list(c(list(Intercept=1), #the Intercept
                                       iset_sc),  #the spatial index
                                     #the covariates
                                     list(covars_sc)),
                        #this is a quick name so you can call upon easily
                        tag='est')

# Local 
covars_sc_hl <- dat_sc_hl[,c(cov_sc,
                         "prov", "set_typ","eps")]; dim(covars_sc_hl)

hist(log(dat_sc_hl$scale))
stk_est_sc_hl<- inla.stack(data=list(y=dat_sc_hl$scale), #the response
                        
                        A=list(A_sc_hl,1),  #the A matrix; the 1 is included to make the list(covariates)
                        
                        effects=list(c(list(Intercept=1), #the Intercept
                                       iset_sc_hl),  #the spatial index
                                     #the covariates
                                     list(covars_sc_hl)),
                        #this is a quick name so you can call upon easily
                        tag='est')

       # Bayesian Hieracrhical models
###---Mod 1
f1_sc <- y ~ x61 + x77  +  f(s, model=spde_sc)+  f(eps, model="iid")

 # Global
mod1_sc <- inla(f1_sc,
                family="gamma",
                control.compute = list(dic=T, waic=T, cpo=T, config=T),
                data = inla.stack.data(stk_est_sc,spde=spde_sc),
                control.predictor = list(A = inla.stack.A(stk_est_sc), compute=T))
summary(mod1_sc)

# Local
f1_sc_hl <- y ~ x61 + x77  +  f(s, model=spde_sc_hl)+  f(eps, model="iid")
mod1_sc_hl <- inla(f1_sc_hl,
                family="gamma",
                control.compute = list(dic=T, waic=T, cpo=T, config=T),
                data = inla.stack.data(stk_est_sc_hl,spde=spde_sc_hl),
                control.predictor = list(A = inla.stack.A(stk_est_sc_hl), compute=T))
summary(mod1_sc_hl)

###---Mod 2
f2_sc <- y ~x61 + x77 + f(set_typ, model="iid")+  f(eps, model="iid")
# Global
mod2_sc <- inla(f2_sc,
                family="gamma",
                control.compute = list(dic=T, waic=T, cpo=T, config=T),
                data = inla.stack.data(stk_est_sc,spde=spde_sc),
                control.predictor = list(A = inla.stack.A(stk_est_sc), compute=T))
summary(mod2_sc)

# Local
mod2_sc_hl <- inla(f2_sc,
                   family="gamma",
                   control.compute = list(dic=T, waic=T, cpo=T, config=T),
                   data = inla.stack.data(stk_est_sc_hl,spde=spde_sc_hl),
                   control.predictor = list(A = inla.stack.A(stk_est_sc_hl), compute=T))
summary(mod2_sc_hl)



###---Mod 3
f3_sc <- y ~x61 + x77 + f(set_typ, model="iid") + f(s, model=spde_sc)+  f(eps, model="iid")
# Global
mod3_sc <- inla(f3_sc,
                family="gamma",
                control.compute = list(dic=T, waic=T, cpo=T, config=T),
                data = inla.stack.data(stk_est_sc,spde=spde_sc),
                control.predictor = list(A = inla.stack.A(stk_est_sc), compute=T))
summary(mod3_sc)

# Local
f3_sc_hl <- y ~x61 + x77 + f(set_typ, model="iid") + f(s, model=spde_sc_hl)+  f(eps, model="iid")
mod3_sc_hl <- inla(f3_sc_hl,
                   family="gamma",
                   control.compute = list(dic=T, waic=T, cpo=T, config=T),
                   data = inla.stack.data(stk_est_sc_hl,spde=spde_sc_hl),
                   control.predictor = list(A = inla.stack.A(stk_est_sc_hl), compute=T))
summary(mod3_sc_hl)

## Fit the population density model
  # select covariates using stepwise regression
# Rename some key variables
dat_sc_hl <- dat_sc_hl %>% mutate(pop = PNLP_Pop, # PNLP-based population count per cluster
                                  bld = sum.COD_bld_count_20240617_unpublished, # CIESIN's building counts per cluster
                                  area = sum.COD_bld_area_20240617_unpublished, # CIESIN's settled area per cluster
                                  dens1 = PNLP_Pop/sum.COD_bld_count_20240617_unpublished, # people per building
                                  dens2 = PNLP_Pop/sum.COD_bld_area_20240617_unpublished) # people per settled area


cov_dens <- c("x12", "x17", "x55", "x68") # best fit covariates from GLM-based stepwise regression

cov1 = sort(union(cov_dens, cov_sc)) # combine with adjustment factor best covariates

  # Prepare data stack for pop density models
covars_dens <-dat_sc_hl[,c(cov_dens,
                         "set_typ","eps")]; dim(covars_dens)


# Set any zero or infinite density values to NA
dat_sc_hl$dens1[dat_sc_hl$dens1==0] = NA
dat_sc_hl$dens1[is.infinite(dat_sc_hl$dens1)] = NA
#dat_sc_hl$dens1[dat_sc_hl$dens1>20]=NA # deal with outliers
boxplot(dat_sc_hl$dens1) # view

stk_est <- inla.stack(data=list(y=dat_sc_hl$dens1), #the response
                                
                                A=list(A_sc_hl,1),  #the A matrix; the 1 is included to make the list(covariates)
                                
                                effects=list(c(list(Intercept=1), #the Intercept
                                               iset_sc_hl),  #the spatial index
                                             #the covariates
                                             list(covars_dens)),
                                #this is a quick name so you can call upon easily
                                tag='est')


      # Bayesian hierarchical models

# mod1a
f1a <-  y ~ -1 + Intercept + x12 + x17 + x55 + x68 +  f(s, model=spde_sc_hl)+ f(eps, model="iid") 

mod1a <- inla(f1a,
              family="gamma",
              control.compute = list(dic=T, waic=T, cpo=T, config=T),
              data = inla.stack.data(stk_est,spde=spde_sc_hl),
              control.predictor = list(A = inla.stack.A(stk_est), compute=T))

summary(mod1a)
index1a <-inla.stack.index(stk_est, "est")$data #---extract the dat_hla location indices
pred1a <- exp(mod1a$summary.linear.predictor[index1a ,"mean"]) #--predicted mean count
sum(fit1a <- pred1a*dat_sc_hl$bld)


# mod2a
f2a <-  y ~ -1 + Intercept + x12 + x17 + x55 + x68  + f(set_typ, model="iid")+  f(eps, model="iid")

mod2a <- inla(f2a,
              family="gamma",
              control.compute = list(dic=T, waic=T, cpo=T, config=T),
              data = inla.stack.data(stk_est,spde=spde_sc_hl),
              control.predictor = list(A = inla.stack.A(stk_est), compute=T))

summary(mod2a)
index2a <-inla.stack.index(stk_est, "est")$data #---extract the dat_hla location indices
pred2a <- exp(mod2a$summary.linear.predictor[index2a ,"mean"]) #--predicted mean count
sum(fit2a <- pred2a*dat_sc_hl$bld)



# mod3a
f3a <-  y ~ -1 + Intercept + x12 + x17 + x55 + x68  + 
  f(set_typ, model="iid")+  f(s, model=spde_sc_hl)+  f(eps, model="iid")

mod3a <- inla(f3a,
              family="gamma",
              control.compute = list(dic=T, waic=T, cpo=T, config=T),
              data = inla.stack.data(stk_est,spde=spde_sc_hl),
              control.predictor = list(A = inla.stack.A(stk_est), compute=T))

summary(mod3a)
index3a <-inla.stack.index(stk_est, "est")$data #---extract the dat_hla location indices
pred3a <- exp(mod3a$summary.linear.predictor[index3a ,"mean"]) #--predicted mean count
sum(fit3a <- pred3a*dat_sc_hl$bld)


#-------------------------------------------------
### 4. Choose the best Bayesian models based on the DIC values
#-------------------------------------------------
t(c(mod1=mod1_sc$dic$dic,mod2=mod2_sc$dic$dic,mod3=mod3_sc$dic$dic)) # Global scaling
#        mod1      mod2      mod3
#[1,] -519.1153 -500.4568 -478.0946 - mod1

t(c(mod1=mod1_sc_hl$dic$dic,mod2=mod2_sc_hl$dic$dic,mod3=mod3_sc_hl$dic$dic)) # Local scaling
#        mod1     mod2      mod3
#[1,] -258.987 92.77991 -253.5914 - mod1
####---Posterior Explorations

t(c(mod1=mod1a$dic$dic,mod2=mod2a$dic$dic,mod3=mod3a$dic$dic)) # Population density
#        mod1     mod2     mod3
#[1,] 1081.152 1085.836 1080.897 - mod3a

# Add predicted adjustment factor to the data
ind1 <-inla.stack.index(stk_est_sc_hl, "est")$data #--estimation indices
fit1 <- exp(mod1_sc_hl$summary.linear.predictor[ind1,"mean"]) #--extract the backtransformed scale_hat
dat_sc_hl$pred_sc <- fit1


#----------------------------------------------
### 5. Predict both the adjustment factor and population density at grid cells
#-------------------------------------------------
   # Read in the prediction grid
pred_covs_hl <-  readRDS(paste0(data_path,"/Haut_Lomami_covs_stack.rds"))
dim(pred_covs_hl)
names(pred_covs_hl)

pred_covs_hl2 <- data.frame(pred_covs_hl) 
dim(pred_covs_hl2)
head(pred_covs_hl2)
names(pred_covs_hl2)

# Rename coordinates variables to lon - lat
pred_covs_hl2$lon <- pred_covs_hl2$x
pred_covs_hl2$lat <- pred_covs_hl2$y

#plot(pred_covs_hl2$lon,pred_covs_hl2$lat)
#plot(mesh_hl,add=T)


# Check grid covariates with NAs
sapply(pred_covs_hl2[,cov1], function(x) length(x[is.na(x)]))

# Standardize grid covariates
pred_covs_hl2[,cov1] <- as.data.frame(apply(pred_covs_hl2[,cov1],2, stdize))

## Build grid projection matrices
  # for population density
Apred <- inla.spde.make.A(mesh = mesh_sc_hl, loc = cbind(pred_covs_hl2$lon, pred_covs_hl2$lat));dim(Apred)

  # for adjustment factor
    # Global
Apred_sc0 <- inla.spde.make.A(mesh = mesh_sc, loc = cbind(pred_covs_hl2$lon, pred_covs_hl2$lat));dim(Apred_sc0)
    # Local
Apred_sc1 <- inla.spde.make.A(mesh = mesh_sc_hl, loc = cbind(pred_covs_hl2$lon, pred_covs_hl2$lat));dim(Apred_sc1)

# Posterior random effects extraction function
str_ranef <- function(dat, strat, st)
{
  uniq <- unique(strat)
  ranef <- rep(0, length(uniq))
  for(i in  1:nrow(dat))
  {
    
    for(j in 1:length(uniq))
    {
      if(strat[i]==uniq[j]) ranef[i] = st[j]
    }
    
  }
  ranef
}


summary(mod3_sc)
#------------------------------------------------------
# Predict scale factor at grid cells
#---------------------------------------------------------------
#pred_covs_hl2$set_typ <-factor(round(pred_covs_hl2$DRC_GHSL_SMOD2))
sim_scale <- function(model, dat, Aprediction, run)
{
  fixedeff  <- scale_hat <- matrix(0, nrow=nrow(dat), ncol = run)
  #inla.seed = as.integer(runif(1)*.Machine$integer.max)
  inla.seed =  2081662007
  set.seed(inla.seed)
  print(inla.seed)
  
  # Simulated from the posterior marginal density
  m1.samp <- inla.posterior.sample(run,
                                   model, seed = inla.seed ,
                                   selection=list(
                                     x61=1,
                                     x77=1),
                                   num.threads="1:1")
  
  sfield_nodes_mean <- model$summary.random$s['mean']
  field_mean <- (Aprediction%*% as.data.frame(sfield_nodes_mean)[, 1])
  
  # Predict response with each set of simulated samples
  for(i in 1:run)
  {
    #
    fixedeff[,i] <-
      model$summary.fixed['(Intercept)', 'mean'] +
      m1.samp[[i]]$latent[1,] * dat[,'x61'] +
      m1.samp[[i]]$latent[2,] * dat[,'x77'] +
      
    
      rnorm(nrow(dat), 0, 1/model$summary.hyperpar$mean[4]) + # IID
      field_mean[,1] # Spatial random effect
    
    
    scale_hat[,i]<- exp(fixedeff[,i])
  }
  
  dat$mean_scale_hat <- apply(scale_hat, 1, mean, na.rm=T) #
  dat$sd_scale_hat <- apply(scale_hat, 1, sd, na.rm=T) #
  dat$lower_scale_hat <- apply(scale_hat, 1, quantile, probs=c(0.025), na.rm=T) #
  dat$upper_scale_hat <- apply(scale_hat, 1, quantile, probs=c(0.975), na.rm=T) #
  dat$cv_scale_hat <- dat$sd_scale_hat/dat$mean_scale_hat#
  
  # Store key outputs
  output <- list(scale_hat = scale_hat,
                 est_data = dat)
  
}
run=30


# Global scaling
sim.scale0 <- sim_scale(mod3_sc, pred_covs_hl2, Apred_sc0, run) 

pred_covs_hl2$scale0 <- sim.scale0$est_data$mean_scale_hat
boxplot(pred_covs_hl2$scale0)

# Local scaling
sim.scale1 <- sim_scale(mod3_sc_hl, pred_covs_hl2, Apred_sc1, run) 
pred_covs_hl2$scale1 <- sim.scale1$est_data$mean_scale_hat
boxplot(pred_covs_hl2$scale1)


# Predict adjusted and unadjusted population density and population counts at grid cells
simDens <- function(model, dat, Aprediction, scale_mat0, scale_mat1, run)
{
  fixedeff  <- dens_hat <- pop_hat <- pop_hat_sc0<- pop_hat_sc1 <- matrix(0, nrow=nrow(dat), ncol = run)
  #inla.seed = as.integer(runif(1)*.Machine$integer.max)
  inla.seed =  1422136236
  
  set.seed(inla.seed)
  print(inla.seed)
  m1.samp <- inla.posterior.sample(run,
                                   model, seed = inla.seed ,
                                   selection=list(
                                     x12=1,#
                                     x17=1,
                                     x55=1,
                                     x68=1
                                   ),
                                   num.threads="1:1")
  
  sfield_nodes_mean <- model$summary.random$s['mean']
  field_mean <- (Aprediction%*% as.data.frame(sfield_nodes_mean)[, 1])
  for(i in 1:run)
  {
    fixedeff[,i] <-
      model$summary.fixed['Intercept', 'mean'] +
      m1.samp[[i]]$latent[1,] * dat[,'x12'] +
      m1.samp[[i]]$latent[2,] * dat[,'x17'] +
      m1.samp[[i]]$latent[3,] * dat[,'x55'] +
      m1.samp[[i]]$latent[4,] * dat[,'x68'] +
      
      
      rnorm(nrow(dat), 0, 1/model$summary.hyperpar$mean[2]) +
      dat$set_ranef +
      field_mean[,1]
    
    dens_hat[,i]<- exp(fixedeff[,i])
    pop_hat[,i] <- dens_hat[,i]*dat$bldp
    
  }
  
  #--Unscaled
  dat$mean_dens_hat <- apply(dens_hat, 1, mean, na.rm=T) #
  dat$mean_pop_hat <- apply(pop_hat, 1, mean, na.rm=T) #
  dat$sd_pop_hat <- apply(pop_hat, 1, sd, na.rm=T) #
  dat$lower_pop_hat <- apply(pop_hat, 1, quantile, probs=c(0.025), na.rm=T) #
  dat$upper_pop_hat <- apply(pop_hat, 1, quantile, probs=c(0.975), na.rm=T) #
  dat$cv_pop_hat <- dat$sd_pop_hat/dat$mean_pop_hat#
  
  
  #--Scaling (global)
  pop_hat_sc0 <- scale_mat0*pop_hat
  dat$mean_pop_hat_sc0 <- apply(pop_hat_sc0, 1, mean, na.rm=T) #
  dat$sd_pop_hat_sc0 <- apply(pop_hat_sc0, 1,sd, na.rm=T) #
  dat$lower_pop_hat_sc0 <- apply(pop_hat_sc0, 1, quantile, probs=c(0.025), na.rm=T) #
  dat$upper_pop_hat_sc0 <- apply(pop_hat_sc0, 1, quantile, probs=c(0.975), na.rm=T) #
  dat$cv_pop_hat_sc0 <- dat$sd_pop_hat_sc0/dat$mean_pop_hat_sc0#
  
  #--Scaling (local)
  pop_hat_sc1 <- scale_mat1*pop_hat
  dat$mean_pop_hat_sc1 <- apply(pop_hat_sc1, 1, mean, na.rm=T) #
  dat$sd_pop_hat_sc1 <- apply(pop_hat_sc1, 1,sd, na.rm=T) #
  dat$lower_pop_hat_sc1 <- apply(pop_hat_sc1, 1, quantile, probs=c(0.025), na.rm=T) #
  dat$upper_pop_hat_sc1 <- apply(pop_hat_sc1, 1, quantile, probs=c(0.975), na.rm=T) #
  dat$cv_pop_hat_sc1 <- dat$sd_pop_hat_sc1/dat$mean_pop_hat_sc1#
  
  output <- list(pop_hat = pop_hat,
                 pop_hat_sc0 = pop_hat_sc0,
                 pop_hat_sc1 = pop_hat_sc1,
                 est_data = dat)
}

run = 30


#------------------------------------------------------------------------------
# ---- Building count-based predictions--------------------------------------
pred_covs_hl2$bldp <- round(pred_covs_hl2$COD_bld_count_20240617_unpublished)

pred_covs_hl2$set_ranef  <- str_ranef(pred_covs_hl2,
                                      pred_covs_hl2$DRC_GHSL_SMOD,
                                      mod3a$summary.random$set_typ$mean)
system.time(str(sim <-  simDens(mod3a,pred_covs_hl2, Apred, 
                                scale_mat0 = sim.scale0$scale_hat,
                                scale_mat1 =sim.scale1$scale_hat, 
                                run))) # model with no settlement type random effect
sum(sim$est_data$mean_pop_hat, na.rm=T) # Base(unadjusted) = 4,616,299
sum(sim$est_data$mean_pop_hat_sc0, na.rm=T) # Proposed (global) = 4,492,673
sum(sim$est_data$mean_pop_hat_sc1, na.rm=T) # Proposed (local) = 4,253,378


# Add to the grid data
# Base (unadjusted)
pred_covs_hl2$pop_hat_count <- sim$est_data$mean_pop_hat
pred_covs_hl2$pop_hat_countL <- sim$est_data$lower_pop_hat
pred_covs_hl2$pop_hat_countU <- sim$est_data$upper_pop_hat
pred_covs_hl2$pop_hat_countCV <- sim$est_data$cv_pop_hat
sum(pred_covs_hl2$pop_hat_count, na.rm=T)


# Global scaling 
pred_covs_hl2$pop_hat_count_sc0 <- sim$est_data$mean_pop_hat_sc0
pred_covs_hl2$pop_hat_count_sc0L <- sim$est_data$lower_pop_hat_sc0
pred_covs_hl2$pop_hat_count_sc0U <- sim$est_data$upper_pop_hat_sc0
pred_covs_hl2$pop_hat_count_sc0CV <- sim$est_data$cv_pop_hat_sc0
sum(pred_covs_hl2$pop_hat_count_sc0, na.rm=T)

# Local scaling
pred_covs_hl2$pop_hat_count_sc1 <- sim$est_data$mean_pop_hat_sc1
pred_covs_hl2$pop_hat_count_sc1L <- sim$est_data$lower_pop_hat_sc1
pred_covs_hl2$pop_hat_count_sc1U <- sim$est_data$upper_pop_hat_sc1
pred_covs_hl2$pop_hat_count_sc1CV <- sim$est_data$cv_pop_hat_sc1
sum(pred_covs_hl2$pop_hat_count_sc1, na.rm=T)

# spatial joins; joining the grid cell predictions to the clusters
grid2_sf <- st_as_sf(pred_covs_hl2, coords=c("lon","lat"), crs=crs(shp_sc_hl))
agg2 <- st_join(shp_sc_hl, grid2_sf) %>% group_by(clstr_d) %>%
  summarise(tot_cnt = sum(pop_hat_count,na.rm=T), # base
            tot_cnt0 = sum(pop_hat_count_sc0,na.rm=T), # global
            tot_cnt1 = sum(pop_hat_count_sc1,na.rm=T)) # local

# Convert to data frame and drop the 'geometry' variable
dat_clstr <- data.frame(shp_sc_hl)%>% dplyr::select(-geometry)
dat_grd <- as.data.frame(agg2)%>% dplyr::select(-geometry)

# Merge the summarised grid data to the observation clusters
merged_all <- merge(dat_clstr, dat_grd, by="clstr_d")


# Fit metrics function
model_metrics <- function(obs, pred, upper, lower)
{
  residual = pred - obs
  MAE = mean(abs(residual), na.rm=T)
  MSE = mean(residual^2, na.rm=T)
  MAPE = 100*sum(abs((obs-pred)/obs))/length(obs)
  RMSE = sqrt(MSE)
  BIAS = mean(residual, na.rm=T)
  corr = cor(obs[!is.na(obs)],pred[!is.na(obs)])
  output <- list(MAE  = MAE,
                 RMSE = RMSE,
                 CC = corr)
  return(output)
}

# Calculate the model fit metrics
dat_hll <- merged_all %>% filter(Microcensu!=0)
met_cnt <- unlist(model_metrics(dat_hll$Microcensu,
                                dat_hll$tot_cnt)) # Base

met_cnt0 <- unlist(model_metrics(dat_hll$Microcensu,
                                 dat_hll$tot_cnt0)) # Global scaling

met_cnt1 <- unlist(model_metrics(dat_hll$Microcensu,
                                 dat_hll$tot_cnt1)) # local scaling

(mets <- rbind(met_cnt, met_cnt0, met_cnt1))

#            MAE     RMSE      corr
#met_cnt  236.3153 373.1521 0.7065717 # unscaled
#met_cnt0 182.2466 313.2704 0.6757389  # Global scaling
#met_cnt1 158.2031 261.3533 0.7119096 # Local scaling

# Calculate reduction in RMAE
# rrmae = (1-158.2/236.3)*100 # 33.1% rrmae # local  scaling
# rrmae = (1-182.3/236.3)*100 # 22.9% rrmae # global scaling



#write.csv(mets, paste0(results_path,"/EA_health_area_metrics_updated.csv"))


# Calculate Zonal Statistics
dim(data.count <- data.frame(cbind(sim$est_data[,c("lon","lat","health_area",
                                                   "health_area_name")], sim$pop_hat)))# unscaled
dim(data.count_sc0 <- data.frame(cbind(sim$est_data[,c("lon","lat","health_area",
                                                       "health_area_name")], sim$pop_hat_sc0)))# scaled - global
dim(data.count_sc1 <- data.frame(cbind(sim$est_data[,c("lon","lat","health_area",
                                                       "health_area_name")], sim$pop_hat_sc1)))# scaled - local


#-----Calculate Provincial Totals and Uncertainties
prov_est <- function(dat, run)
{
  p_hat <- dat[,5:(run+4)]
  tots <- apply(p_hat,2, sum, na.rm=T) #Col sums
  
  tot_sd  <- sd(tots, na.rm=T)
  
  tot_mean  <- mean(tots, na.rm=T)
  
  tot_lower <- quantile(tots, probs=c(0.025))
  tot_median <- quantile(tots, probs=c(0.5))
  tot_upper <- quantile(tots, probs=c(0.975))
  
  return(estimates <- data.frame(estimates=unlist(list(total=tot_mean,
                                                       lower=tot_lower, median=tot_median, upper=tot_upper))))
}

#(prov_count <- t(prov_est(data.count, run))) # Base (unscaled)
#(prov_count_sc0 <- t(prov_est(data.count_sc0, run))) # global scaling
#(prov_count_sc1 <- t(prov_est(data.count_sc1, run))) # local scaling


##----Helath area population estimates

harea_est <- function(datr, run)
{
  names <-  unique(datr$health_area_name)
  uniR <-as.numeric(as.factor(unique(datr$health_area_name)))
  outR <- matrix(0, nrow=length(uniR), ncol=5)
  for(j in 1:length(uniR))
  {
    reg <- datr[datr$health_area_name==names[j],]
    rtots <- apply(reg[,5:(4+run)], 2, sum, na.rm=T)
    
    rtot_mean  <- mean(rtots, na.rm=T)
    rtot_sd <- sd(rtots, na.rm=T)
    
    rtot_lower <- quantile(rtots, probs=c(0.025))
    rtot_median <- quantile(rtots, probs=c(0.5))
    rtot_upper <- quantile(rtots, probs=c(0.975))
    rtot_uncert <- (rtot_upper - rtot_lower)/rtot_mean
    
    restimates <- round(c(rtot_mean, rtot_lower, rtot_median,rtot_upper, rtot_uncert),2)
    outR[j,] <- restimates
  }
  outR <- data.frame(outR)
  return(reg_est <- data.frame(id = uniR,
                               names = names,
                               total = outR[,1],
                               lower = outR[,2],
                               median = outR[,3],
                               upper = outR[,4],
                               uncertainty = outR[,5]))
}


#dim(harea_count <- harea_est(data.count[!is.na(data.count$health_area),], run))
#harea_count  <- harea_count[order(harea_count$names),]

#dim(harea_count_sc0 <- harea_est(data.count_sc0[!is.na(data.count_sc0$health_area),], run))
#harea_count_sc0  <- harea_count_sc0[order(harea_count_sc0$names),]

#dim(harea_count_sc1 <- harea_est(data.count_sc1[!is.na(data.count_sc1$health_area),], run))
#harea_count_sc1  <- harea_count_sc1[order(harea_count_sc1$names),]


# Add variables
#harea_count_sc0$total_cnt <- harea_count$total
#harea_count_sc0total_cnt0 <- harea_count_sc0$total
#harea_count_sc0$total_cnt1 <- harea_count_sc1$total


#write.csv(harea_count_sc0, paste0(results_path, "/health_area_totals_combined.csv"), row.names=F)



#-----------------------------------------------------------------------
# Writing raster files
ref_coords <- cbind(pred_covs_hl2$lon, pred_covs_hl2$lat)
x <- as.matrix(ref_coords)


# Base (unscaled)
# mean
z1a <- as.matrix(pred_covs_hl2$pop_hat_count)
h1a= rasterFromXYZ(cbind(x, z1a))
writeRaster(h1a, filename=paste0(results_path, "/base/Haut-Lomami_mean_count.tif"),
            overwrite=TRUE, options = c('COMPRESS' = 'LZW'))

z1aCV <- as.matrix(pred_covs_hl2$pop_hat_countCV)
h1aCV= rasterFromXYZ(cbind(x, z1aCV))
writeRaster(h1aCV, filename=paste0(results_path, "/base/Haut-Lomami_CV_count.tif"),
            overwrite=TRUE, options = c('COMPRESS' = 'LZW'))

# Global scaling
# mean
z1a_sc0 <- as.matrix(pred_covs_hl2$pop_hat_count_sc0)
h1a_sc0= rasterFromXYZ(cbind(x, z1a_sc0))
writeRaster(h1a_sc0, filename=paste0(results_path, "/scaled_all/Haut-Lomami_mean_count_sc0.tif"),
            overwrite=TRUE, options = c('COMPRESS' = 'LZW'))
# CV
z1a_sc0CV <- as.matrix(pred_covs_hl2$pop_hat_count_sc0CV)
h1a_sc0CV= rasterFromXYZ(cbind(x, z1a_sc0CV))
writeRaster(h1a_sc0CV, filename=paste0(results_path, "/scaled_all/Haut-Lomami_CV_count_sc0.tif"),
            overwrite=TRUE, options = c('COMPRESS' = 'LZW'))



## Local scaling
# mean
z1a_sc1 <- as.matrix(pred_covs_hl2$pop_hat_count_sc1)
h1a_sc1= rasterFromXYZ(cbind(x, z1a_sc1))
writeRaster(h1a_sc1, filename=paste0(results_path, "/scaled_hl/Haut-Lomami_mean_count_sc1.tif"),
            overwrite=TRUE, options = c('COMPRESS' = 'LZW'))
#   CV
z1a_sc1CV <- as.matrix(pred_covs_hl2$pop_hat_count_sc1CV)
h1a_sc1CV= rasterFromXYZ(cbind(x, z1a_sc1CV))
writeRaster(h1a_sc1CV, filename=paste0(results_path, "/scaled_hl/Haut-Lomami_CV_count_sc1.tif"),
            overwrite=TRUE, options = c('COMPRESS' = 'LZW'))


##------------------------------------------------------------------
#    Model Cross Validation (K-fold)
#-------------------------------------------------------------------
# Function: Model fit metrics
mod_metrics2 <- function(obs, pred, upper, lower)
{
  residual = pred - obs
  MAE = mean(abs(residual), na.rm=T)
  MSE = mean(residual^2, na.rm=T)
  MAPE = 100*sum(abs((obs-pred)/obs))/length(obs)
  RMSE = sqrt(MSE)
  corr = cor(obs[!is.na(obs)],pred[!is.na(obs)])
  
  output <- list(MAE  = MAE,
                 RMSE = RMSE,
                 CC= corr)
  return(output)
}

#----Extract settement type effects
set_t <- function(dat, st)
{
  uniq <- unique(dat$set_typ)
  uniq[1]
  for(i in  1:nrow(dat))
  {
    
    for(j in 1:3)
    {
      if(dat$set_typ[i]==uniq[j]) dat$set_typ2[i] = st[j]
    }
    
  }
  dat$set_typ2
}



# Cross-validation for the proposed
cross_validate <- function(dat, n.folds, mod1, formula1,
                           A, mesh, spde,iset, seed)
{
  #--------------------------------------------
  # dat: the input survey data containing the all the variables
  # n.folds: number of test (k) folds to use
  # mod: the best model of the full or reference data
  # A: the projection  matrix used in training the full data model
  # seed: a random sample seed to make results reproducible
  #--------------------------------------------
  #seed = 1325
  set.seed(seed)
  
  N <- nrow(dat)
  
  # train-test sets
  table(ind_train <- factor(sample(x = rep(1:n.folds,c(rep(45,4), 48)),  # Sample IDs for training data
                                   size = N))) 
  
  
  table(as.numeric(ind_train))
  dat$k_fold <- as.numeric(ind_train)
  coords = cbind(dat$lon, dat$lat)
  
  
  k_uniq <-sort(unique(dat$k_fold))
  
  # extract settlement type random effect
  dat$set_ranef1  <- str_ranef(dat,dat$set_typ,
                               mod1$summary.random$set_typ$mean)
  #---------------------------------------------------------------
  #                   in-sample
  #---------------------------------------------------------------
  
  met_list_in <- list()
  pred_list_in <- list()
  for(i in 1:length(k_uniq))
  {
    print(paste0("in-sample cross-validation using fold ", i, sep=""))
    test_ind <- which(dat$k_fold==k_uniq[i])
    
    dim(test <- dat[test_ind, ]) #---test set for fold i
    
    
    train_coords <- coords
    test_coords <- coords[test_ind,]
    
    
    
    # spatial random effects based on the full data best model
    sfield_nodes_mean1 <- mod1$summary.random$s['mean']
    field_mean1 <- (A%*% as.data.frame(sfield_nodes_mean1)[, 1])
    
    
    
    # Mean
    fixed1 <-
      mod1$summary.fixed['Intercept', 'mean'] +
      mod1$summary.fixed['x12', 'mean'] * test[,'x12'] +
      mod1$summary.fixed['x17', 'mean'] * test[,'x17'] +
      mod1$summary.fixed['x55', 'mean'] * test[,'x55'] +
      mod1$summary.fixed['x68', 'mean'] * test[,'x68'] +
      
      test$set_ranef1 +
      mod1$summary.random$eps['mean'][test_ind,1] + #--uncorrelated spatial random effects
      field_mean1[test_ind,1]
    dens_ht1 <- exp(fixed1)
    sum(test$pop_ht1 <- dens_ht1*test$bld, na.rm=T)
    sum(test$pop_ht <- test$pop_ht1*test$scale_prd, na.rm=T)
    
    
    
    par(mfrow =c(1,1))
    plot(test$obs, test$pop_ht, xlab = "Observed",
         ylab = "Predicted", col=c('blue','orange'),
         pch=c(16,16), cex.axis=1.5)
    abline(0,1)
    legend("topleft", c("Observed", "Predicted"), col=c("blue", "orange"), pch=c(16,16),
           bty="n", cex=1.5)
    
    
    
    # calculate fit metrics
    met_in <- mod_metrics2(test$obs,
                           test$pop_ht)
    
    met_list_in[[i]]<- unlist(met_in)
    pred_list_in[[i]] <- data.frame(obs = round(test$obs), pred = round(test$pop_ht),
                                    fold = rep(i, length(test$obs)),
                                    data = rep("insample", length(test$obs)))
  }
  met_list_in_dat <- do.call(rbind,met_list_in)
  metrics_in <- apply(met_list_in_dat, 2, mean,na.rm=T)
  pred_list_in_dat <- do.call(rbind,pred_list_in)
  
  #-----------------------------------------------------------
  #               out - of -sample
  #-----------------------------------------------------------
  
  met_list_out <- list()
  pred_list_out <- list()
  for(i in 1:length(k_uniq))
  {
# i =4
    # dat = dat_sc_hl
    
    print(paste0("out-of-sample cross-validation using fold ", i, sep=""))
    train_ind <- which(dat$k_fold!=k_uniq[i])
    test_ind <- which(dat$k_fold==k_uniq[i])
    dim(train <- dat[train_ind, ])#---train set for fold i
    dim(test <- dat[test_ind, ]) #---test set for fold i
    
    
    train_coords <- coords[train_ind,]
    test_coords <- coords[test_ind,]
    
    # mesh = mesh_sc_hl
    ###---Create projection matrices for training and testing datasets
    Ae<-inla.spde.make.A(mesh=mesh,loc=as.matrix(train_coords));dim(Ae) #training
    
    
    #####------------------------
    # population density
    # iset = iset_sc_hl
    covars_train <- train[,c("x12", "x17", "x55","x68", "set_typ", "eps")]; dim(covars_train)
    stk_train_dens <- inla.stack(data=list(y=train$dens1), #the response
                                 
                                 A=list(Ae,1),  #the A matrix; the 1 is included to make the list(covariates)
                                 
                                 effects=list(c(list(Intercept=1), #the Intercept
                                                iset),  #the spatial index
                                              #the covariates
                                              list(covars_train)
                                 ),
                                 tag='train')
    
    
    # formula1 = f3a
    model1 <-inla(formula1, #the formula
                  data=inla.stack.data(stk_train_dens),  #the data stack
                  family= 'gamma',   #which family the data comes from
                  control.predictor=list(A=inla.stack.A(stk_train_dens),compute=TRUE),  #compute gives you the marginals of the linear predictor
                  control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                  verbose = FALSE) #can include verbose=TRUE to see the log of the model runs
    summary(model1)
    
    
    # Extract Spatial random effects from the full data best model
    sfield_nodes_mean1 <- mod1$summary.random$s['mean']
    field_mean1 <- (A%*% as.data.frame(sfield_nodes_mean1)[, 1])
    
    ##--------
    fixed1 <-
      model1$summary.fixed['Intercept', 'mean'] +
      model1$summary.fixed['x12', 'mean'] * test[,'x12'] +
      model1$summary.fixed['x17', 'mean'] * test[,'x17'] +
      model1$summary.fixed['x55', 'mean'] * test[,'x55'] +
      model1$summary.fixed['x68', 'mean'] * test[,'x68'] +
      
      test$set_ranef1 +
      mod1$summary.random$eps['mean'][test_ind,1] + #--uncorrelated spatial random effects
      field_mean1[test_ind,1]
    dens_ht1 <- exp(fixed1)
    sum(test$pop_ht1 <- dens_ht1*test$bld, na.rm=T)
    # 
    sum(test$pop_ht <- test$pop_ht1*test$scale_prd, na.rm=T)
    
    
    # Visualise scatter plot
    par(mfrow =c(1,1))
    plot(test$obs, test$pop_ht, xlab = "Observed",
         ylab = "Predicted", col=c('blue','orange'),
         pch=c(16,16), cex.axis=1.5)
    abline(0,1)
    legend("topleft", c("Observed", "Predicted"), col=c("blue", "orange"), pch=c(16,16),
           bty="n", cex=1.5)
    
    
    
    # calculate fit metrics
    met_out<- mod_metrics2(test$obs,
                           test$pop_ht)
    
    met_list_out[[i]]<- unlist(met_out)
    pred_list_out[[i]] <- data.frame(obs = round(test$obs), pred = round(test$pop_ht),
                                     fold = rep(i, length(test$obs)),
                                     data = rep("out-of-sample", length(test$obs)))
  }
  met_list_out_dat <- do.call(rbind,met_list_out)
  metrics_out <- apply(met_list_out_dat, 2, mean, na.rm=T) # fit metrics
  pred_list_out_dat <- do.call(rbind,pred_list_out)# predictions
  
  cv_mets <- rbind(metrics_in, metrics_out)
  output <- list( met_list_in_dat = met_list_in_dat,
                  met_list_out_dat = met_list_out_dat,
                  pred_dat = rbind(pred_list_in_dat, pred_list_out_dat),
                  cv_metrics = rbind(metrics_in, metrics_out))
}

dat <- dat_sc_hl
dat$bld <- dat$sum.COD_bld_count_20240617_unpublished
dat$obs <- dat$Microcensu
dat$scale_prd <- dat$pred_sc
(cross_val <- cross_validate(dat,
                             n.folds = 5,
                             mod1 = mod3a,
                             formula1 = f3a,
                             A = A_sc_hl,
                             mesh = mesh_sc_hl,
                             spde = spde_sc_hl,
                             iset = iset_sc_hl,
                             seed = 195))

cross_val$met_list_in_dat  # in-sample metrics per fold
cross_val$met_list_out_dat  # out-of-sample metrics per fold
cross_val$cv_metrics    # combined averaged metrics

met_in1 <- data.frame(cross_val$met_list_in_dat)  # in-sample metrics per fold
#       MAE     RMSE        CC
#[1,] 101.58330 152.4707 0.7720217
#[2,] 138.29738 261.6922 0.5507079
#[3,] 117.80399 193.0491 0.7082482
#[4,] 103.62658 189.5738 0.7624510
#[5,]  67.54239 109.3930 0.8466975
met_in1$cv <- rep("In-sample", nrow(met_in1))
met_in1$Method <- rep("Proposed", nrow(met_in1))

met_out1 <- data.frame(cross_val$met_list_out_dat)  # out-of-sample metrics per fold
#       MAE     RMSE        CC
#[1,]  96.16673 146.8202 0.7944133
#[2,] 131.89954 240.7415 0.5723140
#[3,] 114.41926 187.5725 0.7226567
#[4,] 104.92257 196.5021 0.7642049
#[5,]  73.67346 117.4290 0.8315854
met_out1$cv <- rep("Out-of-sample", nrow(met_out1))
met_out1$Method <- rep("Proposed", nrow(met_out1))
cross_val$cv_metrics    # combined averaged metrics
#               MAE     RMSE        CC
#metrics_in  105.7707 181.2358 0.7280253
#metrics_out 104.2163 177.8131 0.7370349

cross_val$pred_dat  # combined prediction data

meet1 <- rbind(met_in1, met_out1)


# Cross-validation for the base method
#---------------------------------------------
cross_validate2 <- function(dat, n.folds, mod1, formula1,
                            A, mesh, spde,iset, seed)
{
  #--------------------------------------------
  # dat: the input survey data containing the all the variables
  # n.folds: number of test (k) folds to use
  # mod: the best model of the full or reference data
  # A: the projection  matrix used in training the full data model
  # seed: a random sample seed to make results reproducible
  #--------------------------------------------
  #seed = 1325
  set.seed(seed)
  N <- nrow(dat)
  
  
  
  ######
  table(ind_train <- factor(sample(x = rep(1:n.folds,c(rep(45,4), 48)),  # Sample IDs for training data
                                   size = N))) 
  
  
  table(as.numeric(ind_train))
  dat$k_fold <- as.numeric(ind_train)
  coords = cbind(dat$lon, dat$lat)
  
  
  k_uniq <-sort(unique(dat$k_fold))
  
  # mod1 = mod3c
  # A = A_hl
  #dat$bld <- dat$sum.COD_bld_count_20240617_unpublished
  #dat$scale_prd <- dat$pred_sc1
  # dat$obs <- dat$pop
  dat$set_ranef1  <- str_ranef(dat,dat$set_typ,
                               mod1$summary.random$set_typ$mean)
  #---------------------------------------------------------------
  #                   in-sample
  #---------------------------------------------------------------
  
  met_list_in <- list()
  pred_list_in <- list()
  for(i in 1:length(k_uniq))
  {
    #i =3
    
    print(paste0("in-sample cross-validation using fold ", i, sep=""))
    test_ind <- which(dat$k_fold==k_uniq[i])
    
    dim(test <- dat[test_ind, ]) #---test set for fold i
    
    
    train_coords <- coords
    test_coords <- coords[test_ind,]
    
    
    
    # spatial random effects based on the full data best model
    sfield_nodes_mean1 <- mod1$summary.random$s['mean']
    field_mean1 <- (A%*% as.data.frame(sfield_nodes_mean1)[, 1])
    
    
    
    ##--------
    fixed1 <-
      mod1$summary.fixed['Intercept', 'mean'] +
      mod1$summary.fixed['x12', 'mean'] * test[,'x12'] +
      mod1$summary.fixed['x17', 'mean'] * test[,'x17'] +
      mod1$summary.fixed['x55', 'mean'] * test[,'x55'] +
      mod1$summary.fixed['x68', 'mean'] * test[,'x68'] +
      
      test$set_ranef1 +
      mod1$summary.random$eps['mean'][test_ind,1] + #--uncorrelated spatial random effects
      field_mean1[test_ind,1]
    dens_ht1 <- exp(fixed1)
    sum(test$pop_ht <- dens_ht1*test$bld, na.rm=T)
    
    
    
    
    par(mfrow =c(1,1))
    plot(test$obs, test$pop_ht, xlab = "Observed",
         ylab = "Predicted", col=c('blue','orange'),
         pch=c(16,16), cex.axis=1.5)
    abline(0,1)
    legend("topleft", c("Observed", "Predicted"), col=c("blue", "orange"), pch=c(16,16),
           bty="n", cex=1.5)
    
    
    
    # calculate fit metrics
    met_in <- mod_metrics2(test$obs,
                           test$pop_ht)
    
    met_list_in[[i]]<- unlist(met_in)
    pred_list_in[[i]] <- data.frame(obs = round(test$obs), pred = round(test$pop_ht),
                                    fold = rep(i, length(test$obs)),
                                    data = rep("insample", length(test$obs)))
  }
  met_list_in_dat <- do.call(rbind,met_list_in)
  metrics_in <- apply(met_list_in_dat, 2, mean,na.rm=T)
  pred_list_in_dat <- do.call(rbind,pred_list_in)
  
  #-----------------------------------------------------------
  #               out - of -sample
  #-----------------------------------------------------------
  
  met_list_out <- list()
  pred_list_out <- list()
  for(i in 1:length(k_uniq))
  {
    #i=5
    print(paste0("out-of-sample cross-validation using fold ", i, sep=""))
    train_ind <- which(dat$k_fold!=k_uniq[i])
    test_ind <- which(dat$k_fold==k_uniq[i])
    dim(train <- dat[train_ind, ])#---train set for fold i
    dim(test <- dat[test_ind, ]) #---test set for fold i
    
    
    train_coords <- coords[train_ind,]
    test_coords <- coords[test_ind,]
    
    
    ###---Create projection matrices for training and testing datasets
    Ae<-inla.spde.make.A(mesh=mesh,loc=as.matrix(train_coords));dim(Ae) #training
    
    
    #####------------------------
    # population density
    covars_train <- train[,c("x12", "x17", "x55","x68", "set_typ", "eps")]; dim(covars_train)
    stk_train_dens <- inla.stack(data=list(y=train$dens1), #the response
                                 
                                 A=list(Ae,1),  #the A matrix; the 1 is included to make the list(covariates)
                                 
                                 effects=list(c(list(Intercept=1), #the Intercept
                                                iset),  #the spatial index
                                              #the covariates
                                              list(covars_train)
                                 ),
                                 tag='train')
    
    
    
    model1 <-inla(formula1, #the formula
                  data=inla.stack.data(stk_train_dens),  #the data stack
                  family= 'gamma',   #which family the data comes from
                  control.predictor=list(A=inla.stack.A(stk_train_dens),compute=TRUE),  #compute gives you the marginals of the linear predictor
                  control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                  verbose = FALSE) #can include verbose=TRUE to see the log of the model runs
    summary(model1)
    
    
    # Extract Spatial random effects from the full data best model
    sfield_nodes_mean1 <- mod1$summary.random$s['mean']
    field_mean1 <- (A%*% as.data.frame(sfield_nodes_mean1)[, 1])
    
    ##--------
    fixed1 <-
      model1$summary.fixed['Intercept', 'mean'] +
      model1$summary.fixed['x12', 'mean'] * test[,'x12'] +
      model1$summary.fixed['x17', 'mean'] * test[,'x17'] +
      model1$summary.fixed['x55', 'mean'] * test[,'x55'] +
      model1$summary.fixed['x68', 'mean'] * test[,'x68'] +
      
      test$set_ranef1 +
      mod1$summary.random$eps['mean'][test_ind,1] + #--uncorrelated spatial random effects
      field_mean1[test_ind,1]
    dens_ht1 <- exp(fixed1)
    sum(test$pop_ht <- dens_ht1*test$bld, na.rm=T)
    
    
    
    
    par(mfrow =c(1,1))
    plot(test$obs, test$pop_ht, xlab = "Observed",
         ylab = "Predicted", col=c('blue','orange'),
         pch=c(16,16), cex.axis=1.5)
    abline(0,1)
    legend("topleft", c("Observed", "Predicted"), col=c("blue", "orange"), pch=c(16,16),
           bty="n", cex=1.5)
    
    
    
    # calculate fit metrics
    met_out<- mod_metrics2(test$obs,
                           test$pop_ht)
    
    met_list_out[[i]]<- unlist(met_out)
    pred_list_out[[i]] <- data.frame(obs = round(test$obs), pred = round(test$pop_ht),
                                     fold = rep(i, length(test$obs)),
                                     data = rep("out-of-sample", length(test$obs)))
  }
  met_list_out_dat <- do.call(rbind,met_list_out)
  metrics_out <- apply(met_list_out_dat, 2, mean, na.rm=T) # fit metrics
  pred_list_out_dat <- do.call(rbind,pred_list_out)# predictions
  
  cv_mets <- rbind(metrics_in, metrics_out)
  output <- list( met_list_in_dat = met_list_in_dat,
                  met_list_out_dat = met_list_out_dat,
                  pred_dat = rbind(pred_list_in_dat, pred_list_out_dat),
                  cv_metrics = rbind(metrics_in, metrics_out))
}

dat <- dat_sc_hl
dat$bld <- dat$sum.COD_bld_count_20240617_unpublished
dat$obs <- dat$Microcensu
(cross_val2 <- cross_validate2(dat,
                               n.folds = 5,
                               mod1 = mod3a,
                               formula1 = f3a,
                               A = A_sc_hl,
                               mesh = mesh_sc_hl,
                               spde = spde_sc_hl,
                               iset = iset_sc_hl,
                               seed = 165))

cross_val2$met_list_in_dat  # in-sample metrics per fold
cross_val2$met_list_out_dat  # out-of-sample metrics per fold
cross_val2$cv_metrics    # combined averaged metrics


met_in2 <- data.frame(cross_val2$met_list_in_dat)  # in-sample metrics per fold
#       MAE     RMSE        CC
#[1,] 205.8563 316.3214 0.7288364
#[2,] 218.4063 431.1797 0.5720546
#[3,] 230.4920 363.7294 0.7316758
#[4,] 210.2595 339.8089 0.7436682
#[5,] 118.7834 189.7433 0.8343824
met_in2$cv <- rep("In-sample", nrow(met_in2))
met_in2$Method <- rep("Base", nrow(met_in2))


met_out2 <- data.frame(cross_val2$met_list_out_dat)  # out-of-sample metrics per fold
#        MAE     RMSE        CC
#[1,] 188.5227 284.6584 0.7668786
#[2,] 210.6226 392.4255 0.5965018
#[3,] 228.7406 350.2058 0.7524699
#[4,] 209.9067 347.5268 0.7359415
#[5,] 126.6821 205.6746 0.8179419
met_out2$cv <- rep("Out-of-sample", nrow(met_out2))
met_out2$Method <- rep("Base", nrow(met_out2))


meet2 <- rbind(met_in2, met_out2)
metrics_all <- rbind(meet1, meet2)

cross_val2$cv_metrics    # combined averaged metrics
#               MAE     RMSE        CC
# metrics_in  196.7595 328.1565 0.7221235
# metrics_out 192.8949 316.0983 0.7339467

# Base 
Base_mets <- data.frame(cross_val2$cv_metrics)
Base_mets$Method <- rep("Base", nrow(Base_mets))

# Proposed 
Prop_mets <- data.frame(cross_val$cv_metrics)
Prop_mets$Method <- rep("Proposed", nrow(Prop_mets))
# Combined average metrics
rbind(Base_mets, Prop_mets)


#-------------------------------------------------------------
# cross-validation results plots
#------------------------------------------------------------
# Boxplots for fit metrics

# MAE
plot_mae <- metrics_all %>%
  ggplot( aes(x=cv, y=MAE, fill=Method)) +
  #geom_violin(width=1.4) +
  geom_boxplot(width=0.8, 
               alpha=0.8,
               
               notchwidth = 0.2) +
  theme_bw()+
  theme(strip.text = element_text(size = 20),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=14))


rmae <-  ggpar(plot_mae, xlab="Cross-Validation", ylab="MAE",
               legend = "top",
               legend.title = "Method:",size=20,,
               font.legend=c(20),
               font.label = list(size = 15, face = "bold", color ="red"),
               palette = "lancet",
               font.y = c(20),
               font.x = c(20),
               font.main=c(20),
               font.xtickslab =c(18),
               font.ytickslab =c(18),
               ytickslab.rt = 45)
rmae


# RMSE
plot_rmse <- metrics_all %>%
  ggplot( aes(x=cv, y=RMSE, fill=Method)) +

  geom_boxplot(width=0.8, 
               alpha=0.8,
              
               notchwidth = 0.2) +
  theme_bw()+
  theme(strip.text = element_text(size = 20),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=14))


rrmse <-  ggpar(plot_rmse, xlab="Cross-Validation", ylab="RMSE",
                legend = "top",
                legend.title = "Method:",size=20,,
                font.legend=c(20),
                font.label = list(size = 15, face = "bold", color ="red"),
        
                palette = "lancet",
                font.y = c(20),
                font.x = c(20),
                font.main=c(20),
                font.xtickslab =c(18),
                font.ytickslab =c(18),
                ytickslab.rt = 45)
rrmse



# CC
plot_cc <- metrics_all %>%
  ggplot(aes(x=cv, y=CC, fill=Method)) +

  geom_boxplot(width=0.8, 
               alpha=0.8,
               notchwidth = 0.2) +
  theme_bw()+
  theme(strip.text = element_text(size = 20),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=14))


rcc <-  ggpar(plot_cc, xlab="Cross-Validation", ylab="Correlation Coefficient",
              legend = "top",
              legend.title = "Method:",size=20,,
              font.legend=c(20),
              font.label = list(size = 15, face = "bold", color ="red"),
            
              palette = "lancet",
              font.y = c(20),
              font.x = c(20),
              font.main=c(20),
              font.xtickslab =c(18),
              font.ytickslab =c(18),
              ytickslab.rt = 45)
rcc


### Scatter plots of observed versus predicted
dat_hl_ha$obs
harea_count_sc0$total_cnt <- harea_count$total
harea_count_sc0total_cnt0 <- harea_count_sc0$total
harea_count_sc0$total_cnt1 <- harea_count_sc1$total


# Combine the datasets
harea_count$Method <- rep("Base", nrow(harea_count))
harea_count_sc1$Method <- rep("Proposed", nrow(harea_count_sc1))

ha_dat <- rbind(harea_count, harea_count_sc1)
ha_dat$Method <- factor(ha_dat$Method)

# Density plots with means
library(plyr)
cdat <- ddply(ha_dat, "Method", summarise, Mean=mean(total))
cdat

# Overlaid histograms with means

# Mean
hist_mean <- ggplot(ha_dat, aes(x=total, fill=Method)) +
  geom_histogram(bin=30, alpha=.5, color = "black",position="identity") +
  geom_vline(data=cdat, aes(xintercept=Mean,  colour=Method),
             linetype="dashed", size=1)+
  theme_minimal()

rhist_mean <-  ggpar(hist_mean, xlab="Predicted count", ylab="Frequency",
                     legend = "top",
                     legend.title = "Method:",size=20,,
                     font.legend=c(20),
                     font.label = list(size = 18, face = "bold", color ="red"),
              
                     palette = "lancet",
                     font.y = c(20),
                     font.x = c(20),
                     font.main=c(20),
                     font.xtickslab =c(18),
                     font.ytickslab =c(18),
                     ytickslab.rt = 45)
rhist_mean


# Lower
hist_lwr <- ggplot(ha_dat, aes(x=lower, fill=Method)) +
  geom_histogram(bin=30, alpha=.5, color = "black",position="identity") +
  theme_minimal()

rhist_lwr <-  ggpar(hist_lwr, xlab="Predicted count", ylab="Frequency",
                    legend = "top",
                    legend.title = "Method:",size=20,,
                    font.legend=c(20),
                    font.label = list(size = 15, face = "bold", color ="red"),
            
                    palette = "lancet",
                    font.y = c(20),
                    font.x = c(20),
                    font.main=c(20),
                    font.xtickslab =c(18),
                    font.ytickslab =c(18),
                    ytickslab.rt = 45)
rhist_lwr



# Upper
hist_upr <- ggplot(ha_dat, aes(x=upper, fill=Method)) +
  geom_histogram(bin=30, alpha=.5, color = "black",position="identity") +
  theme_minimal()

rhist_upr <-  ggpar(hist_upr, xlab="Predicted count", ylab="Frequency",
                    legend = "top",
                    legend.title = "Method:",size=20,,
                    font.legend=c(20),
                    font.label = list(size = 15, face = "bold", color ="red"),
                 
                    palette = "lancet",
                    font.y = c(20),
                    font.x = c(20),
                    font.main=c(20),
                    font.xtickslab =c(18),
                    font.ytickslab =c(18),
                    ytickslab.rt = 45)
rhist_upr


#---------------------------------------------------
# Posterior maps
#---------------------------------------------------
#--------------------------------------------------
####-----Make spatial maps
#----------------------------------------------------

# rename health area in the shapefile
harea_hl$names <- harea_hl$airesante
harea_hl <- harea_hl[order(harea_hl$names),]
harea_count <- harea_count[order(harea_count$names),]
harea_count_sc1 <- harea_count_sc1[order(harea_count_sc1$names),]

ha_dat <- merge(harea_hl, harea_count_sc1, by="names")

# merge  with the aggregated health area datasets
vars2add <- c("names", "total", "lower",
              "upper")
dim(ha_dat.shp <- ha_dat[,vars2add])



library(ceramic)
library(basemaps)
library(ggmap)

# Create basemap for DRC
basemap <- basemap_terra(ha_dat.shp, map_service = "esri", map_type = "natgeo_world_map")


# Plot Health area maps
tmap_options(check.and.fix = TRUE)


### Scaled
# total
min(ha_dat.shp$lower); max(ha_dat.shp$upper)
total <-
  tm_shape(ha_dat.shp) +
  tm_borders(lwd=1, col="black")+
  tm_polygons(c("total"),
              palette=plasma(100),
              title=c("Population \n count"),
              legend.show=F,
              breaks=c(750, 2500, 5000, 15000, 35000,
                       50000,75000)
  )+
  tm_layout(frame=F, legend.outside = T,
            legend.text.size = 0.5, legend.title.size = 2) #+

# lower
lower <- tm_shape(basemap) +
  tm_rgb() +
  tm_shape(ha_dat.shp) +
  tm_borders(lwd=1, col="black")+
  tm_polygons(c("lower"),
              palette=plasma(100),
              title=c("Population \n count"),
              legend.show=F,
              breaks=c(750, 2500, 5000, 15000, 35000,
                       50000,75000)
  )+
  tm_layout(frame=F, legend.outside = T,
            legend.text.size = 0.5, legend.title.size = 2) 


# upper
upper <- tm_shape(basemap) +
  tm_rgb() +
  tm_shape(ha_dat.shp) +
  tm_borders(lwd=1, col="black")+
  tm_polygons(c("upper"),
              palette=plasma(100),
              title=c("Population \n count"),
              legend.show=F,
              breaks=c(750, 2500, 5000, 15000, 35000,
                       50000,75000)
  )+
  tm_scale_bar(size = 0.9, position = c("right", "bottom"))+
  tm_compass(size=2,position = c("right", "top"))+
  tm_layout(frame=F, legend.outside = T,
            legend.text.size = 0.5, legend.title.size = 2) 



#   Mapping the grid data
# scaled
hl_grd.count_sc <- rast(paste0(results_path, "/scaled_hl/Haut-Lomami_mean_count_sc1.tif"))
plot(hl_sp, lwd=2,col="black")
plot(hl_grd.count_sc, main="Gridded population - count scaled",
     col=viridis(255), add= T)

# CV
hl_grd.count_scCV <- rast(paste0(results_path, "/scaled_hl/Haut-Lomami_CV_count_sc1.tif"))
plot(hl_sp, lwd=2,col="black")
plot(hl_grd.count_scCV, main="Gridded population - count scaled",
     col=viridis(100), add= T)


save.image(paste0(results_path, "/adjustment_factor_Haut-Lomami_example.Rdata"))
