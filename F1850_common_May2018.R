#
#
# ---------------------------------------
# Disclaimer: this script is adapted from Doug McNeall's script famous_common.R,
# located in the famous-git directory at github.com/dougmcneall/famous-git
#
# ---------------------------------------
# cgf: April 20, 2018
# this version works with the revised F1850 simulations run by Bakr
# that have fixed a bug with the x parameter.
# The previous F1850 simulations did not perturb sulfate.
#
# Modifications: added threshold-setting here, so it can be used by multiple routines consistently.
#
# May 4, 2018: re-ran n=350 training cases on Graham. Use those data here.
# ---------------------------------------
#
# Doug McNeall's helper functions
# plus my own data reading directories
library(DiceKriging)
library(RColorBrewer)
library(MASS)
library(ncdf4)
library(RNetCDF)
library(fields)
library(parallel)
library(randtoolbox)
library(dplyr)
library(sensitivity)

pth="/home/mclean80/FKB2018/"
source(paste(pth,'emtools.R',sep=""))
source(paste(pth,'imptools.R',sep=""))
source(paste(pth,'vistools.R',sep=""))

# data directory
datdir <- 'data/'
dd <- paste(pth,datdir,sep="")
set.seed(99)

# where to save the plots
plotpath <- 'plots/'
plotpath <- paste(pth,plotpath,sep="")

# new: Google Drive folder for uploads
# gpth="Projects/Aerosols/FKB2018/figs_May2018/"

# ---------------------------------------
# cgf: Read all data required for CAM UQ analysis
# ---------------------------------------
#
# First, set up thresholds for subsetting parameter space:
ss_thresh=0.85
aod_thresh=0.08 # cgf mod: new larger threshold for AOD
cs_lthresh=0.4175423  ###this is the 2.5th Percentile of plausible emulated outputs for Cess CS
cs_uthresh=0.5382148    ###this is the 97.5th Percentile of plausible emulated outputs for Cess CS

setwd(dd) # move to data directory

# first get the column headings from an existing file that has them in:
# read the (global mean Delta) training data from the simulator outputs:
dum=read.csv("mean.csv",head=T)
colnames=colnames(dum)

# Next read the mean, var and cor for this experiment, and extract only the vars we need
mns=read.table("results_table.2000.FNET.annual_mean.3members.dat",head=T)

# next read the mean SS values (ssm) produced by F1850_PierceSS_Apr2018.R:
ss=read.csv("F1850_PierceSS_output_reduced_May2018.csv",head=T)
ssm=ss$mean # extract mean over all variables

# NEW: read Cess CS (cs)
cs=read.table("results_table.cess.3members.dat",head=T)
cs=cs$CESS_fnet

# lastly, read the parameter values for the training cases, which are the predictors (design matrix):
params=read.table("parameters_F1850_350cases_9par_2017-08-04.txt",head=F,skip=2)

# set up proper col names for params:
colnames(params)=c("x1","x2","x3","x4","x5","x6","x7","x8","x9")

# set up design matrix variable X:
X=params
params.norm=normalize(params) # normalized version of params scaled to between [0,1]

# standard/default parameter values (including NAs for those undefined by default in CAM)
X.standard=c(0.0 , 0.0 , 1.00,  NA,0.88,14.00,1800.00, 0.50,3600.00)

# Bind all data together in a data frame: cols 1-9 are the INPUTS, cols 10-16 are the outputs,  17 is ssm and 18 is cs.
# Note: full_data contains the unstandardized parameter values
full_data=cbind(params,mns$AEROD_V,mns$CLDL,mns$FNET,mns$LWCF,mns$PRECT,mns$QRL,mns$SWCF,ssm,cs)
colnames(full_data)[10:16]=c("AEROD_V","CLDL","FNET","LWCF","PRECT","QRL","SWCF")

# Insert default parameter values from CAM into the full_data frame:
X.standard.append <- c(X.standard, rep(NA, 9))
# also get normalized default predictors
params.stand<-rbind(params, X.standard) # append (unnormalized) default parameters to params.stand
params.stand.norm=normalize(params.stand) # now calculate normalized parameters with defaults included
X.standard.norm.append <- c(params.stand.norm[length(params.stand.norm[,1]),], rep(NA, 9))

# full_data.norm is the same as full_data, but with normalized parameters
full_data.norm=as.data.frame(cbind(params.norm,mns$AEROD_V,mns$CLDL,mns$FNET,mns$LWCF,mns$PRECT,mns$QRL,mns$SWCF,ssm,cs))
colnames(full_data.norm)[10:16]=c("AEROD_V","CLDL","FNET","LWCF","PRECT","QRL","SWCF")

#  version with normalized parameters and added flnt and fsnt (for Covey fig...)
full_data.cc.norm=as.data.frame(cbind(params.norm,mns$AEROD_V,mns$CLDL,mns$FLNT,mns$FNET,mns$FSNT,mns$LWCF,mns$PRECT,mns$QRL,mns$SWCF,ssm,cs))
colnames(full_data.cc.norm)[10:18]=c("AEROD_V","CLDL","FLNT","FNET","FSNT","LWCF","PRECT","QRL","SWCF")

# version with unnormalized parameters and with defaults appended:
full_data.stand <- rbind(full_data, X.standard.append)

# versions with normalized parameters and with defaults appended:
full_data.stand.norm <- rbind(full_data.norm, X.standard.norm.append)

# Curt Covey fig versions (and with defaults appended):
full_data.cc.stand.norm <- rbind(full_data.cc.norm, X.standard.norm.append)

#==========================================

# ---------------------------------------
# Helper functions

f <- function(s){
  strsplit(s, split = "a.pt")[[1]][1]
}


open.field <- function(fn, var){

  # helper function to load a map of var from nc file

  nc <- open.nc(fn)
  nc.var <- var.get.nc(nc, var)
  nc.var
}

load.spatial.ens <- function(fn.list, var){

  # open all nc files in a list, vectorise, and concatenate to
  # an ensemble matrix, each row is a map

  field.list <- lapply(fn.list, FUN = open.field, var = var)

  out <- t(sapply(field.list,cbind)) # should do by columns
  out
}

remap.famous <- function(dat,longs,lats, shift = FALSE){

  # reshape a map in vector form so that image() like functions
  # will plot it correctly

  mat <- matrix(dat, nrow = length(longs), ncol = length(lats))[ ,length(lats):1]

  if(shift){

    block1.ix <- which(longs < shift)
    block2.ix <- which(longs > shift)

    mat.shift <- rbind(mat[ block2.ix, ], mat[block1.ix, ])

    out <- mat.shift
  }

  else{
    out <- mat
  }

  out
}

pc.project <- function(pca,scores.em,Z.em,scale){

  # project principal components

  num.pc <- dim(scores.em)[2]

  if (scale){
    anom <- ((pca$rotation[ ,1:num.pc] %*% t(scores.em))*pca$scale)
    anom.sd <- ((pca$rotation[ ,1:num.pc] %*% t(Z.em))*pca$scale)
  }

  else {
    anom <- pca$rotation[ ,1:num.pc] %*% t(scores.em)
    anom.sd <- pca$rotation[ ,1:num.pc] %*% t(Z.em)
  }

  tens <- t(anom + pca$center)

  return(list(tens = tens, anom.sd = anom.sd))
}


km.pc <- function(Y, X, newdata, num.pc, scale = FALSE, center = TRUE, type = "UK", ...){

  # Base function for emulation of high dimensional data
  # with PCA and Gaussian Process emulator

  if (class(Y)!= 'prcomp'){
    pca <- prcomp(Y,scale = scale, center = center)
  }

  else{
    pca <- Y
  }

  if(is.matrix(newdata)!= TRUE){
    print('matrixifying newdata')
    newdata <- matrix(newdata,nrow = 1)
  }

  scores.em <- matrix(nrow = dim(newdata)[1],ncol = num.pc)
  Z.em <- matrix(nrow = dim(newdata)[1],ncol = num.pc)

  for (i in 1:num.pc){

    # build the GP model

    fit <- km(design = X, response = pca$x[,i])
    pred <- predict(fit, newdata = newdata, type = type, ...)

    scores.em[ ,i] <- pred$mean
    Z.em[ ,i] <- pred$sd

  }

  proj = pc.project(pca, scores.em, Z.em, scale)

  return(list(tens = proj$tens,scores.em = scores.em,Z.em = Z.em,anom.sd = proj$anom.sd))
}


prop.thres <- function(x, thres, above = FALSE){

  # propotion of vector x below a threshold thres

  n <- length(x)

  if(above) bt <- length(x[x > thres])

  else bt <- length(x[x < thres])

  prop <- bt/n

  prop
}



# ---------------------------------------
# pallettes
rb <- brewer.pal(11, "RdBu")
ryg <- brewer.pal(11, "RdYlGn")
pbg <- brewer.pal(9, "PuBuGn")
bg <- brewer.pal(9, "BuGn")
yg <- brewer.pal(9, "YlGn")
byr <- rev(brewer.pal(11,'RdYlBu'))
br <- rev(rb)
blues <-  brewer.pal(9,'Blues')
rblues <-  rev(blues)

greens <-  brewer.pal(9,'Greens')
ygb <- brewer.pal(9, "YlGnBu")
brbg <- brewer.pal(11, "BrBG")
yob <- brewer.pal(9, "YlOrBr")
yor <- brewer.pal(9, "YlOrRd")

acc <- brewer.pal(8,'Paired')

col.amaz <- acc[1]
col.namerica <- acc[2]
col.seasia <- acc[3]
col.congo <- acc[4]
col.global <- acc[5]


pch.global <- 3
pch.amaz <- 1
pch.congo <- 2
pch.seasia <- 5
pch.namerica <- 4

# ---------------------------------------

# ---------------------------------------
# input space

inputs.set <- function(X, y, thres, obs, obs.sd = 0, disc = 0, disc.sd = 0, n = 100000, abt = FALSE){

  # find a set of inputs that are consistent with a particular
  # set of implausibility (either below or above)

  X.mins <- apply(X,2,min)
  X.maxes <- apply(X,2,max)

  X.unif <- samp.unif(n, mins = X.mins, maxes = X.maxes)
  colnames(X.unif) <- colnames(X)

  fit <- km(~., design = X, response = y, control = list(trace = FALSE))
  pred <- predict(fit, newdata = X.unif, type = 'UK')
  pred.impl <- impl(em = pred$mean, em.sd = pred$sd,
                    disc = disc, obs = obs, disc.sd = disc.sd, obs.sd = obs.sd)

  if(abt){
    # choose those above the threshold
    ix.bt <- pred.impl > thres
  }

  else{
    ix.bt <- pred.impl < thres
  }

  X.out <- X.unif[ix.bt, ]

  return(list(X.out = X.out, fit = fit, X.unif = X.unif, pred = pred,pred.impl = pred.impl))

}


dfunc.up <- function(x,y,...){
  require(MASS)
  require(RColorBrewer)

  #  br <- brewer.pal(9, 'Blues')
  br <- brewer.pal(9, 'Greys')
  # function for plotting 2d kernel density estimates in pairs() plot.
  kde <- kde2d(x,y)
  image(kde, col = br, add = TRUE)

}

dfunc.up.line <- function(x,y,...){
  # cgf mod to plot reference lines based on last row (default parameters)
  require(MASS)
  require(RColorBrewer)
  nt=length(x)

    br <- brewer.pal(9, 'Blues')
  #  br <- brewer.pal(9, 'Greys')
  # function for plotting 2d kernel density estimates in pairs() plot.
  kde <- kde2d(x[1:nt-1],y[1:nt-1])
  image(kde, col = br, add = TRUE,xlim=c(0,1),ylim=c(0,1))
  abline(v=x[nt],col="red",lwd=1.5)
  abline(h=y[nt],col="red",lwd=1.5)

}

dfunc.up.pt <- function(x,y,...){
  # cgf mod to plot reference points based on last row (default parameters)
  require(MASS)
  require(RColorBrewer)
  nt=length(x)

  br <- brewer.pal(9, 'Blues')

  # function for plotting 2d kernel density estimates in pairs() plot.
  kde <- kde2d(x[1:nt-1],y[1:nt-1])
  image(kde, col = br, add = TRUE,xlim=c(0,1),ylim=c(0,1))
  points(x[nt],y[nt],col="red",pch=19)

}

dfunc.image <- function(x,y,...){
  # produce paired contour images of CS by x/y in a pairs plot
  require(MASS)
  require(RColorBrewer)
  require(akima)
  # z-values: would be nice to find a way to generalize this, by passing it through panel
  yn<-normalize(y)
  xn<-normalize(x)
  z<-X.out$cs
  # colours now specified in colarr.

  # create contour surface for this pair
  fld<-interp(xn,yn,z)
  image(fld, col = colarr, add = TRUE,zlim=c(zl1,zl2),breaks=breaks,nlevel=nl)
}
