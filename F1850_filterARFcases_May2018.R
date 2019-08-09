#
# ---------------------------------------
# cgf: April 20, 2018
# this version works with the revised F1850 simulations run by Bakr
# that have fixed a bug with the x parameter.
# The previous F1850 simulations did not perturb sulfate.
#
# Modifications: added threshold-setting in common, so it can be used by multiple routines consistently.
# ---------------------------------------
#
# Apr 26: Added code to make a pairs image plot of CS by input parameters
#
# May 4, 2018: re-ran n=350 training cases on Graham. Use those data here.
# ---------------------------------------
# ---------------------------------------
#   IMPORTANT NOTE: this requires using normalized inputs, whereas the subsetting relies on non-normalized data.
#   So you'll need to change full_data.norm to full_data to do the subsetting.
# ---------------------------------------
# ---------------------------------------
#
#
# cgf: Run the emulator for outputs CS, AOD and SS based on Input parameters,
# then use predict() to expand the probability space (to n=100,000)
# then use SS and AOD, and aerosol parameters, to extract ARF "scenario" cases:
#   Sh.Bh, Sm.Bm, Sl.Bl, Sh.Bl and Sl.Bh
# while also looking for end-members in Cess CS.
#
#
# Adapted from Doug McNeall's script famous_direct.R
# Direct emulation of Forest fraction of the FAMOUS
# model
#
# NOTE: not using normalized predictors, for simplicity.

# -----------------------------------------------------------------
# 0. Packages, functions and data
# these are sourced from the F1850_common_May2018.R setup script
# -----------------------------------------------------------------

source('F1850_common_May2018.R')

# -----------------------------------------------------------------
# 1. Define thresholds for subsetting parameter space:
# -----------------------------------------------------------------
# thresholds are now set in F1850_common_May2018.R, not here

# -----------------------------------------------------------------
# 2. Constrain on SS > ss_thresh and dAOD < aod_thresh
# -----------------------------------------------------------------
# fit models: first emulate SS

# Change this to full_data[,1:9] to run the subsetting
# ---------------------------------------
#X=full_data.norm[,1:9]  # normalized data required to make the image plots of CS by param.
X=full_data[,1:9]
# ---------------------------------------
# build emulator for skill score
fit.ss <- km(~., design = X, response = full_data[,"ssm"], control = list(trace = FALSE))
# build emulator for dAOD
fit.aod <- km(~., design = X, response = full_data[,"AEROD_V"], control = list(trace = FALSE))
# build emulator for Cess CS
fit.cs <- km(~., design = X, response = full_data[,"cs"], control = list(trace = FALSE))

# set up full parameter space for emulation: 10^5 samples gives reasonable resolution:
n=100000
# input parameter ranges based on mins/maxs of X
X.mins <- apply(X,2,min)
X.maxes <- apply(X,2,max)
# sample from uniform distribution n times
X.unif <- samp.unif(n, mins = X.mins, maxes = X.maxes)
colnames(X.unif) <- colnames(X)
# emulate AOD
print("Emulating AOD...")
pred.aod <- predict(fit.aod, newdata = X.unif, type = 'UK')
# emulate SS
print("Emulating SS...")
pred.ss <- predict(fit.ss, newdata = X.unif, type = 'UK')
# emulate CS
print("Emulating Cess CS...")
pred.cs <- predict(fit.cs, newdata = X.unif, type = 'UK')

# for each output, select plausible data from X parameters, pred.ss, pred.aod and pred.cs
X.out<-as.data.frame(subset(X.unif,pred.ss$mean > ss_thresh & abs(pred.aod$mean) < aod_thresh))
cs.out<-as.data.frame(subset(pred.cs$mean,pred.ss$mean > ss_thresh & abs(pred.aod$mean) < aod_thresh))
aod.out<-as.data.frame(subset(pred.aod$mean,pred.ss$mean > ss_thresh & abs(pred.aod$mean) < aod_thresh))
ss.out<-as.data.frame(subset(pred.ss$mean,pred.ss$mean > ss_thresh & abs(pred.aod$mean) < aod_thresh))
names(cs.out)="cs"
names(ss.out)="ss"
names(aod.out)="AEROD_V"
X.out=cbind(X.out,aod.out,ss.out,cs.out)

print(paste("SS>thresh and AOD<thresh: ",100*dim(X.out)[1]/n," % of prob space."))


# -----------------------------------------------------------------
# 3. extract ARF "scenario" cases from X.out:
#   Sh.Bh, Sm.Bm, Sl.Bl, Sh.Bl and Sl.Bh
# while also looking for end-members in Cess CS.
# -----------------------------------------------------------------

DefaultVals <- c(0.88, 14.00, 1800, 0.50, 3600)
DefaultLowerBounds <- c(0.8,8.4,900,0.50,1800)
DefaultUpperBounds <- c(0.99,19.6,14400,0.85,28800)

# Black Carbon and Sulphate Definitions
# sulfate hygroscopic (x1)
# cgf mod: new thresholds for SO4 based on (now correct) effect of x1:
highS <- 0.5
midS <- 0.25
lowS <- 0.1

# BC gamma (x2)
highG <- 0.80
midG <- 0.5
lowG <- 0.2
# BC delta (x3)
highD <- 20
midD <- 10
lowD <- 5
# BC epsilon (x4)
lowE <- 5
midE <- 10
highE <- 15

# Filter high BC Scenarios
Bh <- X.out[(X.out$x2>highG) & (X.out$x3>highD) & (X.out$x4<lowE),]
ShBh<- Bh[Bh$x1 > highS,]
SlBh<- Bh[Bh$x1 < lowS,]

# Filter low BC Scenarios
Bl <- X.out[(X.out$x2<lowG) & (X.out$x3<lowD) & (X.out$x4>highE),]
ShBl<- Bl[Bl$x1 > highS,]
SlBl<- Bl[Bl$x1 < lowS,]

# Filter mid BC Scenarios (use between() function from dplyr)
Bm <- X.out[between(X.out$x2,midG-0.15,midG+0.15) & between(X.out$x3,midD-2,midD+2) & between(X.out$x4,midE-2,midE+2),]
SmBm<- Bm[between(Bm$x1,midS-0.15,midS+0.15),]

# -----------------------------------------------------------------
# Select cases based on CS: highest for Sh.B* and lowest for Sl.B* (F = "final" selection)
# -----------------------------------------------------------------
outtable=list()
outtable=rbind(outtable,ShBl[ShBl$cs==max(ShBl$cs),])
outtable=rbind(outtable,ShBh[ShBh$cs==max(ShBh$cs),])
outtable=rbind(outtable,SmBm[which.min(abs(SmBm$cs-median(SmBm$cs))),])
outtable=rbind(outtable,SlBl[SlBl$cs==min(SlBl$cs),])
outtable=rbind(outtable,SlBh[SlBh$cs==min(SlBh$cs),])
rownames(outtable)=c("ShBh","ShBl","SmBm","SlBh","SlBl")
#
file="F1850_filterARFcases_finalScenarios_May2018.csv"
write.csv(outtable,file=paste(dd,file,sep=""))
# Google Drive: clean up old version, and save new one:
# drive_trash(paste(gpth,file,sep=""))
# drive_upload(paste(dd,file,sep=""),path=gpth)




# -----------------------------------------------------------------
# now make a scatter plot for each subset: SS vs CS
# -----------------------------------------------------------------
pfile='F1850_ARFscenarios_SS_vs_CS.pdf'
pdf(file = paste(plotpath,pfile,sep=""))
par(mfrow= c(3,2))
xlm=c(0.85,1.00)
ylm=c(0.25,0.75)
plot(ShBh$ss,ShBh$cs,pch=19,main="Sh.Bh",xlim=xlm,ylim=ylm)
plot(ShBl$ss,ShBl$cs,pch=19,main="Sh.Bl",xlim=xlm,ylim=ylm)
plot(SmBm$ss,SmBm$cs,pch=19,main="Sm.Bm",xlim=xlm,ylim=ylm)
plot(SlBh$ss,SlBh$cs,pch=19,main="Sl.Bh",xlim=xlm,ylim=ylm)
plot(SlBl$ss,SlBl$cs,pch=19,main="Sl.Bl",xlim=xlm,ylim=ylm)

dev.off()
# Google Drive: clean up old version, and save new one:
# drive_trash(paste(gpth,pfile,sep=""))
# drive_upload(paste(plotpath,pfile,sep=""),path=gpth)

#
# -----------------------------------------------------------------
# now make a density plot showing spread of all plausible cases, and then spread of ARF subsets.
# -----------------------------------------------------------------
pfile='F1850_ARFscenarios_CessCSdistributions.pdf'
pdf(file = paste(plotpath,pfile,sep="")) #, width = 7, height = 7)
layout(1)
xl=c(0.35,0.65)
# first plot density for all plausible cases (in X.out)
plot(density(X.out$cs)$x,density(X.out$cs)$y,col="darkgrey",type='h',xlim=c(0.3,0.65),ylim=c(0,10),xlab="Cess CS (K / Wm**2)",ylab="Density")
# now add coloured lines to show the different cases: Sh are red, Sl are blue
lines(density(ShBl$cs)$x,density(ShBl$cs)$y,col="red",lwd=3)
lines(density(ShBh$cs)$x,density(ShBh$cs)$y,col="red",lwd=2,lty=2)
lines(density(SmBm$cs)$x,density(SmBm$cs)$y,col="darkgreen",lwd=3)
lines(density(SlBl$cs)$x,density(SlBl$cs)$y,col="blue",lwd=3)
lines(density(SlBh$cs)$x,density(SlBh$cs)$y,col="blue",lwd=2,lty=2)
legend('topright', legend = c('Sh.Bl', 'Sh.Bh','Sm.Bm','Sl.Bl','Sl.Bh'),
       col=c("red","red","darkgreen","blue","blue"),lty=c(1,2,1,1,2),lwd=2,
       text.col = 'black', cex = 1.3)
dev.off()
# Google Drive: clean up old version, and save new one:
# drive_trash(paste(gpth,pfile,sep=""))
# drive_upload(file,path=gpth)

# Apr 25, 2018
# cgf: new method here is to extract parameters associated with percentiles of the CS distribution:
outdat=list()
for(i in c(0.0,0.1,0.25,0.5,0.75,0.9,1.0)){
thisrow=X.out[X.out$cs==quantile(X.out$cs,i,type=1),]
  outdat=rbind(outdat,thisrow)
}
# write outdat to file, and upload to Google Drive
file="F1850_May2018_emulatedInputsOutputs_for_CessCSquantiles.csv"
write.csv(outdat,file=paste(dd,file,sep=""))
# Google Drive: clean up old version, and save new one:
# drive_trash(paste(gpth,file,sep=""))
# drive_upload(paste(dd,file,sep=""),path=gpth)

#
#
# Apr 26, 2018: image plot showing CS as function of input params
# NEEDS NORMALIZED DATA, so re-run the emulator here:
X=full_data.norm[,1:9]  # normalized data required to make the image plots of CS by param.
# ---------------------------------------
fit.ss <- km(~., design = X, response = full_data[,"ssm"], control = list(trace = FALSE))
# build emulator for dAOD
fit.aod <- km(~., design = X, response = full_data[,"AEROD_V"], control = list(trace = FALSE))
# build emulator for Cess CS
fit.cs <- km(~., design = X, response = full_data[,"cs"], control = list(trace = FALSE))

# set up full parameter space for emulation: 10^5 samples gives reasonable resolution:
n=100000
# input parameter ranges based on mins/maxs of X
X.mins <- apply(X,2,min)
X.maxes <- apply(X,2,max)
# sample from uniform distribution n times
X.unif <- samp.unif(n, mins = X.mins, maxes = X.maxes)
colnames(X.unif) <- colnames(X)
# emulate AOD
print("Emulating AOD...")
pred.aod <- predict(fit.aod, newdata = X.unif, type = 'UK')
# emulate SS
print("Emulating SS...")
pred.ss <- predict(fit.ss, newdata = X.unif, type = 'UK')
# emulate CS
print("Emulating Cess CS...")
pred.cs <- predict(fit.cs, newdata = X.unif, type = 'UK')

#
# select plausible data from X parameters, pred.ss, pred.aod and pred.cs
X.out<-as.data.frame(subset(X.unif,pred.ss$mean > ss_thresh & abs(pred.aod$mean) < aod_thresh))
cs.out<-as.data.frame(subset(pred.cs$mean,pred.ss$mean > ss_thresh & abs(pred.aod$mean) < aod_thresh))
aod.out<-as.data.frame(subset(pred.aod$mean,pred.ss$mean > ss_thresh & abs(pred.aod$mean) < aod_thresh))
ss.out<-as.data.frame(subset(pred.ss$mean,pred.ss$mean > ss_thresh & abs(pred.aod$mean) < aod_thresh))
names(cs.out)="cs"
names(ss.out)="ss"
names(aod.out)="AEROD_V"
X.out=cbind(X.out,aod.out,ss.out,cs.out)

print(paste("SS>thresh and AOD<thresh: ",100*dim(X.out)[1]/n," % of prob space."))

# 2D pairs plot of INPUTs for plausible CS:
pfile='CS_for_Inputs_Pairs_May2018_Nov2018resubmit.pdf'
file = paste(plotpath,pfile,sep="")
pdf(width = 7, height = 7, file=file)
# force labels
par(lab=c(3,3,6))
# shading range for CS = zl1-->zl2:
zl1=0.40
zl2=0.60
breaks=seq(zl1,zl2,0.025)
nl=length(breaks)-1 # number of colours
colarr <- rev(brewer.pal(nl, 'RdBu'))

pairs(X.out[,1:9], panel = dfunc.image, gap = 0.75, upper.panel = NULL,xlim=c(0,1),ylim=c(0,1))
# Nov 2018: add legend bar
require(MASS)
require(RColorBrewer)
require(akima)
# z-values: would be nice to find a way to generalize this, by passing it through panel
yn<-normalize(X.out[,8])
xn<-normalize(X.out[,9])
z<-X.out$cs

br <- rev(brewer.pal(nl, 'RdYlBu'))
# create contour surface for this pair
fld<-interp(xn,yn,z)
image.plot(fld, col = colarr, zlim=c(zl1,zl2),legend.only = T,legend.shrink = 0.66,legend.lab = expression(paste(lambda," (K / Wm-2)",sep="")),legend.line = -2.1,nlevel=nl,breaks=breaks)
dev.off()
# Google Drive: clean up old version, and save new one:
# drive_trash(paste(gpth,pfile,sep=""))
# drive_upload(file,path=gpth)


stop() # stop script execution here, unless black carbon csv files desired

# Plot distribution of Cess CS in different scenarios
  hBChShSS <- hBChS[hBChS$SS >= minSS,]
  hBClShSS <- hBClS[hBClS$SS >= minSS,]
  lBChShSS <- lBChS[lBChS$SS >= minSS,]
  lBClShSS <- lBClS[lBClS$SS >= minSS,]
  hBChSSScount <- length(hBChShSS[["SS"]])
  hBClSSScount <- length(hBClShSS[["SS"]])
  lBChSSScount <- length(lBChShSS[["SS"]])
  lBClSSScount <- length(lBClShSS[["SS"]])

  text(x=7, y=5.5, labels="# High SS Cases within Climate Scenario Definition")
  text(x=7, y=5.0, labels=paste("hBChS # cases: ", hBChSSScount))
  text(x=7, y=4.5, labels=paste("hBClS # cases: ", hBClSSScount))
  text(x=7, y=4.0, labels=paste("lBChS # cases: ", lBChSSScount))
  text(x=7, y=3.5, labels=paste("lBClS # cases: ", lBClSSScount))
  text(x=7, y=2.5, labels=paste("High SS cases defined as cases with SS >", minSS))

  # Boxplot of Full Scenarios
  boxplot(x = c(hBChS["SS"], hBClS["SS"], lBChS["SS"], lBChS["SS"]),
          names = c("hBChS", "hBClS", "lBChS", "lBClS"), range=0,
          xlab="Aerosol Climate Scenarios", ylab="Skill Score",
          main="Spread of Skill Score in Different Aerosol Scenarios")

  # Boxplot of High Skill Score within the Full Scenarios
  boxplot(x = c(hBChShSS["SS"], hBClShSS["SS"], lBChShSS["SS"], lBChShSS["SS"]),
          names = c("hBChS", "hBClS", "lBChS", "lBClS"), range=0,
          xlab="Aerosol Climate Scenarios", ylab="Skill Score",
          main=paste("Spread of Skill Score in Different Aerosol Scenarios\n with SS >=", minSS))

  # Boxplots of parameter values of the High Skill Score cases
  for(x in names(paramDescriptions)){
    boxplot(hBChShSS[[x]],hBClShSS[[x]],lBChShSS[[x]],lBClShSS[[x]],
            names = c("hBChS", "hBClS", "lBChS", "lBClS"),
            range = 0, xlab = "Aerosol Climate Scenarios", ylab = paramDescriptions[x],
            main = paste(x, "Parameter Values in High SS Aerosol Climate Scenarios"))
    points(x=1, y= mean(hBChShSS[[x]]), col="green", pch=19)
    points(x=2, y= mean(hBClShSS[[x]]), col="green", pch=19)
    points(x=3, y= mean(lBChShSS[[x]]), col="green", pch=19)
    points(x=4, y= mean(lBClShSS[[x]]), col="green", pch=19)
    abline(h=DefaultVals[x], col="red")
  }

# Writing high Skill Score cases
write.csv(hBChShSS, "hBChS.csv")
write.csv(hBClShSS, "hBClS.csv")
write.csv(lBChShSS, "lBChS.csv")
write.csv(lBClShSS, "lBClS.csv")

# Summary Table
hBChSMax <- c(max(hBChShSS["x5"]), max(hBChShSS["x6"]), max(hBChShSS["x7"]), max(hBChShSS["x8"]), max(hBChShSS["x9"]))
hBChSMin <- c(min(hBChShSS["x5"]), min(hBChShSS["x6"]), min(hBChShSS["x7"]), min(hBChShSS["x8"]), min(hBChShSS["x9"]))
hBChSMean <- c(mean(hBChShSS[["x5"]]), mean(hBChShSS[["x6"]]), mean(hBChShSS[["x7"]]), mean(hBChShSS[["x8"]]), mean(hBChShSS[["x9"]]))
hBClSMax <- c(max(hBClShSS["x5"]), max(hBClShSS["x6"]), max(hBClShSS["x7"]), max(hBClShSS["x8"]), max(hBClShSS["x9"]))
hBClSMin <- c(min(hBClShSS["x5"]), min(hBClShSS["x6"]), min(hBClShSS["x7"]), min(hBClShSS["x8"]), min(hBClShSS["x9"]))
hBClSMean <- c(mean(hBClShSS[["x5"]]), mean(hBClShSS[["x6"]]), mean(hBClShSS[["x7"]]), mean(hBClShSS[["x8"]]), mean(hBClShSS[["x9"]]))
lBChSMax <- c(max(lBChShSS["x5"]), max(lBChShSS["x6"]), max(lBChShSS["x7"]), max(lBChShSS["x8"]), max(lBChShSS["x9"]))
lBChSMin <- c(min(lBChShSS["x5"]), min(lBChShSS["x6"]), min(lBChShSS["x7"]), min(lBChShSS["x8"]), min(lBChShSS["x9"]))
lBChSMean <- c(mean(lBChShSS[["x5"]]), mean(lBChShSS[["x6"]]), mean(lBChShSS[["x7"]]), mean(lBChShSS[["x8"]]), mean(lBChShSS[["x9"]]))
lBClSMax <- c(max(lBClShSS["x5"]), max(lBClShSS["x6"]), max(lBClShSS["x7"]), max(lBClShSS["x8"]), max(lBClShSS["x9"]))
lBClSMin <- c(min(lBClShSS["x5"]), min(lBClShSS["x6"]), min(lBClShSS["x7"]), min(lBClShSS["x8"]), min(lBClShSS["x9"]))
lBClSMean <- c(mean(lBClShSS[["x5"]]), mean(lBClShSS[["x6"]]), mean(lBClShSS[["x7"]]), mean(lBClShSS[["x8"]]), mean(lBClShSS[["x9"]]))

names(hBChSMax) <- paramHeads
names(hBChSMin) <- paramHeads
names(hBClSMax) <- paramHeads
names(hBClSMin) <- paramHeads
names(lBChSMax) <- paramHeads
names(lBChSMin) <- paramHeads
names(lBClSMax) <- paramHeads
names(lBClSMin) <- paramHeads

sumTab <- rbind(paramDescriptions, DefaultUpperBounds, DefaultLowerBounds, DefaultVals,
                hBChSMax, hBChSMin, hBChSMean, hBClSMax, hBClSMin, hBClSMean, lBChSMax, lBChSMin, lBChSMean, lBClSMax, lBClSMin, lBClSMean)
