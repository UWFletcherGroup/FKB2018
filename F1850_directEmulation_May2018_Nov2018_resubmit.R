# ---------------------------------------
# cgf: April 20, 2018
# this version works with the revised F1850 simulations run by Bakr
# that have fixed a bug with the x parameter.
# The previous F1850 simulations did not perturb sulfate.
#
# Modifications: - variable names are now all consistent with CAM, so "aod" becomes "AEROD_V"
#                - dAOD is much more variable now that x1 is working correctly,
#                  so AOD threshold needs to be increased to 0.08 (max possible AOD = 0.2 (0.12 + 0.08))
#                  dAOD can also be (slightly) negative, so modify selection for abs(aod)
#                - threshold setting moved to common routine.
#
# May 4, 2018:
# re-ran n=350 training cases on Graham. Use those data here.
#
# Apr 25, 2018:
# use 10th/90th percentiles for CS threshold (taken from emulated output, and defined in common)
# also changed the x range for CS to include more range.
#
# Mar 21, 2018:
# cgf mod: use normalized predictors (full_data.norm)
#
# ---------------------------------------
#
# cgf: Run the emulator to predict CS and AOD and SS,
# use predict() to expand the probability space,
# then use SS to constrain parameter ranges.
#
# Adapted from Doug McNeall's script famous_direct.R
# Direct emulation of Forest fraction of the FAMOUS model

# -----------------------------------------------------------------
# 0. Packages, functions and data
# these are sourced from the F1850_common_May2018.R setup script
# -----------------------------------------------------------------

source('F1850_common_May2018.R')

# -----------------------------------------------------------------
# 1. Find a set of inputs consistent with SS > thresh
# - fit km() model for SS given inputs X, and predict (emulate) SS for full probability space
# - then select areas of prob space where SS > thresh,
# - and make pairs density plot of X.
# -----------------------------------------------------------------

# fit model

# define input parameters
X=full_data.norm[,1:9]

# fit GP/kriging model
fit <- km(~., design = X, response = full_data.norm[,"ssm"], control = list(trace = FALSE))

# set up full parameter space
n=100000 # number of prediction points

# ranges based on mins/maxs of normalized INPUTS
X.mins <- apply(X,2,min)
X.maxes <- apply(X,2,max)

# sample from 9-dimensional uniform distribution n times
X.unif <- samp.unif(n, mins = X.mins, maxes = X.maxes)
colnames(X.unif) <- colnames(X)
pred.ss <- predict(fit, newdata = X.unif, type = 'UK')

# select plausible data from pred
#ix.at<-pred.ss$mean > ss_thresh
#X.out <- X.unif[ix.at, ]
# append default parameters as bottom row
#X.out.stand <- rbind(X.out,params.stand.norm[length(params.stand.norm[,1]),])


# -----------------------------------------------------------------
# 2. Next, add an additional constraint not in the SS: dAOD < aod_thresh
# -----------------------------------------------------------------
# this is the threshold for dAOD: the acceptable difference in AOD between perturbed and default CAM.
# note that because of epsilon parameter, all perturbed cases have higher AOD than default, i.e. dAOD is always positive.
# build emulator for dAOD and emulate parameter space
fit <- km(~., design = X, response = full_data.norm[,"AEROD_V"], control = list(trace = FALSE))
pred.aod <- predict(fit, newdata = X.unif, type = 'UK')

# apply plausibility thresholds (for skill score and dAOD)
# select plausible data from pred.ss and pred.aod
X.out<-subset(X.unif,pred.ss$mean > ss_thresh & abs(pred.aod$mean) < aod_thresh)
# append default parameters as bottom row
X.out.stand <- rbind(X.out,params.stand.norm[length(params.stand.norm[,1]),])

print(paste("SS>thresh and AOD<thresh: ",100*dim(X.out)[1]/n," % of prob space."))

# 2D pairs plot of INPUTs for cases with SS>thresh and dAOD < aod_thresh:
pfile='best_inputs_SS_and_AOD.pdf'
file = paste(plotpath,pfile,sep="")
pdf(width = 7, height = 7, file = file)
par(lab=c(3,3,7)) # constrain 3 labels on each plot
# in response to Reviewer #2, removed the red line, replaced by point.
pairs(X.out.stand, panel = dfunc.up.pt, gap = 0.75, upper.panel = NULL,xlim=c(0,1),ylim=c(0,1))

dev.off()

# Google Drive: clean up old version, and save new one:
# drive_trash(paste(gpth,pfile,sep=""))
# drive_upload(file,path=gpth)

# -----------------------------------------------------------------
# 2b. Next, show pairs density for high (>0.5) and low (<0.5) CS from constrained sample:
# -----------------------------------------------------------------

# build emulator for climate sensitivity (cs) and emulate parameter space
fit <- km(~., design = X, response = full_data.norm[,"cs"], control = list(trace = FALSE))
pred.cs <- predict(fit, newdata = X.unif, type = 'UK')

# apply plausibility thresholds (for ss and dAOD), and high/low cs thresholds
# select plausible data from pred.ss and pred.aod
X.out.hcs<-subset(X.unif,pred.ss$mean > ss_thresh & abs(pred.aod$mean) < aod_thresh & pred.cs$mean > cs_uthresh)
X.out.lcs<-subset(X.unif,pred.ss$mean > ss_thresh & abs(pred.aod$mean) < aod_thresh & pred.cs$mean < cs_lthresh)
# append default parameters as bottom row
X.out.hcs.stand <- rbind(X.out.hcs,params.stand.norm[length(params.stand.norm[,1]),])
X.out.lcs.stand <- rbind(X.out.lcs,params.stand.norm[length(params.stand.norm[,1]),])

print(paste("SS>thresh and AOD<thresh, High CS: ",100*dim(X.out.hcs)[1]/n," % of prob space."))
print(paste("SS>thresh and AOD<thresh, Low CS: ",100*dim(X.out.lcs)[1]/n," % of prob space."))

#
# Nov 2018: attempt to show HighCS and LowCS on the same plot
# first rbind the two sets of data together:
pfile='best_inputs_SS_and_AOD_High_AND_LowCS_Nov2018resubmit.pdf'
file = paste(plotpath,pfile,sep="")
pdf(width = 7, height = 7, file = file)
XC<-X.out.lcs.stand
YC<-X.out.hcs.stand
XY<-rbind(XC,YC)
par(lab=c(2,2,7))
cx=1.0
pairs(XY,
      lower.panel=function(x, y, ...) {
        Xx <- x[seq_len(nrow(XC))] # corresponds to X subset
        Xy <- y[seq_len(nrow(XC))] # corresponds to X subset

        nt=length(Xx)
        br <- brewer.pal(9, 'Blues')

        # function for plotting 2d kernel density estimates in pairs() plot.
        kde <- kde2d(Xx[1:nt-1],Xy[1:nt-1])
        image(kde, col = br, add = TRUE)
        points(Xx[nt],Xy[nt],col="black",pch=19)

        if(par('mfg')[2] == 1) axis(2,cex.axis=cx) # if left plot, add left axis
        if(par('mfg')[1] == ncol(XC)) axis(1,cex.axis=cx) # if bottom plot add bottom axis
      },
      upper.panel=function(x, y, ...) {
        Yx <- x[(nrow(XC) + 1):length(x)] # Y subset
        Yy <- y[(nrow(XC) + 1):length(y)] # Y subset

        nt=length(Yx)
        br <- brewer.pal(9, 'Reds')

        # function for plotting 2d kernel density estimates in pairs() plot.
        kde <- kde2d(Yx[1:nt-1],Yy[1:nt-1])
        image(kde, col = br, add = TRUE)
        points(Yx[nt],Yy[nt],col="black",pch=19)

        if(par('mfg')[2] == ncol(YC)) axis(4,cex.axis=cx) # if right plot, add right axis
        if(par('mfg')[1] == 1) axis(3,cex.axis=cx) # if top plot, add top axis
      },xlim=c(0,1),ylim=c(0,1),gap=0.75,cex.axis=cx) # move the default tick labels off the plot
dev.off()
# Google Drive: clean up old version, and save new one:
# drive_trash(paste(gpth,pfile,sep=""))
# drive_upload(file,path=gpth)


# -----------------------------------------------------------------
# 3. Next, plot histograms for OUTPUTs for this subset of "best cases"
# -----------------------------------------------------------------

pfile='OUTPUT_histograms_bestCases_SS_and_AOD_Nov2018resubmit.pdf'
file = paste(plotpath,pfile,sep="")
pdf(file = file, width = 5, height = 7)
par(mfrow= c(4,2) , las = 1, mar = c(4,3,1,2))

# cgf Aug 16, 2018: try including it again...
usedata<-full_data.norm
varlabs=names(usedata)
varlabs[c(10,17,18)]=c("AOD","SS","lambda")
unitstrs=c("(%)","(Wm-2)","(Wm-2)","(mm/day)","(Wm-2/s)","(Wm-2)","(unitless)","(K/Wm-2)")
# different ranges for each variable, same number of bins...
br <- 51
# loop over OUTPUTS, build emulator and emulate parameter space, then subset by SS and AOD

for (i in 1:8){
 print(paste("processing",names(usedata)[i+10],sep=""))
  # indexing starts at 11 for CLDL
 fit <- km(~., design = X, response = usedata[,i+10], control = list(trace = FALSE))
 pred <- predict(fit, newdata = X.unif, type = 'UK')
 # get density for the training set (CAM4) which will be overlaid as a black line
 histt <- density(usedata[,i+10], n=512)
 # define two density objects for each OUTPUT (one is all cases, the other is subset)
 hist1 <- density(pred$mean, n=256)
 hist <- density(subset(pred$mean,pred.ss$mean > ss_thresh & abs(pred.aod$mean) < aod_thresh),n=256)
 yl=c(0,max(hist$y,hist1$y,histt$y))
 x1=min(hist$x,hist1$x,histt$x)
 x2=max(abs(hist$x),abs(hist1$x),abs(histt$x))
 xl=c(x1,x2)
 xtit=paste(varlabs[i+10]," ",unitstrs[i],sep="")
 if(i==8){xtit=expression(paste(lambda," (K / Wm-2)",sep=""))}
 plot(hist1$x,hist1$y,xlim = xl, main = '',xlab = xtit, ylab = '', axes = T,lwd=2,ylim=yl,col="darkgrey",type="h")
 lines(histt$x,histt$y,col="black",lwd=2)
 lines(hist$x,hist$y,col="red",lwd=2)
 legend('topleft',legend = paste("(",letters[i],")",sep=""),bty='n',cex=1.2,xjust=1)
 # cgf mod (Aug 16, 2018): add line to other plots at 0 (default model)
 if(i<7){abline(v=0,col="darkgreen",lwd=2)}
 # cgf mod: add line to SS plot showing default model SS = 1.0
 if(i==7){abline(v=1.0,col="darkgreen",lwd=2)}
 # cgf mod: add line to CS plot showing default model CS = 0.45
 if(i==8){abline(v=0.45,col="darkgreen",lwd=2)}
}

dev.off()
# Google Drive: clean up old version, and save new one:
# drive_trash(paste(gpth,pfile,sep=""))
# drive_upload(file,path=gpth)
