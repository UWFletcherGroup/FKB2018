# ---------------------------------------
# Disclaimer: this script is adapted from Doug McNeall's script famous_direct.R,
# for direct emulation of Forest fraction of the FAMOUS climate model.
# Source is located in the famous-git directory at github.com/dougmcneall/famous-git
#
# ---------------------------------------
# cgf: April 20, 2018
# this version works with the revised F1850 simulations run by Bakr
# that have fixed a bug with the x parameter.
# The previous F1850 simulations did not perturb sulfate.
#
# Modifications: added threshold-setting in common, so it can be used by multiple routines consistently.
#
# May 4, 2018: re-ran n=350 training cases on Graham. Use those data here.
#
# May 24, 2018: fixed non-reproducibility bug when running validation twice.
# ---------------------------------------
#
# cgf: Run the emulator to predict CS and AOD and SS,
# use predict() to expand the probability space,
# then use SS to constrain parameter ranges.
#
# Adapted from Doug McNeall's script famous_direct.R
# Direct emulation of Forest fraction of the FAMOUS
# model


# -----------------------------------------------------------------
# 0. Packages, functions and data
# these are sourced from the F1850_common_May2018.R setup script
# -----------------------------------------------------------------
library(plotrix)
source('F1850_common_May2018.R')

# -----------------------------------------------------------------
# 1. Compare GP (kriging) with MLR.
# -----------------------------------------------------------------

varlabs=c("SS","CS")

pfile='F1850_validateEmulation_SS_AOD_CS_compareMLR_Nov2018.pdf'
file = paste(plotpath,pfile,sep="")

pdf(file = file, width = 5, height = 5)
par(mfrow= c(2,2),mar = c(4,4,0.5,1)) #, las = 1, mar = c(4,4,0.5,1), oma = c(5,2,0,0), cex.axis = 1.0)
# do training/test splits? No, I think this is done within LOO.km function...

# build emulator on training data; the output is prediction for each data point.
varnames=c("ssm","cs") #"AEROD_V",

let=c("(a) ","(b) ","(c) ")
for(i in 1:length(varnames)){
  varname=varnames[i]
  print(paste("Processing: ",varname))
  fit<-km(~., design = X, response = full_data[,varname], control = list(trace = FALSE))
  pred<-leaveOneOut.km(fit,trend.reestim = F,type='UK')

  #
  # make CI plot:
  xl=c(0.2,1.0)
  xla=""
  varlab=varlabs[i]
  if(varname=="AEROD_V"){xl=c(0,0.15)}
  if(varname=="cs"){
    xl=c(0.3,0.7)
    xla="CAM4"
    varlab=expression(paste(lambda))
  }
  x=full_data[,varname]
  y=pred$mean
  RMSE_test=round(sqrt(mean((y-x)^2)),3)
  corr=round(cor(y,x),3)
  plotCI(x,y,xlim=xl,ylim=xl,ui=y+(2*pred$sd),li=y-(2*pred$sd),ylab="GP Emulator",xlab=xla,pch=21,col="black",scol="darkgrey",cex=0.25, pt.bg=par("bg"))
  abline(a=0,b=1)
  pform=as.formula(paste("y~",varname))
  model1=lm(pform,data=full_data)
  abline(model1,col="red")
  legend('topleft',legend=varlab,bty='n',cex=1.1)
  legend('bottomright',legend=c(paste("RMSE=",RMSE_test),paste(" r=",corr)),bty='n',xjust=1)

  # simple test using MLR (tested, we don't really need LOOCV here...)
  form=as.formula(paste(varname,"~x1+x2+x3+x4+x5+x6+x7+x8+x9",sep=""))
  fit.mlr<-lm(form,data=full_data)
  # get 95% prediction intervals for MLR
  fit.mlr.pi=predict(fit.mlr, newdata=full_data, interval="predict")
  x=full_data[,varname]
  y=fit.mlr$fitted.values
  RMSE_test=round(sqrt(mean((y-x)^2)),3)
  corr=round(cor(y,x),3)
  plotCI(x,y,xlim=xl,ylim=xl,ui=fit.mlr.pi[,3],li=fit.mlr.pi[,2],ylab="MLR Emulator",xlab=xla,pch=1,col="black",scol="darkgrey",cex=0.25)
  abline(a=0,b=1)
  pform=as.formula(paste("y~",varname))
  model1=lm(pform,data=full_data)
  abline(model1,col="red")
  legend('topleft',legend=varlab,bty='n',cex=1.1)
  legend('bottomright',legend=c(paste("RMSE=",RMSE_test),paste(" r=",corr)),bty='n',xjust=1)


}

dev.off()

# Google Drive: clean up old version, and save new one:
# drive_trash(paste(gpth,pfile,sep=""))
# drive_upload(file,path=gpth)
