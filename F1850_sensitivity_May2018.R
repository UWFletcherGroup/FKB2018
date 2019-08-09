# ---------------------------------------
# cgf: April 20, 2018
# this version works with the revised F1850 simulations run by Bakr
# that have fixed a bug with the x parameter.
# The previous F1850 simulations did not perturb sulfate.
#
# Modifications: variable names are now all consistent with CAM, so "aod" becomes "AEROD_V"
#
# May 4, 2018: re-ran n=350 training cases on Graham. Use those data here.
# ---------------------------------------
#
# cgf: Run the emulator and conduct variance-based sensitivity analysis.
#
# Adapted from Doug McNeall's script famous_sensitivity.R
#

# -----------------------------------------------------------------
# 0. Packages, functions and data
# -----------------------------------------------------------------

source('F1850_common_May2018.R')

# -----------------------------------------------------------------
# 1. Sensitivity analysis using the extended FAST algorithm of Saltelli et al,
# using the R package "sensitivity".
# ------------------------------------------------------------------------------

#
# We should check whether we need to first normalize the predictors.
X=full_data.norm[,1:9]
fit.ssm <- km(~., design = X, response = full_data.norm[,"ssm"], control = list(trace = FALSE))
fit.cs <- km(~., design = X, response = full_data.norm[,"cs"], control = list(trace = FALSE))
fit.aod <- km(~., design = X, response = full_data.norm[,"AEROD_V"], control = list(trace = FALSE))
#fit.seasia.norm <- km(~.,design = X.norm, response = SEASIA_MOD_FRAC)
#fit.congo.norm <- km(~.,design = X.norm, response = CONGO_MOD_FRAC)
#fit.namerica.norm <- km(~.,design = X.norm, response = NAMERICA_MOD_FRAC)
#fit.global.norm <- km(~.,design = X.norm, response = GLOB_MOD_FRAC)

# generate the design to run the emulator at, using fast99
x <- fast99(model = NULL, factors = colnames(X), n = 1000,
            q = "qunif", q.arg = list(min = 0, max = 1))

# run the emulator at the sensitivity analysis design points
fast.pred.ssm <- predict(fit.ssm, newdata = x$X, type = 'UK')
fast.pred.cs <- predict(fit.cs, newdata = x$X, type = 'UK')
fast.pred.aod <- predict(fit.aod, newdata = x$X, type = 'UK')


# Calculate the sensitivity indices
fast.ssm <- tell(x, fast.pred.ssm$mean)
fast.cs <- tell(x, fast.pred.cs$mean)
fast.aod <- tell(x, fast.pred.aod$mean)

bp.convert <- function(fastmodel){
  # get the FAST summary into an easier format for barplot
  fast.summ <- print(fastmodel)
  fast.diff <- fast.summ[ ,2] - fast.summ[ ,1]
  fast.bp <- t(cbind(fast.summ[ ,1], fast.diff))
  fast.bp
}

# Plot the sensitivity indices
#x11(width = 5.8, height = 10)
pfile='FAST_sensitivityHistograms_Nov2018.pdf'
file = paste(plotpath,pfile,sep="")

pdf(width = 5, height = 7, file = file)
par(mfrow = c(3,1), mar = c(4,4,0,2), las = 1, oma = c(5,4,1,2), fg = 'darkgrey', xaxs = 'i', cex.axis = 1.1)
barplot(bp.convert(fast.aod), ylim = c(0,1), col = c('blue', 'lightgrey'),border='black',ylab='Fraction of total variance')
mtext(side = 3, adj = 0, line = -0.5, text = '(a)', cex = 1.2, col = 'black')
legend('top', legend = c('Main effect', 'Interaction'),
       fill = c('blue', 'lightgrey'), bty = 'n', text.col = 'black', cex = 1.3)

barplot(bp.convert(fast.ssm), ylim = c(0,1), col = c('blue', 'lightgrey'),border='black',ylab='Fraction of total variance')
mtext(side = 3, adj = 0, line = -0.5, text = '(b)', cex = 1.2, col = 'black')
#legend('top', legend = c('Main effect', 'Interaction'),
#       fill = c('blue', 'lightgrey'), bty = 'n', text.col = 'black', cex = 1.3)

barplot(bp.convert(fast.cs), ylim = c(0,1), col = c('blue', 'lightgrey'),border='black',ylab='Fraction of total variance',xlab='Input parameter')
mtext(side = 3, adj = 0, line = -0.5, text = '(c)', cex = 1.2, col = 'black')
#legend('top', legend = c('Main effect', 'Interaction'),
#       fill = c('blue', 'lightgrey'), bty = 'n', text.col = 'black', cex = 1.3)


dev.off()
# Google Drive: clean up old version, and save new one:
# drive_trash(paste(gpth,pfile,sep=""))
# drive_upload(file,path=gpth)
