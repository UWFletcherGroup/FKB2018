# ---------------------------------------
# Disclaimer: this script is adapted from Doug McNeall's script famous_intro.R,
# for creating some descriptive plots of the ensemble of the FAMOUS climate model.
# (figure BL_obs_ensemble_mean_sd and figure frac_pairs)
# Source is located in the famous-git directory at github.com/dougmcneall/famous-git
#
# ---------------------------------------
# cgf: April 20, 2018
# this version works with the revised F1850 simulations run by Bakr
# that have fixed a bug with the x parameter.
# The previous F1850 simulations did not perturb sulfate.
#
# May 4, 2018: re-ran n=350 training cases on Graham. Use those data here.
# ---------------------------------------
#
# Nov 2018: new version of the pairs plot with shaded background colour based on correlation

# -------------------------------------------------------------
# 0. Packages, functions and data
# these are sourced from the F1850_common_May2018.R setup script
# -------------------------------------------------------------
source('F1850_common_May2018.R')


# -------------------------------------------------------------------
# 1. Pairs plot of input space and All OUTPUT data
# -------------------------------------------------------------------
# extract parameter and variable names
old.names=names(full_data)
new.names<-c("x1\n(IN)","x2\n(IN)","x3\n(IN)","x4\n(IN)","x5\n(IN)","x6\n(IN)","x7\n(IN)","x8\n(IN)","x9\n(IN)","AOD\n(OUT)","CLDL\n(OUT)","FNET\n(OUT)","LCWF\n(OUT)",
           "PRECT\n(OUT)","QRL\n(OUT)","SCWF\n(OUT)","SS\n(OUT)","CS\n(OUT)")
colnames(full_data) <- new.names

library(RColorBrewer)
# get array of colours for correlation matrix:
cols = brewer.pal(11, "RdBu")   # goes from red to white to blue
pal = colorRampPalette(cols)
cor_colors = data.frame(correlation = seq(-1,1,0.01),
                        correlation_color = pal(201)[1:201])  # assigns a color for each r correlation value
cor_colors$correlation_color = as.character(cor_colors$correlation_color)

mypanel <- function(x,y,...){
  # pairs plot panel function
  ll <- par("usr")
  r <- cor(x, y,method="spearman",use="complete.obs")
  test <- cor.test(x,y)
  bgcolor = cor_colors[2+(-r+1)*100,2]    # converts correlation into a specific color
  rect(ll[1], ll[3], ll[2], ll[4], col=bgcolor , border = NA)
  points(x, y, ... )
}

t.p <- function(x, y, labels, cex, font, ...){
  ll <- par("usr")
  text(x, y, labels, cex, font, col = 'black', ...)

}

# create pairwise scatter plot
pfile='F1850_INPUTs_OUTPUTs_pairs_Nov2018_newPairs.pdf'
file = paste(plotpath,pfile,sep="")

pdf(file = file, width = 9, height = 9)

par(fg = 'grey90')
ndt=length(full_data[,1])
plotdata=full_data.stand
labnames<-names(plotdata)
labnames[c(10,17,18)]<-c("AOD","SS",expression(paste(lambda)))
par(lab=c(2,2,4))
pairs(plotdata, gap = 0.5,
      lower.panel = mypanel,
      upper.panel = NULL,
      label.pos = 0.7,
      text.panel = t.p,
      col = c(rep('black', ndt), 'red'),
      cex = c(rep(0.1,ndt), 0.75),
      pch = c(rep(19, ndt),19),
      las = 2,cex.axis=1.2,labels=labnames
)
dev.off()

# Google Drive: clean up old version, and save new one:
# drive_trash(paste(gpth,pfile,sep=""))
# drive_upload(file,path=gpth)
