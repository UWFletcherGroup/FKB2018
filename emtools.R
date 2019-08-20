# emtools.R
# A collection of tools for working with emulators.
# Doug McNeall
# dougmcneall@gmail.com


oaat.design <- function(design, n, med = TRUE, hold = NULL){
  # function for creating one-at-a-time design matrix
  # INPUTS:
  # design .... original design (e.g. a latin hypercube or output from expand.grid)
  # n ......... number of design points in each dimension 
  # med ....... Use median rather than mean?
  # hold ...... Use supplied value to hold the non-changing points
  #  
  # OUTPUTS:
  # ........... (n x nd) rows, nd columns design matrix, sweeping through parameter space
  
  oamat <- NULL
  
  nd <- ncol(design)
  
  
  if(med){
    meandes <- apply(design,2,median)	
  }
  
  else{
    meandes <- apply(design,2,mean)
  }
  
  if(is.null(hold) == FALSE){
    
    meandes <- hold
  }
  
  
  mindes <- apply(design,2,min)
  maxdes <- apply(design,2,max)
  
  for (j in 1:nd){
    # base matrix of 'best' values to put the sweep values into
    basemat <- matrix(meandes, nrow = n, ncol = nd , byrow = TRUE)
    # use seq and length.out
    vec <- seq(from = mindes[j], to = maxdes[j], length.out = n)
    basemat[ ,j] <- vec
    oamat <- rbind(oamat,basemat)
    
  }
  
  oamat
  
}



normalize <- function(X, wrt = NULL){
  # Normalize a matrix to [0,1] on a per-column basis.
  # Normalize relative to matrix wrt if included.
  
  f <- function(X){
    (X-min(X))/(max(X)-min(X))
  }
  
  # test to see if we have a matrix, array or data frame
  if(length(dim(X))==2){
    out <- apply(X,2,f)
  }
  
  else{	
    out <- f(X)
  }
  
  if(is.null(wrt) == FALSE){
    # if argument wrt is given
    
    n <- nrow(X)
    mmins <- t(kronecker(apply(wrt,2,min),t(rep(1,n))))
    mmaxs <- t(kronecker(apply(wrt,2,max),t(rep(1,n))))
    
    out <- (X-mmins)/(mmaxs-mmins)
    
  }
  
  out
}



remap <- function(dat,xvec,yvec){
  # reshape climate data for plotting
  
  out <- t(matrix(dat, nrow = length(yvec),ncol = length(xvec)))
  out
  
}


samp.beta <- function(n,mins,maxes,wrt,shape1 = 2,shape2 = 2){
  # n is desired number or samples
  # mins, maxes are desired limits of samples in original space
  # wrt is : what should samples be normalized to (typically design)
  # shapes are shape parameters for the beta dist.
  # OUTPUTS
  # x.unn is unnormalized beta sample (defined by mins, maxes)
  # x.nwrt is normalized with respect to matrix wrt
  
  n.beta <- rbeta(n = length(mins)*n, shape1 = shape1,shape2 = shape2)
  dim(n.beta) <- c(length(mins),n)
  rownames(n.beta) <- names(mins)
  n.beta <- (t(n.beta))
  
  
  # beta distributions in the space defined by mins, maxes
  x.unn <- unnormalize(n.beta,mins,maxes)
  
  # Then, normalized to the space desribed by wrt
  x.nwrt <- n.wrt(x.unn,wrt)
  
  return(list(x.unn = x.unn,x.nwrt = x.nwrt))
  
}


samp.norm <- function(n, means, sds){
  # Sample from a normal distribution, and place in
  # a m = length(mins or maxes) x n matrix.
  
  out <- rnorm(n=length(means)*n, mean=means , sd = sds)
  dim(out) <- c(length(means),n)
  rownames(out) <- names(means)
  return(t(out))
}


samp.unif <- function(n, mins, maxes){
  # Sample from a uniform distribution, and place in
  # a m = length(mins or maxes) x n matrix.
  out <- runif(n=length(mins)*n, min=mins , max = maxes)
  dim(out) <- c(length(mins),n)
  rownames(out) <- names(mins)
  return(t(out))
}


taat.design <- function(X, n, means = NULL){
  # Build a two at a time emulator design
  # hold all of the other parameters at their mid values
  
  maxes <- apply(X, 2, max)
  mins  <- apply(X, 2, min)
  
  if(is.null(means)) means <- apply(X, 2, mean)
  
  nip <- ncol(X) # number of input parameters
  
  col.ix <- combn(1:nip,2)
  
  em.vec <- seq(from = 0, to = 1, length.out = n)
  
  des.cols <- expand.grid(em.vec, em.vec)
  
  holder <- matrix(means, ncol = nip, nrow = nrow(des.cols), byrow = TRUE)
  
  out <- NULL
  
  for(i in 1:ncol(col.ix)){
    
    mat.part <- holder
    
    colu <- col.ix[,i]
    
    mat.part[, as.matrix(colu[1])] <- des.cols[,1]
    mat.part[, as.matrix(colu[2])] <- des.cols[,2]
    
    out <- rbind(out, mat.part)
    
  }
  
  return(list(des = out, ix = col.ix))
  
}



unnormalize <- function(n,un.mins,un.maxes){
  # Return a normalized matrix to it's
  # un-normalized state, given a vector of
  # mins and a vector of maxes.
  
  un <- sweep(n,2,(un.maxes - un.mins),FUN = '*')
  unnormed <- sweep(un,2,un.mins,FUN ="+")
  unnormed
  
}