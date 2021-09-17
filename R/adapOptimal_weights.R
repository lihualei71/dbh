library(locfdr)
library(EbayesThresh)
adaptive_optimal.weights <- function(groups, 
                                     zvals, 
                                     alpha, 
                                     side, 
                                     type = "max", 
                                     pi0Est = F) {
  weights <- c()
  lfdr_g <- list()
  ngroups <- length(unique(groups))
  pi0_g <- c()
  n_g <- c()
  for(i in 1: ngroups) {
    g <- unique(groups)[i]
    zvals_g <- zvals[groups == g]
    n_g[i] <- sum(groups == g)
    lfdrres <- locfdr(zvals_g, nulltype = 0, plot = 0)
    if(side == "one") {
      lfdr_g[[i]] <- lfdrres$fdr
      lfdr_g[[i]][which(sign(zvals_g) < 0)] <- 1
    } else {
      #lfdr_g[[i]] <- lfdrres$fdr
    }
    pi0_g[i] <- ifelse(pi0Est, lfdrres$fp0[1,3], 1)
  }
  n <- sum(n_g)
  pi_g <- n_g/n
  
  lfdr <- unlist(lfdr_g, F, F)
  lfdr.order <- order(lfdr)
  st.lfdr <- lfdr[lfdr.order]
  k = max(which(cumsum(st.lfdr)/(1:n) <= alpha))
  thr <- st.lfdr[k]
  
  st.groups <- groups[lfdr.order]
  st.pvals <- (zvals_pvals(zvals, side))[lfdr.order]
  
  t_g <- c()
  thr_diff <- c()
  for(j in 1:ngroups) {
    g <- unique(groups)[j]
    st.lfdr_g <- st.lfdr[st.groups == j]
    st.pvals_g <- st.pvals[st.groups == j]
    if(type == "lin") {
      t_g[j] <- approx(x = st.lfdr_g, y = st.pvals_g, xout = thr, yleft = 0, yright = 1, method="linear", ties = mean)$y
    } else if(type == "lin_exp") {
      t_g[j] <- exp(approx(x = st.lfdr_g, y = log(st.pvals_g), xout = thr, yleft = -99999, yright = 0, method="linear", ties = mean)$y)
    } else if(type == "max") {
      if(!any(which(st.lfdr_g < thr))) {
        t_g[j] <- 0
      } else {
        cut_index <- max(which(st.lfdr_g < thr))
        thr_diff[j] <- thr - st.lfdr_g[cut_index]
        t_g[j] <- st.pvals_g[cut_index]
      }
    } else {
      stop("Specify your type of estimation of t_g.")
    }
  }
  
  if(sum(t_g*pi_g*pi0_g) < 1e-6) {
    return(rep(1/sum(pi_g*pi0_g), length(zvals)))
  }
  #?
  t_g.init <- t_g
  
  t_g <- t_g/sum(t_g*pi_g*pi0_g)
  
  for(j in 1: ngroups) {
    g <- unique(groups)[j]
    weights[groups == g] <- t_g[j]
  }
  
  return(list(weights = weights, thr = thr, t_g = t_g, t_g.init = t_g.init, thr_diff = thr_diff))
}


oracle.weights <- function(alpha, n_g, pi0_g, mu1_g, pi0Est = T) {
  pi_g <- n_g/sum(n_g)
  #lambda_alpha <- solve.mFDR(alpha, n_g, pi0_g, mu1_g)
  mFDR1 <- function(c) {
    V <- c()
    R <- c()
    for(i in 1:length(pi0_g)) {
      pi0 = pi0_g[i]
      mu1 = mu1_g[i]
      n = n_g[i]
      t = lfdr.inverse(c, pi0, mu1)
      V[i] = n*pi0*t
      R[i] = V[i] + n*(1-pi0)*nonnull.cdf(t, mu1)
    }
    return(sum(V)/sum(R))
  }
  lambda_alpha <- vecbinsolv(zf = alpha, fun = mFDR1 , tlo = 0, thi = 1, nits = 30)
  w <- c()
  t <- c()
  for(i in 1:length(n_g)) {
    t[i] <- lfdr.inverse(lambda_alpha, pi0 = pi0_g[i], mu1 = mu1_g[i])
  }
  if(pi0Est) {
    t <- t/sum(t*pi0_g*pi_g)
  } else{
    t <- t/mean(t)
  }
  
  for(i in 1:length(n_g)) {
    w <- c(w, rep(t[i], n_g[i]))
  }
  return(w)
}


unique(oracle.weights(alpha, n_g, pi0_g, mu1_g, pi0Est = F))


lfdr.inverse <- function(lfdr, pi0, mu1){
  1 - pnorm(1/mu1*log((pi0/lfdr - pi0)/(1-pi0)) + mu1/2)
}

mFDR <- function(c, n_g, pi0_g, mu1_g) {
  V <- c()
  R <- c()
  for(i in 1:length(pi0_g)) {
    pi0 = pi0_g[i]
    mu1 = mu1_g[i]
    n = n_g[i]
    t = lfdr.inverse(c, pi0, mu1)
    V[i] = n*pi0*t
    R[i] = V[i] + n*(1-pi0)*nonnull.cdf(t, mu1)
  }
  return(sum(V)/sum(R))
}

pval.ds <- function(p, mu1, pi0){
  pi0+(1-pi0)*(dnorm(qnorm(p, lower.tail = F)-mu1)/dnorm(qnorm(p, lower.tail = F)))
}

nonnull.cdf <- function(p, mu1) {
  pnorm(qnorm(1-p) - mu1,lower.tail = F)
}