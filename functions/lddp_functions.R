require(MASS)

dpm <- function(y, prior, mcmc, standardise = FALSE) {
  
  set.seed(124)
  multinom <- function(probs) {
    probs <- t(apply(probs,1,cumsum)) 
    res <- rowSums(probs - runif(nrow(probs)) < 0) + 1 
    return(res)  
  }
  
  n <- length(y)
  yt <- y
  if(standardise) {
    yt <- (y-mean(y))/sd(y)
  }
  
  m0 <- prior$m0
  S0 <- prior$S0
  a <- prior$a
  b <- prior$b
  aalpha <- prior$aalpha
  balpha <- prior$balpha
  L <- prior$L
  
  nburn <- mcmc$nburn
  nsave <- mcmc$nsave
  nskip <- mcmc$nskip
  nsim <- nsave*nskip+nburn
  
  p <- ns <- rep(0,L)
  v <- rep(1/L,L)
  v[L] <- 1
  
  prop <- matrix(NA_real_, nrow = n, ncol = L)
  
  z <- matrix(NA_real_, nrow = nsim, ncol = n)
  z_tmp <- vector(length = n)
  
  z[1,] <- rep(1,n)
  
  P <- Mu <- Sigma2 <- matrix(NA_real_, nrow = nsim, ncol = L)
  
  Mu_tmp <- Sigma2_tmp  <- vector(length = L)
  
  Mu[1,] <- rep(mean(yt), L)
  Sigma2[1,] <- rep(var(yt), L)
  
  Mu_tmp <- Mu[1,]
  Sigma2_tmp <- Sigma2[1,]
  
  alpha <- numeric(nsim)
  alpha[1] <- 1
  
  for(i in 2:nsim) {
    Sigma2_tmp <- Sigma2[i-1,]
    
    cumv <- cumprod(1-v)
    p[1] <- v[1]
    p[2:L] <- v[2:L]*cumv[1:(L-1)]
    
    for(l in 1:L){
      prop[,l] <- p[l]*dnorm(yt, mean = Mu_tmp[l], sd = sqrt(Sigma2_tmp[l]))
    }
    prob <- prop/rowSums(prop)
    
    z_tmp <- multinom(prob)
    ns <- sapply(1:L, function(x, v) sum(v == x), v = z_tmp)
    yt_z_l <- sapply(1:L, function(x, v, y) sum(y[v == x]), v = z_tmp, y = yt)
    
    v[1:(L-1)] <- rbeta(L-1, 1 + ns[1:(L-1)], alpha[i-1] + rev(cumsum(rev(ns[-1]))))
    alpha[i] <- rgamma(1, shape = aalpha + L - 1, balpha - sum(log(1 - v[1:(L-1)])))
    
    varmu <- 1/((1/S0) + (ns/Sigma2_tmp))
    meanmu <- ((yt_z_l/Sigma2_tmp) + (m0/S0))/((1/S0) + (ns/Sigma2_tmp))
    Mu_tmp <- rnorm(L, mean = meanmu, sd = sqrt(varmu))
    
    yt_z_l_mu <- sapply(1:L, function(x, v, y, mu) sum((y[v == x] - mu[x])^2), v = z_tmp, y = yt, mu = Mu_tmp)
    
    Sigma2_tmp <- 1/rgamma(L, a + ns/2, b + 0.5*yt_z_l_mu)
    
    P[i,] <- p
    z[i,] <- z_tmp
    Mu[i,] <- Mu_tmp
    Sigma2[i,] <- Sigma2_tmp
  }
  if(standardise){
    Mu <- sd(y)*Mu + mean(y)
    Sigma2 <- var(y)*Sigma2
  }
  
  res <- list()
  res$z <- z[seq(nburn+1, nsim, by = nskip),]
  res$P <- P[seq(nburn+1, nsim, by = nskip),]
  res$Mu <- Mu[seq(nburn+1, nsim, by = nskip),]
  res$Sigma2 <- Sigma2[seq(nburn+1, nsim, by = nskip),]
  return(res)
}


dpm_survival <- function(y, survival, prior, mcmc, standardise = FALSE) {
  ##survival=0: living; survival=1: dead
  set.seed(124)
  multinom <- function(probs) {
    probs <- t(apply(probs,1,cumsum)) 
    res <- rowSums(probs - runif(nrow(probs)) < 0) + 1 
    return(res)  
  }
  
  n <- length(y)
  yt <- y
  if(standardise) {
    yt <- (y-mean(y))/sd(y)
  }
  yt_stand=yt
  
  
  m0 <- prior$m0
  S0 <- prior$S0
  a <- prior$a
  b <- prior$b
  aalpha <- prior$aalpha
  balpha <- prior$balpha
  L <- prior$L
  
  nburn <- mcmc$nburn
  nsave <- mcmc$nsave
  nskip <- mcmc$nskip
  nsim <- nsave*nskip+nburn
  
  p <- ns <- rep(0,L)
  v <- rep(1/L,L)
  v[L] <- 1
  
  prop <- matrix(NA_real_, nrow = n, ncol = L)
  
  z <- matrix(NA_real_, nrow = nsim, ncol = n)
  
  z_tmp <- rep(1,n)
  z[1,]=z_tmp
  
  P <- Mu <- Sigma2 <- matrix(NA_real_, nrow = nsim, ncol = L)
  
  Mu_tmp <- Sigma2_tmp  <- vector(length = L)
  
  Mu[1,] <- rep(mean(yt), L)
  Sigma2[1,] <- rep(var(yt), L)
  
  Mu_tmp <- Mu[1,]
  Sigma2_tmp <- Sigma2[1,]
  
  alpha <- numeric(nsim)
  alpha[1] <- 1
  
  y_impute_survival=function(yt_stand,survival,mu_all,sd_all,z){
    yt_temp=yt_stand
    n=length(yt_stand)
    
    missing_index=!survival
    z_mis=z[missing_index]
    
    region_upper=rep(Inf,sum(missing_index))
    region_lower=yt_stand[missing_index]
    #print(region_lower)
    #print(sum(region_lower>region_upper))
    yt_temp[missing_index]=TruncatedNormal::rtnorm(n=1, mu=mu_all[z_mis], sd=sd_all[z_mis],
                                                   lb=region_lower, ub=region_upper)
    return(yt_temp)
  }
  
  for(i in 2:nsim) {
    yt=y_impute_survival(yt_stand=yt_stand,survival=survival,
                         mu_all=Mu_tmp,sd_all=sqrt(Sigma2_tmp),z=z_tmp)
    
    cumv <- cumprod(1-v)
    p[1] <- v[1]
    p[2:L] <- v[2:L]*cumv[1:(L-1)]
    
    for(l in 1:L){
      prop[,l] <- p[l]*dnorm(yt, mean = Mu_tmp[l], sd = sqrt(Sigma2_tmp[l]))
    }
    prob <- prop/rowSums(prop)
    
    z_tmp <- multinom(prob)
    ns <- sapply(1:L, function(x, v) sum(v == x), v = z_tmp)
    yt_z_l <- sapply(1:L, function(x, v, y) sum(y[v == x]), v = z_tmp, y = yt)
    
    v[1:(L-1)] <- rbeta(L-1, 1 + ns[1:(L-1)], alpha[i-1] + rev(cumsum(rev(ns[-1]))))
    alpha[i] <- rgamma(1, shape = aalpha + L - 1, balpha - sum(log(1 - v[1:(L-1)])))
    
    varmu <- 1/((1/S0) + (ns/Sigma2_tmp))
    meanmu <- ((yt_z_l/Sigma2_tmp) + (m0/S0))/((1/S0) + (ns/Sigma2_tmp))
    Mu_tmp <- rnorm(L, mean = meanmu, sd = sqrt(varmu))
    
    yt_z_l_mu <- sapply(1:L, function(x, v, y, mu) sum((y[v == x] - mu[x])^2), v = z_tmp, y = yt, mu = Mu_tmp)
    
    Sigma2_tmp <- 1/rgamma(L, a + ns/2, b + 0.5*yt_z_l_mu)
    
    P[i,] <- p
    z[i,] <- z_tmp
    Mu[i,] <- Mu_tmp
    Sigma2[i,] <- Sigma2_tmp
  }
  if(standardise){
    Mu <- sd(y)*Mu + mean(y)
    Sigma2 <- var(y)*Sigma2
  }
  
  res <- list()
  res$z <- z[seq(nburn+1, nsim, by = nskip),]
  res$P <- P[seq(nburn+1, nsim, by = nskip),]
  res$Mu <- Mu[seq(nburn+1, nsim, by = nskip),]
  res$Sigma2 <- Sigma2[seq(nburn+1, nsim, by = nskip),]
  return(res)
}

ols.function <- function(X, y, vcov = FALSE) {
  res <- list()
  if(vcov) {
    res$vcov <- solve(crossprod(X)) 
    res$coeff <- res$vcov%*%crossprod(X,y)
  } else {
    res$coeff <- try(solve(crossprod(X), crossprod(X,y)), silent = TRUE) 
  }
  res 
}

lddp_moves <- function(y, X, prior, mcmc, standardise = TRUE) {
  multinom <- function(prob) {
    probs <- t(apply(prob,1,cumsum)) 
    res <- rowSums(probs - runif(nrow(probs)) < 0) + 1 
    return(res)  
  }
  
  swap <- function(vec, from, to) {
    tmp <- to[ match(vec, from) ]
    tmp[is.na(tmp)] <- vec[is.na(tmp)]
    return(tmp)
  }
  
  yt <- y
  if(standardise == TRUE) {
    yt <- (y-mean(y))/sd(y)
  }
  n <- length(y)
  k <- ncol(X)
  
  m <- prior$m0
  S <- prior$S0
  nu <- prior$nu
  psi <- prior$Psi
  a <- prior$a
  b <- prior$b
  aalpha <- prior$aalpha
  balpha <- prior$balpha
  L <- prior$L
  
  nburn <- mcmc$nburn
  nsave <- mcmc$nsave
  nskip <- mcmc$nskip
  nsim <- nburn + nsave*nskip
  
  p <- ns <- rep(0, L)
  v <- rep(1/L,L)
  v[L] <- 1
  
  z <- matrix(NA_real_, nrow = nsim, ncol = n, dimnames = list(1:nsim, 1:n))
  z_tmp <- vector(length = n)
  
  z[1,] <- rep(1,n)
  
  beta <- matrix(0, nrow = L, ncol = k)
  aux <- ols.function(X, yt)$coeff
  if(!inherits(aux, "try-error")) {
    for(l in 1:L) {
      beta[l,] <- aux
    }
  }
  
  tau <- rep(1/var(yt),L)
  prop <- prob <- matrix(NA_real_, nrow = n, ncol = L)
  
  P <- Tau <- matrix(NA_real_, nrow = nsim, ncol = L, dimnames = list(1:nsim, 1:L))
  Beta <- array(NA_real_,c(nsim,L,k), dimnames = list(1:nsim, 1:L, colnames(X)))
  alpha <- numeric(nsim)
  
  Beta[1,,] <- beta
  Tau[1,] <- tau
  alpha[1] <- 1
  
  mu <- mvrnorm(1, mu = m, Sigma = S)
  Sigmainv <- rWishart(1, df = nu, solve(nu*psi))[,,1] 
  
  for(i in 2:nsim) {
    cumv <- cumprod(1-v)
    p[1] <- v[1]
    p[2:L] <- v[2:L]*cumv[1:(L-1)]
    
    for(l in 1:L) {
      prop[,l] <- p[l]*dnorm(yt, mean = X%*%beta[l,], sd = sqrt(1/tau[l]))
    }
    
    prob <- prop/rowSums(prop)
    
    z_tmp <- multinom(prob)
    
    ns <- sapply(1:L, function(x, v) sum(v == x), v = z_tmp)
    
    v[1:(L-1)] <- rbeta(L-1, 1 + ns[1:(L-1)], alpha[i-1] + rev(cumsum(rev(ns[-1]))))
    alpha[i] <- rgamma(1, shape = aalpha + L - 1, balpha - sum(log(1 - v[1:(L-1)])))
    
    Sigmainv_mu <- Sigmainv%*%mu
    
    for(l in 1:L) {
      X_l  <- matrix(X[z_tmp == l, ], ncol = k, nrow = ns[l])
      yt_l <- yt[z_tmp == l]
      V <- solve(Sigmainv + tau[l]*crossprod(X_l))
      mu1 <- V%*%(Sigmainv_mu + tau[l]*crossprod(X_l, yt_l))
      beta[l,] <- mvrnorm(1, mu = mu1, Sigma = V)
      
      aux <- yt_l - X_l%*%beta[l,]
      tau[l] <- rgamma(1, shape = a + (ns[l]/2), rate = b + 0.5*crossprod(aux))
    }
    
    #Label-switching moves
    if(L > 1) {
      # Move 1
      act_comp <- (1:L)[ns != 0]
      if(length(act_comp) > 1) {
        comp_sel <- sample(act_comp, 2, replace = FALSE)
        prob_change <- min(1, (p[comp_sel[1]]/p[comp_sel[2]])^(ns[comp_sel[2]]-ns[comp_sel[1]]))
        if(!is.na(prob_change) && runif(1) <= prob_change) {
          beta[comp_sel,]  <- beta[rev(comp_sel),]
          tau[comp_sel]  <- tau[rev(comp_sel)]
          p[comp_sel] <- p[rev(comp_sel)]	        	
          z_tmp <- swap(z_tmp, comp_sel, rev(comp_sel))
        }
      }
    
      # Move 2
      comp_sel <- sample(2:L, 1, replace = FALSE)
      prob_change <- min(1, ((1-v[comp_sel])^ns[comp_sel - 1])/((1-v[comp_sel - 1])^ns[comp_sel]))
      if(!is.na(prob_change) && runif(1) <= prob_change) {	        	
        comp_sel <- c(comp_sel, comp_sel - 1)
        beta[comp_sel,]  <- beta[rev(comp_sel),]
        tau[comp_sel]  <- tau[rev(comp_sel)]
        p[comp_sel] <- p[rev(comp_sel)]
        v[comp_sel] <- v[rev(comp_sel)]
        z_tmp <- swap(z_tmp, comp_sel, rev(comp_sel))
      }
    }  
      
    S_inv <- solve(S)
    Vaux <- solve(S_inv + L*Sigmainv)
    if(k == 1) {
      meanmu <- Vaux%*%(S_inv%*%m + Sigmainv%*%sum(beta))
    } else {
      meanmu <- Vaux%*%(S_inv%*%m + Sigmainv%*%colSums(beta))
    }
    mu <- mvrnorm(1, mu = meanmu, Sigma = Vaux)
    
    Vaux1 <- 0
    for(l in 1:L) {
      Vaux1 <- Vaux1 + tcrossprod(beta[l,] - mu)
    }
    Sigmainv <- rWishart(1, nu + L, solve(nu*psi + Vaux1))[,,1]
    
    P[i,] <- p
    z[i,] <- z_tmp
    Beta[i,,] <- beta
    Tau[i,] <- tau
  }
  
  if (standardise == TRUE) {
    Beta[,,1] <- sd(y)*Beta[,,1] + mean(y)
    if(k > 1) {
      Beta[,,2:k] <- sd(y)*Beta[,,2:k]
    }
    Sigma2 <- var(y)*(1/Tau)
  } else {
    Sigma2 <- 1/Tau
  }
  
  res <- list()
  #latent component indicator
  res$z <- z[seq(nburn+1, nsim, by = nskip),]
  #weights
  res$P <- P[seq(nburn+1, nsim, by = nskip),]
  #regression coefficients
  res$Beta <- Beta[seq(nburn+1, nsim, by = nskip),,,drop = FALSE]
  #Normal components variances
  res$Sigma2 <- Sigma2[seq(nburn+1, nsim, by = nskip),]
  res
}

lddp <- function(y, X, prior, mcmc, standardise = TRUE) {
  multinom <- function(prob) {
    probs <- t(apply(prob,1,cumsum)) 
    res <- rowSums(probs - runif(nrow(probs)) < 0) + 1 
    return(res)  
  }
  
  yt <- y
  if(standardise == TRUE) {
    yt <- (y-mean(y))/sd(y)
  }
  n <- length(y)
  k <- ncol(X)
  
  m <- prior$m0
  S <- prior$S0
  nu <- prior$nu
  psi <- prior$Psi
  a <- prior$a
  b <- prior$b
  aalpha <- prior$aalpha
  balpha <- prior$balpha
  L <- prior$L
  
  nburn <- mcmc$nburn
  nsave <- mcmc$nsave
  nskip <- mcmc$nskip
  nsim <- nburn + nsave*nskip
  
  p <- ns <- rep(0, L)
  v <- rep(1/L,L)
  v[L] <- 1
  
  z <- matrix(NA_real_, nrow = nsim, ncol = n, dimnames = list(1:nsim, 1:n))
  z_tmp <- vector(length = n)
  
  z[1,] <- rep(1,n)
  
  beta <- matrix(0, nrow = L, ncol = k)
  aux <- ols.function(X, yt)$coeff
  if(!inherits(aux, "try-error")) {
    for(l in 1:L) {
      beta[l,] <- aux
    }
  }
  
  tau <- rep(1/var(yt),L)
  prop <- prob <- matrix(NA_real_, nrow = n, ncol = L)
  
  P <- Tau <- matrix(NA_real_, nrow = nsim, ncol = L, dimnames = list(1:nsim, 1:L))
  Beta <- array(NA_real_,c(nsim,L,k), dimnames = list(1:nsim, 1:L, colnames(X)))
  alpha <- numeric(nsim)
  
  Beta[1,,] <- beta
  Tau[1,] <- tau
  alpha[1] <- 1
  
  mu <- mvrnorm(1, mu = m, Sigma = S)
  Sigmainv <- rWishart(1, df = nu, solve(nu*psi))[,,1] 
  
  for(i in 2:nsim) {
    cumv <- cumprod(1-v)
    p[1] <- v[1]
    p[2:L] <- v[2:L]*cumv[1:(L-1)]
    
    for(l in 1:L) {
      prop[,l] <- p[l]*dnorm(yt, mean = X%*%beta[l,], sd = sqrt(1/tau[l]))
    }
    
    prob <- prop/rowSums(prop)
    
    z_tmp <- multinom(prob)
    
    ns <- sapply(1:L, function(x, v) sum(v == x), v = z_tmp)
    
    v[1:(L-1)] <- rbeta(L-1, 1 + ns[1:(L-1)], alpha[i-1] + rev(cumsum(rev(ns[-1]))))
    alpha[i] <- rgamma(1, shape = aalpha + L - 1, balpha - sum(log(1 - v[1:(L-1)])))
    
    Sigmainv_mu <- Sigmainv%*%mu
    
    for(l in 1:L) {
      X_l  <- matrix(X[z_tmp == l, ], ncol = k, nrow = ns[l])
      yt_l <- yt[z_tmp == l]
      V <- solve(Sigmainv + tau[l]*crossprod(X_l))
      mu1 <- V%*%(Sigmainv_mu + tau[l]*crossprod(X_l, yt_l))
      beta[l,] <- mvrnorm(1, mu = mu1, Sigma = V)
      
      aux <- yt_l - X_l%*%beta[l,]
      tau[l] <- rgamma(1, shape = a + (ns[l]/2), rate = b + 0.5*crossprod(aux))
    }
    
    S_inv <- solve(S)
    Vaux <- solve(S_inv + L*Sigmainv)
    if(k == 1) {
      meanmu <- Vaux%*%(S_inv%*%m + Sigmainv%*%sum(beta))
    } else {
      meanmu <- Vaux%*%(S_inv%*%m + Sigmainv%*%colSums(beta))
    }
    mu <- mvrnorm(1, mu = meanmu, Sigma = Vaux)
    
    Vaux1 <- 0
    for(l in 1:L) {
      Vaux1 <- Vaux1 + tcrossprod(beta[l,] - mu)
    }
    
    Sigmainv <- rWishart(1, nu + L, solve(nu*psi + Vaux1))[,,1]
    
    P[i,] <- p
    z[i,] <- z_tmp
    Beta[i,,] <- beta
    Tau[i,] <- tau
  }
  
  if (standardise == TRUE) {
    Beta[,,1] <- sd(y)*Beta[,,1] + mean(y)
    if(k > 1) {
      Beta[,,2:k] <- sd(y)*Beta[,,2:k]
    }
    Sigma2 <- var(y)*(1/Tau)
  } else {
    Sigma2 <- 1/Tau
  }
  
  res <- list()
  #latent component indicator
  res$z <- z[seq(nburn+1, nsim, by = nskip),]
  #weights
  res$P <- P[seq(nburn+1, nsim, by = nskip),]
  #regression coefficients
  res$Beta <- Beta[seq(nburn+1, nsim, by = nskip),,,drop = FALSE]
  #Normal components variances
  res$Sigma2 <- Sigma2[seq(nburn+1, nsim, by = nskip),]
  res
}


lddp_moves_inte <- function(y, X, prior, mcmc, k_min=0,k_max=30,standardise = FALSE) {
  # y = con_mmse2_com$m12
  # X = X
  # prior = prior
  # mcmc = mcmc
  # standardise = FALSE
  multinom <- function(prob) {
    probs <- t(apply(prob,1,cumsum)) 
    res <- rowSums(probs - runif(nrow(probs)) < 0) + 1 
    return(res)  
  }
  
  swap <- function(vec, from, to) {
    tmp <- to[ match(vec, from) ]
    tmp[is.na(tmp)] <- vec[is.na(tmp)]
    return(tmp)
  }
  
  yt <- y
  if(standardise == TRUE) {
    yt <- (y-mean(y))/sd(y)
  }
  
  yt_stand=yt
  n <- length(y)
  k <- ncol(X)
  
  m <- prior$m0
  S <- prior$S0
  nu <- prior$nu
  psi <- prior$Psi
  a <- prior$a
  b <- prior$b
  aalpha <- prior$aalpha
  balpha <- prior$balpha
  L <- prior$L
  
  nburn <- mcmc$nburn
  nsave <- mcmc$nsave
  nskip <- mcmc$nskip
  nsim <- nburn + nsave*nskip
  
  p <- ns <- rep(0, L)
  v <- rep(1/L,L)
  v[L] <- 1
  
  beta <- matrix(0, nrow = L, ncol = k)
  aux <- ols.function(X, yt)$coeff
  if(!inherits(aux, "try-error")) {
    for(l in 1:L) {
      beta[l,] <- aux
    }
  }
  
  tau <- rep(1/var(yt),L)
  prop <- prob <- matrix(NA_real_, nrow = n, ncol = L)
  
  P <- Tau <- matrix(NA_real_, nrow = nsim, ncol = L, dimnames = list(1:nsim, 1:L))
  Beta <- array(NA_real_,c(nsim,L,k), dimnames = list(1:nsim, 1:L, colnames(X)))
  alpha <- numeric(nsim)
  
  Beta[1,,] <- beta
  Tau[1,] <- tau
  alpha[1] <- 1
  
  y_post=matrix(0,nrow=nsim,ncol = length(yt))
  
  mu <- mvrnorm(1, mu = m, Sigma = S)
  Sigmainv <- rWishart(1, df = nu, solve(nu*psi))[,,1] 
  
  # cumv <- cumprod(1-v)
  # p[1] <- v[1]
  # p[2:L] <- v[2:L]*cumv[1:(L-1)]
  # 
  # P[1,]=p
  
  # z <- matrix(NA_real_, nrow = nsim, ncol = n, dimnames = list(1:nsim, 1:n))
  # 
  # for(l in 1:L) {
  #   prop[,l] <- p[l]*dnorm(yt, mean = X%*%beta[l,], sd = sqrt(1/tau[l]))
  # }
  # 
  # prob <- prop/rowSums(prop)
  # 
  # z_tmp <- multinom(prob)
  # z[1,] <- z_tmp
  
  z <- matrix(NA_real_, nrow = nsim, ncol = n, dimnames = list(1:nsim, 1:n))
  z_tmp <- rep(1,n)
  
  z[1,] <- z_tmp
  
  y_impute=function(yt_stand,beta,tau,z_tmp,k_min,k_max){
    n=length(yt_stand)
    for(j in 1:n){#j=201
      region_upper=yt_stand[j]+1
      region_lower=yt_stand[j]
      if(region_lower<k_min|is.na(region_lower)){
        region_lower=-Inf
      }
      
      if(region_upper>k_max|is.na(region_upper)){
        region_upper=Inf
      }
      
      mu_bar=sum(X[j,]*beta[z_tmp[j],]);sigma_bar=1/tau[z_tmp[j]]
      
      yt_stand[j]=TruncatedNormal::rtnorm(n=1, mu=mu_bar, sd=sqrt(sigma_bar),
                                          lb=region_lower, ub=region_upper)
      
      
      #print(j)
      
    }
    return(yt_stand)
  }
  
  y_impute2=function(yt_stand,beta,tau,z_tmp,k_min,k_max){
    n=length(yt_stand)
    
    region_upper=yt_stand+1
    region_lower=yt_stand
    
    region_upper[region_upper>k_max|is.na(region_upper)]=Inf
    region_lower[region_lower<k_min|is.na(region_lower)]=-Inf
    
    mu_bar=diag(X %*% t(beta[z_tmp,]));sigma_bar=1/tau[z_tmp]
    
    yt_new=TruncatedNormal::rtnorm(n=1, mu=mu_bar, sd=sqrt(sigma_bar),
                                   lb=region_lower, ub=region_upper)
    return(yt_new)
  }
  
  for(i in 2:nsim) {
    #i=2
    cumv <- cumprod(1-v)
    p[1] <- v[1]
    p[2:L] <- v[2:L]*cumv[1:(L-1)]
    
    yt=y_impute2(yt_stand=yt_stand,
                 beta=beta,tau=tau,z_tmp=z_tmp,k_min=k_min,k_max=k_max)
    
    for(l in 1:L) {
      prop[,l] <- p[l]*dnorm(yt, mean = X%*%beta[l,], sd = sqrt(1/tau[l]))
    }
    
    prob <- prop/rowSums(prop)
    
    z_tmp <- multinom(prob)
    
    ns <- sapply(1:L, function(x, v) sum(v == x), v = z_tmp)
    
    v[1:(L-1)] <- rbeta(L-1, 1 + ns[1:(L-1)], alpha[i-1] + rev(cumsum(rev(ns[-1]))))
    alpha[i] <- rgamma(1, shape = aalpha + L - 1, balpha - sum(log(1 - v[1:(L-1)])))
    
    Sigmainv_mu <- Sigmainv%*%mu
    
    for(l in 1:L) {
      X_l  <- matrix(X[z_tmp == l, ], ncol = k, nrow = ns[l])
      yt_l <- yt[z_tmp == l]
      V <- solve(Sigmainv + tau[l]*crossprod(X_l))
      mu1 <- V%*%(Sigmainv_mu + tau[l]*crossprod(X_l, yt_l))
      beta[l,] <- mvrnorm(1, mu = mu1, Sigma = V)
      
      aux <- yt_l - X_l%*%beta[l,]
      tau[l] <- rgamma(1, shape = a + (ns[l]/2), rate = b + 0.5*crossprod(aux))
    }
    
    #Label-switching moves
    if(L > 1) {
      # Move 1
      act_comp <- (1:L)[ns != 0]
      if(length(act_comp) > 1) {
        comp_sel <- sample(act_comp, 2, replace = FALSE)
        prob_change <- min(1, (p[comp_sel[1]]/p[comp_sel[2]])^(ns[comp_sel[2]]-ns[comp_sel[1]]))
        if(!is.na(prob_change) && runif(1) <= prob_change) {
          beta[comp_sel,]  <- beta[rev(comp_sel),]
          tau[comp_sel]  <- tau[rev(comp_sel)]
          p[comp_sel] <- p[rev(comp_sel)]	        	
          z_tmp <- swap(z_tmp, comp_sel, rev(comp_sel))
        }
      }
      
      # Move 2
      comp_sel <- sample(2:L, 1, replace = FALSE)
      prob_change <- min(1, ((1-v[comp_sel])^ns[comp_sel - 1])/((1-v[comp_sel - 1])^ns[comp_sel]))
      if(!is.na(prob_change) && runif(1) <= prob_change) {	        	
        comp_sel <- c(comp_sel, comp_sel - 1)
        beta[comp_sel,]  <- beta[rev(comp_sel),]
        tau[comp_sel]  <- tau[rev(comp_sel)]
        p[comp_sel] <- p[rev(comp_sel)]
        v[comp_sel] <- v[rev(comp_sel)]
        z_tmp <- swap(z_tmp, comp_sel, rev(comp_sel))
      }
    }  
    
    S_inv <- solve(S)
    Vaux <- solve(S_inv + L*Sigmainv)
    if(k == 1) {
      meanmu <- Vaux%*%(S_inv%*%m + Sigmainv%*%sum(beta))
    } else {
      meanmu <- Vaux%*%(S_inv%*%m + Sigmainv%*%colSums(beta))
    }
    mu <- mvrnorm(1, mu = meanmu, Sigma = Vaux)
    
    Vaux1 <- 0
    for(l in 1:L) {
      Vaux1 <- Vaux1 + tcrossprod(beta[l,] - mu)
    }
    Sigmainv <- rWishart(1, nu + L, solve(nu*psi + Vaux1))[,,1]
    
    P[i,] <- p
    z[i,] <- z_tmp
    Beta[i,,] <- beta
    Tau[i,] <- tau
    y_post[i,]<-yt
    print(i)
  }
  
  if (standardise == TRUE) {
    Beta[,,1] <- sd(y)*Beta[,,1] + mean(y)
    if(k > 1) {
      Beta[,,2:k] <- sd(y)*Beta[,,2:k]
    }
    Sigma2 <- var(y)*(1/Tau)
  } else {
    Sigma2 <- 1/Tau
  }
  
  res <- list()
  #latent component indicator
  res$z <- z[seq(nburn+1, nsim, by = nskip),]
  #weights
  res$P <- P[seq(nburn+1, nsim, by = nskip),]
  #regression coefficients
  res$Beta <- Beta[seq(nburn+1, nsim, by = nskip),,,drop = FALSE]
  #Normal components variances
  res$Sigma2 <- Sigma2[seq(nburn+1, nsim, by = nskip),]
  
  res$y_post <- y_post[seq(nburn+1, nsim, by = nskip),]
  res
}

lddp_inte <- function(y, X, prior, mcmc, k_min=0,k_max=30,standardise = TRUE) {
  multinom <- function(prob) {
    probs <- t(apply(prob,1,cumsum)) 
    res <- rowSums(probs - runif(nrow(probs)) < 0) + 1 
    return(res)  
  }
  
  yt <- y
  if(standardise == TRUE) {
    yt <- (y-mean(y))/sd(y)
  }
  yt_stand=yt
  n <- length(y)
  k <- ncol(X)
  
  m <- prior$m0
  S <- prior$S0
  nu <- prior$nu
  psi <- prior$Psi
  a <- prior$a
  b <- prior$b
  aalpha <- prior$aalpha
  balpha <- prior$balpha
  L <- prior$L
  
  nburn <- mcmc$nburn
  nsave <- mcmc$nsave
  nskip <- mcmc$nskip
  nsim <- nburn + nsave*nskip
  
  p <- ns <- rep(0, L)
  v <- rep(1/L,L)
  v[L] <- 1
  
  z <- matrix(NA_real_, nrow = nsim, ncol = n, dimnames = list(1:nsim, 1:n))
  z_tmp <- rep(1,n)
  
  z[1,] <- rep(1,n)
  
  beta <- matrix(0, nrow = L, ncol = k)
  aux <- ols.function(X, yt)$coeff
  if(!inherits(aux, "try-error")) {
    for(l in 1:L) {
      beta[l,] <- aux
    }
  }
  
  tau <- rep(1/var(yt),L)
  prop <- prob <- matrix(NA_real_, nrow = n, ncol = L)
  
  P <- Tau <- matrix(NA_real_, nrow = nsim, ncol = L, dimnames = list(1:nsim, 1:L))
  Beta <- array(NA_real_,c(nsim,L,k), dimnames = list(1:nsim, 1:L, colnames(X)))
  alpha <- numeric(nsim)
  
  Beta[1,,] <- beta
  Tau[1,] <- tau
  alpha[1] <- 1
  
  mu <- mvrnorm(1, mu = m, Sigma = S)
  Sigmainv <- rWishart(1, df = nu, solve(nu*psi))[,,1] 
  
  y_post=matrix(0,nrow=nsim,ncol = length(yt))
  
  y_impute=function(yt_stand,beta,tau,z_tmp,k_min,k_max){
    n=length(yt_stand)
    for(j in 1:n){#j=201
      region_upper=yt_stand[j]+1
      region_lower=yt_stand[j]
      if(region_lower<k_min|is.na(region_lower)){
        region_lower=-Inf
      }
      
      if(region_upper>k_max|is.na(region_upper)){
        region_upper=Inf
      }
      
      mu_bar=sum(X[j,]*beta[z_tmp[j],]);sigma_bar=1/tau[z_tmp[j]]
      
      yt_stand[j]=TruncatedNormal::rtnorm(n=1, mu=mu_bar, sd=sqrt(sigma_bar),
                                          lb=region_lower, ub=region_upper)
      
      
      #print(j)
      
    }
    return(yt_stand)
  }
  
  y_impute2=function(yt_stand,beta,tau,z_tmp,k_min,k_max){
    n=length(yt_stand)
    
    region_upper=yt_stand+1
    region_lower=yt_stand
    
    region_upper[region_upper>k_max|is.na(region_upper)]=Inf
    region_lower[region_lower<k_min|is.na(region_lower)]=-Inf
    
    mu_bar=diag(X %*% t(beta[z_tmp,]));sigma_bar=1/tau[z_tmp]
    
    yt_new=TruncatedNormal::rtnorm(n=1, mu=mu_bar, sd=sqrt(sigma_bar),
                                   lb=region_lower, ub=region_upper)
    return(yt_new)
  }
  
  for(i in 2:nsim) {
    cumv <- cumprod(1-v)
    p[1] <- v[1]
    p[2:L] <- v[2:L]*cumv[1:(L-1)]
    
    yt=y_impute2(yt_stand=yt_stand,
                 beta=beta,tau=tau,z_tmp=z_tmp,k_min=k_min,k_max=k_max)
    
    for(l in 1:L) {
      prop[,l] <- p[l]*dnorm(yt, mean = X%*%beta[l,], sd = sqrt(1/tau[l]))
    }
    
    prob <- prop/rowSums(prop)
    
    z_tmp <- multinom(prob)
    
    ns <- sapply(1:L, function(x, v) sum(v == x), v = z_tmp)
    
    v[1:(L-1)] <- rbeta(L-1, 1 + ns[1:(L-1)], alpha[i-1] + rev(cumsum(rev(ns[-1]))))
    alpha[i] <- rgamma(1, shape = aalpha + L - 1, balpha - sum(log(1 - v[1:(L-1)])))
    
    Sigmainv_mu <- Sigmainv%*%mu
    
    for(l in 1:L) {
      X_l  <- matrix(X[z_tmp == l, ], ncol = k, nrow = ns[l])
      yt_l <- yt[z_tmp == l]
      V <- solve(Sigmainv + tau[l]*crossprod(X_l))
      mu1 <- V%*%(Sigmainv_mu + tau[l]*crossprod(X_l, yt_l))
      beta[l,] <- mvrnorm(1, mu = mu1, Sigma = V)
      
      aux <- yt_l - X_l%*%beta[l,]
      tau[l] <- rgamma(1, shape = a + (ns[l]/2), rate = b + 0.5*crossprod(aux))
    }
    
    S_inv <- solve(S)
    Vaux <- solve(S_inv + L*Sigmainv)
    if(k == 1) {
      meanmu <- Vaux%*%(S_inv%*%m + Sigmainv%*%sum(beta))
    } else {
      meanmu <- Vaux%*%(S_inv%*%m + Sigmainv%*%colSums(beta))
    }
    mu <- mvrnorm(1, mu = meanmu, Sigma = Vaux)
    
    Vaux1 <- 0
    for(l in 1:L) {
      Vaux1 <- Vaux1 + tcrossprod(beta[l,] - mu)
    }
    
    Sigmainv <- rWishart(1, nu + L, solve(nu*psi + Vaux1))[,,1]
    
    P[i,] <- p
    z[i,] <- z_tmp
    Beta[i,,] <- beta
    Tau[i,] <- tau
    y_post[i,]<-yt
  }
  
  if (standardise == TRUE) {
    Beta[,,1] <- sd(y)*Beta[,,1] + mean(y)
    if(k > 1) {
      Beta[,,2:k] <- sd(y)*Beta[,,2:k]
    }
    Sigma2 <- var(y)*(1/Tau)
  } else {
    Sigma2 <- 1/Tau
  }
  
  res <- list()
  #latent component indicator
  res$z <- z[seq(nburn+1, nsim, by = nskip),]
  #weights
  res$P <- P[seq(nburn+1, nsim, by = nskip),]
  #regression coefficients
  res$Beta <- Beta[seq(nburn+1, nsim, by = nskip),,,drop = FALSE]
  #Normal components variances
  res$Sigma2 <- Sigma2[seq(nburn+1, nsim, by = nskip),]
  res$y_post=y_post[seq(nburn+1, nsim, by = nskip),]
  
  res
}

predic_lddp<-function(fit_lddp,X){
  n=dim(X)[1]
  nsim=dim(fit_lddp$Beta)[1]
  L=dim(fit_lddp$Beta)[2]
  means= array(NA_real_,c(nsim,L,n))
  sigma2=fit_lddp$Sigma2
  beta=fit_lddp$Beta
  y_post=matrix(0,nrow = n,ncol = nsim)
  
  for(i in 1:nsim){
    for (l in 1:L){
      means[i,l,]=X%*%beta[i,l,]
    }
    
  }
  
  for(i in 1:nsim){
    for (j in 1:n){
      y_post[j,i]=mixAK::rMVNmixture(n=1,weight = fit_lddp$P[i,],mean = means[i,,j],Sigma = sigma2[i,])
    }
    
  }
  return(y_post)
}

density_lddp<-function(fit_lddp,X,grid){
  n=dim(X)[1]
  nsim=dim(fit_lddp$Beta)[1]
  L=dim(fit_lddp$Beta)[2]
  means= array(NA_real_,c(nsim,L,n))
  sigma2=fit_lddp$Sigma2
  beta=fit_lddp$Beta
  y_dens=array(NA_real_,c(nsim,length(grid),n))
  
  for(i in 1:nsim){
    for (l in 1:L){
      means[i,l,]=X%*%beta[i,l,]
    }
    
  }
  
  for(i in 1:nsim){
    for (j in 1:n){
      y_dens[i,,j]=mixAK::dMVNmixture(x=grid,weight = fit_lddp$P[i,],mean = means[i,,j],Sigma = sigma2[i,])
    }
    
  }
  return(y_dens)
}
