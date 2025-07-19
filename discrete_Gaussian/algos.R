require(pracma)
require(Rcpp)
require(RcppArmadillo)
sourceCpp('distou_revise.cpp')

#########################################################
## Codes for "Discrete Hamiltonian-Assisted Metropolis Sampling" for discrete Gaussian
#########################################################

#--- Implements the following MCMC algorithms
#--- Metropolis, GWG, NCG, AVG, V-DHAMS and O-DHAMS. 
#--- gradient: function for querying gradient of f(s)
#--- energy: function for querying f(s)
#--- 'distou_revise.cpp': cpp file for implementing the over-relaxation algorithm
#--- In all samplers below:
#         n         : Number of draws
#         seed      : Random seed for reproducibility
#         param_list: List of sampler-specific parameters

rho = 0.9
sigma= 5
nsize = 10
ndim = 8
grid = -nsize:nsize 
cov_mat = sigma^2*(rho*matrix(rep(1, ndim^2), ndim, ndim)+(1-rho)*diag(rep(1, ndim)))
inv_cov_mat = solve(cov_mat)


gradient = function(x){
  return(-inv_cov_mat%*%x)
}

energy = function(x) {
  return(-0.5*t(x)%*%inv_cov_mat%*%x)
}


#-- Metropolis with a neighborhood
#-- param_list: a list with parameter
#    r, corresponding to parameter r
MH_uniform = function(n, param_list, seeds){
  set.seed(seeds)
  r = param_list$r
  xs = matrix(vector(length=(n*ndim)), ndim, n)
  accs = vector(length = n)
  x = sample(-nsize:nsize, ndim, replace = TRUE)
  ngrids = matrix(rep(-r:r, ndim), ndim, (2*r+1), byrow = TRUE)
  old_energy = energy(x)
  for(i in 1:n){
    condition_matrix <- x + ngrids
    valid_moves <- (condition_matrix >= -nsize) & (condition_matrix <= nsize)
    probs <- ifelse(valid_moves, 1, 0)
    sampled_index <- sample(which(as.vector(probs) == 1), 1)
    
    entry <- ((sampled_index - 1) %% ndim) + 1
    change <- ((sampled_index - 1) %/% ndim) - r  # Adjust change correctly
    
    xnew <- x
    xnew[entry] <- x[entry] + change
    if((-10 <= xnew[entry]) & (xnew[entry] <= 10)){
        new_energy = energy(xnew)
        ratio = exp(new_energy-old_energy)
        if (runif(1)<=ratio){
          x = xnew
          xs[,i] = x
          accs[i] = 1
          old_energy = new_energy
        }
        else {
          xs[,i] = x
          accs[i] = 0
        }
    }else{
      xs[,i] = x
      accs[i] = 0
    }
  }
  return(list(xs = xs, accs=accs))
}

#-- NCG
#-- param_list: a list with elements 
#     a, corresponding to parameter delta
NCG = function(n, param_list, seeds){
  set.seed(seeds)
  delta = param_list$a
  xs = matrix(vector(length=(n*ndim)), ndim, n)
  accs = vector(length = n)
  x = sample(-nsize:nsize, ndim, replace = TRUE)
  dimsize = 2*nsize+1
  ngrids = matrix(rep(grid, ndim), ndim , length(grid), byrow=TRUE)
  old_energy = energy(x)
  oldgrad = gradient(x)
  for(i in 1:n){
    energys = 0.5*matrix(rep(oldgrad, dimsize), ndim, dimsize, byrow = FALSE)*ngrids-(1/(2*delta))* sweep(-ngrids, 1, -as.vector(x))^2
    ps = exp(energys)/rowSums(exp(energys))
    xnew = vector(length = ndim)
    for(j in 1:ndim){
      w = runif(1)
      xnew[j] = utodis(w, cumsum(ps[j,]))-(nsize+1)
    }
    newgrad = gradient(xnew)
    new_energy = energy(xnew)
    qdown = prod((exp(0.5*oldgrad*xnew-(1/(2*delta))*(x-xnew)^2))/rowSums(exp(cbind(0.5*matrix(rep(oldgrad, dimsize), ndim, dimsize, byrow = FALSE)*ngrids-(1/(2*delta))* sweep(-ngrids, 1, -as.vector(x))^2))))
    qup = prod((exp(0.5*newgrad*x-(1/(2*delta))*(x-xnew)^2))/rowSums(exp(cbind(matrix(0.5*rep(newgrad, dimsize), ndim, dimsize, byrow = FALSE)*ngrids-(1/(2*delta))* sweep(-ngrids, 1, -as.vector(xnew))^2))))
    ratio = exp(new_energy-old_energy)
    Q = (qup/qdown)*ratio
    if (runif(1)<=Q){
      x = xnew
      xs[,i] = xnew
      accs[i] = 1
      old_energy = new_energy
      oldgrad = newgrad
    }
    else{
      xs[,i] = x
      accs[i] = 0
    }
  }
  return(list(xs=xs, accs=accs))
}

#-- Ordinal GWG
#-- param_list: a list with parameter
#    r, corresponding to parameter r
GWG = function(n, param_list, seed){
  set.seed(seed)
  r = param_list$r
  xs = matrix(vector(length=(n*ndim)), ndim, n)
  accs = vector(length = n)
  x = sample(-nsize:nsize, ndim, replace = TRUE)
  ngrids = matrix(rep(-r:r, ndim), ndim, (2*r+1), byrow = TRUE)
  old_energy = energy(x)
  oldgrad = gradient(x)
  for(i in 1:n){
      probs = exp(matrix(rep(oldgrad, (2*r+1)), ndim, (2*r+1), byrow = FALSE)*ngrids)
      condition_matrix = x+ ngrids
      probs[condition_matrix < -nsize | condition_matrix > nsize] = 0
      probs_vector <- as.vector(probs)
      cum_probs <- cumsum(probs_vector)
      random_value <- runif(1, min = 0, max = tail(cum_probs, 1))
      sampled_index <- which(cum_probs >= random_value)[1]
      entry <- ((sampled_index - 1) %% ndim) + 1
      change <- ((sampled_index - 1) %/% ndim) + 1
      if((-10 <= x[entry]+change-(r+1)) & (x[entry]+change-(r+1) <= 10)){
        xnew = x
        xnew[entry] = x[entry] +change - (r+1)
        newgrad = gradient(xnew)
        new_energy = energy(xnew)
        back_probs = exp(matrix(rep(newgrad, (2*r+1)), ndim, (2*r+1), byrow = FALSE)*ngrids)
        condition_matrix = xnew + ngrids
        back_probs[condition_matrix < -nsize | condition_matrix > nsize] = 0
        back_probs_vector <- as.vector(back_probs)
        ratio = exp(new_energy-old_energy+(newgrad[entry]+oldgrad[entry])*((r+1)-change))*sum(probs_vector)/sum(back_probs_vector)
        if (runif(1)<=ratio){
          x = xnew
          xs[,i] = xnew
          accs[i] = 1
          old_energy = new_energy
          oldgrad = newgrad
        }
        else{
          xs[,i] = x
          accs[i] = 0
        }
      }
      else{
        xs[,i] = x
        accs[i] = 0
      }
  }
  return(list(xs=xs, accs=accs))
}

#-- AVG
#-- param_list: a list with elements 
#     a, corresponding to parameter delta
AVG = function(n, param_list, seeds){
  set.seed(seeds)
  delta = param_list$a
  xs = matrix(vector(length=(n*ndim)), ndim, n)
  accs = vector(length = n)
  x = sample(-nsize:nsize, ndim, replace = TRUE)
  dimsize = 2*nsize+1
  ngrids = matrix(rep(grid, ndim), ndim , length(grid), byrow=TRUE)
  old_energy = energy(x)
  oldgrad = gradient(x)
  for(i in 1:n){
    y = sqrt(2/delta)*x+rnorm(ndim)
    energys = matrix(rep(oldgrad, dimsize), ndim, dimsize, byrow = FALSE)*ngrids-0.5* sweep(-sqrt(2/delta)*ngrids, 1, -as.vector(y))^2
    ps = (1/rowSums(exp(energys)))*exp(energys)
    xnew = vector(length = ndim)
    for(j in 1:ndim){
      w = runif(1)
      xnew[j] = utodis(w, cumsum(ps[j,]))-(nsize+1)
    }
    newgrad = gradient(xnew)
    new_energy = energy(xnew)
    qdown = prod((exp(oldgrad*xnew-0.5*(y-sqrt(2/delta)*xnew)^2))/rowSums(exp(cbind(matrix(rep(oldgrad, dimsize), ndim, dimsize, byrow = FALSE)*ngrids-0.5* sweep(-sqrt(2/delta)*ngrids, 1, -as.vector(y))^2))))
    qup = prod((exp(newgrad*x-0.5*(y-sqrt(2/delta)*x)^2))/rowSums(exp(cbind(matrix(rep(newgrad, dimsize), ndim, dimsize, byrow = FALSE)*ngrids-0.5* sweep(-sqrt(2/delta)*ngrids, 1, -as.vector(y))^2))))
    
    ratio = exp(new_energy-old_energy-0.5*sum((y-sqrt(2/delta)*xnew)^2)+0.5*sum((y-sqrt(2/delta)*x)^2))
    Q = (qup/qdown)*ratio
    if (runif(1)<=Q){
      x = xnew
      xs[,i] = xnew
      accs[i] = 1
      old_energy = new_energy
      oldgrad = newgrad
    }
    else{
      xs[,i] = x
      accs[i] = 0
    }
  }
  return(list(xs=xs, accs=accs))
}

#-- Discrete Hams
#-- param_list: a list with elements 
#     a, corresponding to parameter delta
#     b, corresponding to parameter epsilon
#     gamma, corresponding to parameter beta
#     phi, corresponding to parameter phi
DHAMS =  function(n, param_list, seeds){
  set.seed(seeds)
  delta = param_list$a
  eps = param_list$b
  gamma = param_list$gamma
  phi = param_list$c
  dimsize = (2*nsize)+1
  xs = matrix(vector(length=(n*ndim)), ndim, n)
  us = matrix(vector(length=(n*ndim)), ndim, n)
  accs = vector(length = n)
  s = sample(-nsize:nsize, ndim, replace = TRUE)
  u = rnorm(ndim)
  ngrids = matrix(rep(grid, ndim), ndim , length(grid), byrow=TRUE)
  oldgrad = gradient(s)
  for(i in 1:n){
    u12 = eps*u+sqrt(1-eps^2)*rnorm(ndim)
    y = s - delta*u12
    energys = matrix(rep(oldgrad, dimsize), ndim, dimsize, byrow = FALSE)*ngrids-(1/(2*delta^2))* sweep(-ngrids, 1, -as.vector(y))^2
    pfore = (1/rowSums(exp(energys)))*exp(energys)
    snew = vector(length = ndim)
    for(j in 1:ndim){
      cumd = cumsum(pfore[j,])
      snew[j] = utodis((-distou((s[j]+(nsize+1)),cumd)+gamma*runif(1))%%1,cumd)-(nsize+1)
    }
    newgrad = gradient(snew)
    unew = (y-snew)/delta+as.vector(phi*(newgrad-oldgrad))
    energys = matrix(rep(newgrad, dimsize), ndim, dimsize, byrow = FALSE)*ngrids-(1/(2*delta^2))* sweep(-ngrids, 1, -as.vector(snew+delta*unew))^2
    pback = (1/rowSums(exp(energys)))*exp(energys)
    qprod = 1
    for (j in 1:ndim) {
      qprod = qprod*trans_prob(snew[j]+nsize+1, s[j]+nsize+1, cumsum(pback[j,]), -1, gamma)/trans_prob(s[j]+nsize+1, snew[j]+nsize+1, cumsum(pfore[j,]), -1, gamma)
    }
    log_ratio = energy(snew)-energy(s)-0.5*(sum(unew^2)-sum(u12^2))
    q  = qprod*exp(log_ratio)
    if (is.nan(q)){
      u = -u12
      xs[,i] = s
      us[,i] = -u12
      accs[i] = 0
    }else if (runif(1)<=q){
      s = snew
      u = unew
      xs[,i] = snew
      us[,i] = unew
      accs[i] = 1
      oldgrad = newgrad
    }
    else{
      u = -u12
      xs[,i] = s
      us[,i] = -u12
      accs[i] = 0
    }
  }
  return(list(xs=xs, us=us, accs=accs))
}
