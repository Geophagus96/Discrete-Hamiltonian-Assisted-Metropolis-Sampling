require(pracma)
require(parallel)
require(doParallel)
require(Rcpp)
require(RcppArmadillo)
sourceCpp('distou_revise.cpp')
params = readRDS('header_highdim.RData')
sourceCpp('sparse_linear_algos.cpp')

#--- params: the RData file that stores design matrix X, response y and hyper-parameter values

x = params$x
z = params$z
api = params$api
bpi = params$bpi
asigma = params$asigma
bsigma = params$bsigma
g = params$g
l = params$l

#--- nsample: number of draws per chain
#    nburn_in: number of burn-in draws
#    nplot: number of draws for trace plots
#    ngaps: interval (gap) between draws to calculate metrics
#    nrepeat: number of repeated chains
#    ncont: number of draws per chain continuing from previous draws

nburn_in = 5000
nsample = 29000
ngaps = 1000
nrepeat = 50
ndim = ncol(x)
ncont = 25000

#--- Evaluation Metrics calculation functions, including
#      average_flips: average flips
#      inclusion_prob_single: PIP

average_flips = function(xss){
    diffs = sapply(2:ncol(xss), function(i) sum(abs(xss[, i] - xss[, i - 1])))
    return(mean(diffs))
}

inclusion_prob_single = function(xss, dim, sample_interval){
    sumi = vector(length = length(sample_interval))
    for(k in 1:length(sample_interval)){
        if(k == 1){
            sumi[k] = sum(xss[dim, 1:sample_interval[1]])
        }
        else {
            sumi[k] = sumi[k-1]+sum(xss[dim, (sample_interval[k-1]+1):sample_interval[k]])
        }
    }
    return(sumi/sample_interval)
}

ESS_2 = function(xss_multi){
  ndim = dim(xss_multi)[2]
  n = dim(xss_multi)[3]
  m = dim(xss_multi)[1]
  ess_2 = vector(length = ndim)
  for(i in 1:ndim){
    xss = xss_multi[,i,]
    row_mean = rowMeans(xss)
    W = 1/(m*(n-1))*sum(sweep(xss, 1, row_mean)^2)
    B = n/(m-1)*sum((row_mean-mean(xss))^2)
    ess_2[i] = n*W/B
  }
  return(ess_2)
}

#--- The following functions are used for generating repeated chains
#--- adaptive_tuning_grad_single: a function for repeated chains with random initilization
#    Arguments:  
#       method: sampler
#       method_name: display name for that sampler
#       param_list: list of parameters
#       ndim: dimension of the distribution (dimension of design matrix)
#       nsample: number of draws
#       ngaps: interval (gap) between number of draws
#--- adaptive_tuning_grad_continue: a function for repeated chains continuing from previous simulation results
#    Arguments:  
#       method: sampler
#       method_name: display name for that sampler
#       param_list: list of parameters
#       ndim: dimension of the distribution (dimension of design matrix)
#       nsample: number of draws
#       ngaps: interval (gap) between number of draws
#       further: indicating whether continuing from the first part of simulation results or the second part


adaptive_tuning_grad_single = function(method, method_name, param_list, ndim, nsample, ngaps){
    sample_interval = ngaps*(1:(nsample/ngaps)) 

    no_cores <- detectCores() - 1

    print('MCMC RUNNING')
    
    registerDoParallel(core=no_cores)
    set.seed(123)
    s0s = matrix(rbinom(nrepeat * ncol(x), size = 1, prob = 0.5), nrow = nrepeat, ncol = ncol(x))

    results = foreach(j = 1:nrepeat) %dopar% method(nsample, param_list, j, x, z, s0s[j,], api, bpi, asigma, bsigma, g, l)
    print('MCMC COMPLETED')

    xss = array(1, dim = c(nrepeat, ndim, nsample))
    accs  =  0*vector(length = nrepeat)
    for(j in 1:nrepeat){
        xss[j,,] = results[[j]]$ss
        accs[j] = mean(results[[j]]$accs)
    }
    
    registerDoParallel(core=no_cores)
    include_prob1 = simplify2array(foreach(i = 1:nrepeat) %dopar% inclusion_prob_single(xss[i,,], 1, sample_interval))
    ave_include_prob1 = rowMeans(include_prob1)

    registerDoParallel(core=no_cores)
    include_prob601 = simplify2array(foreach(i = 1:nrepeat) %dopar% inclusion_prob_single(xss[i,,], 601, sample_interval))
    ave_include_prob601 = rowMeans(include_prob601)
    
    registerDoParallel(core=no_cores)
    ess2s = ESS_2(xss)

    cat(paste(method_name,'ESS for x1'))
    print(ess2s[1])

    cat(paste(method_name,'ESS for x601'))
    print(ess2s[601])

    cat(paste(method_name, 'accs'))
    print(mean(accs))
 
    cat(paste(method_name, 'pip of x1'))
    print(ave_include_prob1)

    cat(paste(method_name, 'pip of x601'))
    print(ave_include_prob601)

    saveRDS(list(xss=xss, ess2s=ess2s, accs = mean(accs), pip1 = ave_include_prob1, pip2 = ave_include_prob601, pip4=ave_include_prob4, pip5=ave_include_prob604), file = paste('adaptive_tuning', method_name, '.RData', sep=''))
}

adaptive_tuning_grad_continue = function(method, method_name, param_list, ndim, nsample, ngaps, further=FALSE, good= 'good'){
    sample_interval = ngaps*(1:(nsample/ngaps)) 
    if(further == FALSE){
        prev_data = readRDS(paste('adaptive_tuning',method_name, '.RData', sep=''))
    }
    else{
        prev_data = readRDS(paste('adaptive_tuning',method_name, 'cont', '.RData', sep=''))
        nsample = nsample
    }
    s0s = matrix(vector(length = nrepeat*ncol(x)), nrow = nrepeat, ncol = ncol(x))
    for(i in 1:nrepeat){
        s0s[i, ] = prev_data$xs[i,,dim(prev_data$xs)[3]]
    }
    no_cores <- detectCores() - 1

    print('Reinitialization Completed')
    print('MCMC RUNNING')
    
    registerDoParallel(core=no_cores)


    results = foreach(j = 1:nrepeat) %dopar% method(nsample, param_list, j, x, z, s0s[j,], api, bpi, asigma, bsigma, g, l)
    print('MCMC COMPLETED')

    xss = array(1, dim = c(nrepeat, ndim, nsample))
    accs  =  0*vector(length = nrepeat)
    for(j in 1:nrepeat){
        xss[j,,] = results[[j]]$ss
        accs[j] = mean(results[[j]]$accs)
    }
    
    registerDoParallel(core=no_cores)
    include_prob1 = simplify2array(foreach(i = 1:nrepeat) %dopar% inclusion_prob_single(xss[i,,], 1, sample_interval))
    ave_include_prob1 = rowMeans(include_prob1)

    registerDoParallel(core=no_cores)
    include_prob601 = simplify2array(foreach(i = 1:nrepeat) %dopar% inclusion_prob_single(xss[i,,], 601, sample_interval))
    ave_include_prob601 = rowMeans(include_prob601)

    registerDoParallel(core=no_cores)
    ess2s = ESS_2(xss)

    cat(paste(method_name,'ESS for x1'))
    print(ess2s[1])

    cat(paste(method_name,'ESS for x601'))
    print(ess2s[601])

    cat(paste(method_name, 'accs'))
    print(mean(accs))
 
    cat(paste(method_name, 'pip of x1'))
    print(ave_include_prob1)

    cat(paste(method_name, 'pip of x601'))
    print(ave_include_prob601)

    if(further == FALSE){
        saveRDS(list(xss=xss, ess2s=ess2s, accs = mean(accs), pip1 = ave_include_prob1, pip2 = ave_include_prob601, pip4=ave_include_prob4, pip5=ave_include_prob604), file = paste('adaptive_tuning', method_name, 'cont.RData', sep=''))
    }
    else {
        saveRDS(list(xss=xss, ess2s=ess2s, accs = mean(accs), pip1 = ave_include_prob1, pip2 = ave_include_prob601, pip4=ave_include_prob4, pip5=ave_include_prob604), file = paste('adaptive_tuning', method_name, 'fcont_', good, '.RData', sep=''))
    }
}

#--- Given the size of the simulation, we separate it into three parts,
#     each containing one part of the sample

param_lists_ncg = list()
param_ncg = c(0.136)
for(i in 1:length(param_ncg)){
  param_lists_ncg[[i]] = list(a=param_ncg[i])
}
adaptive_tuning_grad_single(NCG, 'NCG', param_lists_ncg[[1]], ndim, nsample, ngaps)
adaptive_tuning_grad_continue(NCG, 'NCG', param_lists_ncg[[1]], ndim, nsample, ngaps, further = FALSE)
adaptive_tuning_grad_continue(NCG, 'NCG', param_lists_ncg[[1]], ndim, nsample, ngaps, TRUE, 'good')

param_lists_avg = list()
param_avg = c(0.159)
for(i in 1:length(param_avg)){
  param_lists_avg[[i]] = list(a=param_avg[i])
}
adaptive_tuning_grad_single(AVG, 'AVG', param_lists_avg[[1]], ndim, nsample, ngaps)
adaptive_tuning_grad_continue(AVG, 'AVG', param_lists_avg[[1]], ndim, nsample, ngaps, further = FALSE)
adaptive_tuning_grad_continue(AVG, 'AVG', param_lists_avg[[1]], ndim, nsample, ngaps, TRUE, 'good')

param_lists_hams = list()
param_lists_hams[[1]] = list(a=0.283, b=0.85, c = 1.0, d=0.0)
adaptive_tuning_grad_single(DHAMS, 'OutsideHams', param_lists_hams[[1]], ndim, nsample, ngaps)
adaptive_tuning_grad_continue(DHAMS, 'OutsideHams', param_lists_hams[[1]], ndim, nsample, ngaps, further=FALSE)
adaptive_tuning_grad_continue(DHAMS, 'OutsideHams', param_lists_hams[[1]], ndim, nsample, ngaps, TRUE, 'good')

param_lists_hams = list()
param_lists_hams[[1]] = list(a=0.260, b=0.85, c = 0.1, d= 0.0)
adaptive_tuning_grad_single(DHAMS, 'OverrelaxedOutsideHams', param_lists_hams[[1]], ndim, nsample, ngaps)
adaptive_tuning_grad_continue(DHAMS, 'OverrelaxedOutsideHams',  param_lists_hams[[1]], ndim, nsample, ngaps, further=FALSE)
adaptive_tuning_grad_continue(DHAMS, 'OverrelaxedOutsideHams', param_lists_hams[[1]], ndim, nsample, ngaps, TRUE, 'good')
