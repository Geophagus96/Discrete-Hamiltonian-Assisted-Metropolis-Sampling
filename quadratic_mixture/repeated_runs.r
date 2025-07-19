require(pracma)
require(parallel)
require(doParallel)
source('algos.R')
prob_tables = readRDS('probtable.RData')

#--- nsample: number of draws per chain
#    nburn_in: number of burn-in draws
#    nplot: number of draws for trace plots
#    ngaps: interval (gap) between draws to calculate metrics
#    nrepeat: number of repeated chains
#    nsize: (2*nsize+1) is the lattice size

nsample = 24000
nburn_in = 500
nplot = 150
ngaps = 400
nrepeat = 100
nsize = 10

prob_table_one = prob_tables[[1]]
prob_table_two = prob_tables[[2]]
prob_table_four = prob_tables[[3]]
mus = prob_tables[[4]]
Ws = prob_tables[[5]]

grid_int_transformer2 = function(x, dimsize){
  bins = nsize+x
  powers = (2*nsize+1)^(1:dimsize-1)
  return(sum(bins*powers)+1)
}

int_grid_transformer2 = function(intx, dimsize){
  x = vector(length = dimsize)
  for(i in 1:dimsize){
    x[i] = intx%%(2*nsize+1)
    intx = intx%/%(2*nsize+1)
  }
  return(x)
}

e2 = 0
e12 = 0
for(i in 1:(2*nsize+1)^2){
  x = int_grid_transformer2(i-1, 2)-nsize
  e2 = e2+x[1]^2*prob_table_two[i]
  e12 = e12+x[1]*x[2]*prob_table_two[i]
}

RowVar <- function(x, ...) {
  rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
}

# --- Functions to calculate metrics for evaluating samplers
#
# --- tv1_single_param:
#     Computes the total variation (TV) distance between the one-dimensional
#     true distribution and its corresponding empirical distribution from samples.
#
# --- Arguments:
#     xss             : Input sampling data
#     idcomb          : Indices of the selected pair of dimensions for the marginal
#     sample_interval : Interval of number of draws used for calculating the TV distance
#
# --- tv2_single_param:
#     Computes the total variation (TV) distance between the two-dimensional
#     true distribution and its corresponding empirical distribution from samples.
#
# --- Arguments:
#     xss             : Input sampling data
#     idcomb          : Indices of the selected pair of dimensions for the marginal
#     sample_interval : Interval of number of draws used for calculating the TV distance
#
# --- tv4_single_param:
#     Computes the total variation (TV) distance between the four-dimensional
#     true distribution and its corresponding empirical distribution from samples.
#
# --- Arguments:
#     xss             : Input sampling data
#     idcomb          : Indices of the selected combination of four dimensions for the marginal
#     sample_interval : Interval of number of draws used for calculating the TV distance



tv1_single_param = function(xss, iddim, sample_interval){
  tv1 = vector(length=length(sample_interval))
  int_xs = nsize+xss[iddim,]
  for (k in 1:length(sample_interval)) {
    freq = table(int_xs[1:sample_interval[k]])
    b = rep(0, (2*nsize+1))
    b[1+as.numeric(names(freq))] = freq/sample_interval[k]
    tv1[k]=0.5*sum(abs(prob_table_one-b))
  }
  return(tv1)
}

tv2_single_param = function(xss, idcomb, sample_interval){
  tv2 = vector(length=length(sample_interval))
  int_xs = apply(xss[idcomb, ], 2, function(x){grid_int_transformer2(x, 2)})
  for (k in 1:length(sample_interval)) {
      freq = table(int_xs[1:sample_interval[k]])
      b = rep(0, (2*nsize+1)^2)
      b[as.numeric(names(freq))] = freq/sample_interval[k]
      tv2[k] = 0.5*sum(abs(prob_table_two-b))
  }
  return(tv2)
}

tv4_single_param = function(xss, idcomb, sample_interval){
  tv4 = vector(length = length(sample_interval))
  int_xs = apply(xss[idcomb, ], 2, function(x){grid_int_transformer2(x, 4)})
  for (k in 1:length(sample_interval)) {
      freq = table(int_xs[1:sample_interval[k]])
      b = rep(0, (2*nsize+1)^4)
      b[as.numeric(names(freq))] = freq/sample_interval[k]
      tv4[k] = 0.5*sum(abs(prob_table_four-b))
  }
  return(tv4)
}

ex_single_param_custom = function(xss, dimid, sample_interval){
  ex = vector(length = length(sample_interval))
  for(k in 1:length(sample_interval)){
    ex[k] = mean(xss[dimid,1:sample_interval[k]])
  }
  return(ex)
}

ex2_single_param_custom = function(xss, dimid, sample_interval){
  ex2 = vector(length = length(sample_interval))
  for(k in 1:length(sample_interval)){
    ex2[k] = mean(xss[dimid,1:sample_interval[k]]^2)
  }
  return(ex2)
}

ex12_single_param_custom = function(xss, dimids, sample_interval){
  ex12 = vector(length = length(sample_interval))
  for(k in 1:length(sample_interval)){
    ex12[k] = mean(xss[dimids[1],1:sample_interval[k]]*xss[dimids[2],1:sample_interval[k]])
  }
  return(ex12)
}

# --- Function to execute repeated runs and return a dataset of sampler outputs with evaluation metrics
#
# --- Returned:
#     xss_best       : Array of samples with shape 
#                      (number of repeated runs) × (distribution dimension) × (number of draws)
#     accs_best      : Average acceptance rates
#     ave_tv1        : Mean total variation (TV) distance for one-dimensional marginal distributions
#     ave_tv2        : Mean TV distance for two-dimensional marginal distributions
#     ave_tv4        : Mean TV distance for four-dimensional marginal distributions
#     ave_tv1_var    : Standard deviation of TV distance for one-dimensional marginals
#     ave_tv2_var    : Standard deviation of TV distance for two-dimensional marginals
#     ave_tv4_var    : Standard deviation of TV distance for four-dimensional marginals
#     ave_ex1_bias   : Mean bias in estimating E[x_i]
#     ave_ex1_var    : Variance in estimating E[x_i]
#     ave_ex2_bias   : Mean bias in estimating E[x_i²]
#     ave_ex2_var    : Variance in estimating E[x_i²]
#     ave_ex12_bias  : Mean bias in estimating E[x_i x_j]
#     ave_ex12_var   : Variance in estimating E[x_i x_j]
#
# --- Arguments:
#     method         : The sampler method to run
#     method_name    : Display name of the chosen sampler
#     best_param     : List of input parameters for the sampler
#     ndim           : Dimension of the target distribution
#     nsample        : Number of draws after burn-in
#     nburn_in       : Number of burn-in draws to discard
#     ngaps          : Interval (gap) between samples used for metric calculations

repeated_runs = function(method, method_name, param_lists, ndim, nsample, nburn_in, ngaps){

    no_cores = detectCores()-1
    sample_interval = ngaps*(1:(nsample/ngaps)) 

    xss_best = array(1, dim = c(nrepeat, length(param_lists), ndim, nsample))
    accs_best  =  matrix(rep(0, nrepeat*length(param_lists)), nrow = nrepeat, ncol = length(param_lists))
    results = list()
    for(i in 1:length(param_lists)){
        registerDoParallel(cores=no_cores)
        results[[i]] = foreach(j = 1:nrepeat) %dopar% method((nburn_in+nsample), param_lists[[i]], j)
    }
    for(j in 1:nrepeat){
        for(i in 1:length(param_lists)){
            xss_best[j,i,,] = results[[i]][[j]]$xs[, ((nburn_in+1):(nburn_in+nsample))]
            accs_best[j,i] = mean(results[[i]][[j]]$accs[((nburn_in+1):(nburn_in+nsample))])
        }
    }


    combinations = combn(1:ndim, 2)
    ncomb = ncol(combinations)

    ave_ex1_bias = 0.0*vector(length = length(sample_interval))
    ave_ex1_var = 0.0*vector(length = length(sample_interval))
    ave_ex2_bias = 0.0*vector(length = length(sample_interval))
    ave_ex2_var = 0.0*vector(length = length(sample_interval))
    ave_ex12_bias = 0.0*vector(length = length(sample_interval))
    ave_ex12_var = 0.0*vector(length = length(sample_interval))
    for(j in 1:ndim){
        ex1_best = simplify2array(foreach(i = 1:nrepeat) %dopar% ex_single_param_custom(xss_best[i,,,], j, sample_interval))
        ave_ex1_bias_best = rowMeans(ex1_best)
        ave_ex1_bias_best = (ave_ex1_bias_best)^2
        ave_ex1_var_best = RowVar(ex1_best)
        ave_ex1_bias = ave_ex1_bias+ave_ex1_bias_best
        ave_ex1_var = ave_ex1_var+ave_ex1_var_best
    }
    ave_ex1_bias = ave_ex1_bias/ndim
    ave_ex1_var = ave_ex1_var/ndim

    for(j in 1:ndim){
        ex2_best = simplify2array(foreach(i = 1:nrepeat) %dopar% ex2_single_param_custom(xss_best[i,,,], j, sample_interval))
        ave_ex2_bias_best = rowMeans(ex2_best)
        ave_ex2_bias_best = (ave_ex2_bias_best-e2)^2
        ave_ex2_var_best = RowVar(ex2_best)
        ave_ex2_bias = ave_ex2_bias+ave_ex2_bias_best
        ave_ex2_var = ave_ex2_var+ave_ex2_var_best
    }
    ave_ex2_bias = ave_ex2_bias/ndim
    ave_ex2_var = ave_ex2_var/ndim

    for(j in 1:ncomb){
        ex12_best = simplify2array(foreach(i = 1:nrepeat) %dopar% ex12_single_param_custom(xss_best[i,,,], combinations[,j], sample_interval))
        ave_ex12_bias_best = rowMeans(ex12_best)
        ave_ex12_bias_best = (ave_ex12_bias_best-e12)^2
        ave_ex12_var_best = RowVar(ex12_best)
        ave_ex12_bias = ave_ex12_bias+ave_ex12_bias_best
        ave_ex12_var = ave_ex12_var+ave_ex12_var_best
    }
    ave_ex12_bias = ave_ex12_bias/ncomb
    ave_ex12_var = ave_ex12_var/ncomb

    tv1s = array(0.0, dim=c(ndim, length(sample_interval), nrepeat))
    ave_tv1 = 0.0*vector(length = length(sample_interval))
    ave_tv1_var = 0.0*vector(length = length(sample_interval))
    for(j in 1:ndim){
        tv1_best = simplify2array(foreach(i = 1:nrepeat) %dopar% tv1_single_param(xss_best[i,,,], j, sample_interval))
        tv1s[j,,] = tv1_best
        ave_tv1 = ave_tv1 + rowMeans(tv1_best)
        ave_tv1_var = ave_tv1_var + sqrt(RowVar(tv1_best))
    }
    ave_tv1 = ave_tv1/ndim
    ave_tv1_var = ave_tv1_var/ndim

    tv2s = array(0.0, dim=c(ncomb, length(sample_interval), nrepeat))
    ave_tv2 = 0.0*vector(length = length(sample_interval))
    ave_tv2_var = 0.0*vector(length = length(sample_interval))
    for(j in 1:ncomb){
        tv2_best = simplify2array(foreach(i = 1:nrepeat) %dopar% tv2_single_param(xss_best[i,,,], combinations[,j], sample_interval))
        tv2s[j,,] = tv2_best
        ave_tv2 = ave_tv2 + rowMeans(tv2_best)
        ave_tv2_var = ave_tv2_var + sqrt(RowVar(tv2_best))
    }
    ave_tv2 = ave_tv2/ncomb
    ave_tv2_var = ave_tv2_var/ncomb

    combinations_4 = combn(1:ndim, 4)
    ncomb_4 = ncol(combinations_4)
    tv4s = array(0.0, dim=c(ncomb_4, length(sample_interval), nrepeat))
    ave_tv4 = 0.0*vector(length = length(sample_interval))
    ave_tv4_var = 0.0*vector(length = length(sample_interval))
    for(j in 1:ncomb_4){
      tv4_best = simplify2array(foreach(i = 1:nrepeat) %dopar% tv4_single_param(xss_best[i,,,], combinations_4[,j], sample_interval))
      tv4s[j,,] = tv4_best
      ave_tv4 = ave_tv4 + rowMeans(tv4_best)
      ave_tv4_var = ave_tv4_var + sqrt(RowVar(tv4_best))
    }
    ave_tv4 = ave_tv4/ncomb_4
    ave_tv4_var = ave_tv4_var/ncomb_4

    saveRDS(list(xss_best, colMeans(accs_best), ave_tv1, ave_tv2, ave_ex1_bias, ave_ex1_var, ave_ex2_bias, ave_ex2_var, ave_ex12_bias, ave_ex12_var, ave_tv1_var, ave_tv2_var, ave_tv4, ave_tv4_var, tv1s, tv2s, tv4s, param_lists), file=paste('disgau', method_name, '.RData', sep=''))

}

MHU_param = list(list(r=4))
repeated_runs(MH_uniform, 'MHU', MHU_param, 8, nsample, nburn_in, ngaps)

GWG_param = list(list(r=2))
repeated_runs(GWG, 'GWG', GWG_param, 8, nsample, nburn_in, ngaps)

NCG_param = list(list(a= 3.3))
repeated_runs(NCG, 'NCG', NCG_param, 8, nsample, nburn_in, ngaps)

AVG_param = list(list(a=1.86))
repeated_runs(AVG, 'AVG', AVG_param, 8, nsample, nburn_in, ngaps)

VDHAMS_param = list(list(a=1.07, b=0.9, gamma = 1.0, c = -0.5))
repeated_runs(DHAMS, 'V-DHAMS', VDHAMS_param, 8, nsample, nburn_in, ngaps)

ODHAMS_param = list(list(a=0.77, b=0.9, gamma=0.1, c=-0.7))
repeated_runs(DHAMS, 'O-DHAMS', ODHAMS_param, 8, nsample, nburn_in, ngaps)