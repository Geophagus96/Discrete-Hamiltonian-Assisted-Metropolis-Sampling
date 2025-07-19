require(abind)
require(pracma)
require(parallel)
require(doParallel)
require(Rcpp)
require(RcppArmadillo)
ngaps = 1000
nrepeat = 50
no_cores <- detectCores() - 1

sample_interval = ngaps*(1:(65000/ngaps)) 
method_names = c('AVG', 'NCG', 'OutsideHams', 'OverrelaxedOutsideHams')
display_names = c('AVG', 'NCG', 'V-DHAMS','O-DHAMS')



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

#--- Computation of PIP related metrics
#    ave_include_prob_cov1: standard deviation of PIP for s_1
#    ave_include_prob_cov1: standard deviation of PIP for s_1
#    saved to PIP_ds.RData
ave_include_prob_cov1 = matrix(vector(length = 4*length(sample_interval)), 4, length(sample_interval))
ave_include_prob_cov2 = matrix(vector(length = 4*length(sample_interval)), 4, length(sample_interval))
for(i in 1:4){
    print(paste(method_names[i], ' good results'))
    
    repart1 = readRDS(paste('adaptive_tuning', method_names[i], '.RData', sep=''))
    repart2 = readRDS(paste('adaptive_tuning', method_names[i], 'cont.RData', sep=''))
    xss = abind(repart1$xss, repart2$xss, along = 3)
    accs = (repart1$accs+repart2$accs)

    rm(repart1)
    rm(repart2)
    repart3 = readRDS(paste('adaptive_tuning', method_names[i], 'fcont_good.RData', sep=''))
    xss = abind(xss, repart3[[1]][,,1:29000], along = 3)

    xss = xss[,,8001:73000]
    accs = accs+repart3$accs

    print(paste('acceptance probability',accs/3))
    registerDoParallel(core=no_cores)
    include_prob1 = simplify2array(foreach(i = 1:nrepeat) %dopar% inclusion_prob_single(xss[i,,], 1, sample_interval))
    ave_include_prob1_sd = apply(include_prob1, 1, sd)

    registerDoParallel(core=no_cores)
    include_prob601 = simplify2array(foreach(i = 1:nrepeat) %dopar% inclusion_prob_single(xss[i,,], 601, sample_interval))
    ave_include_prob601_sd = apply(include_prob601, 1, sd)

    rm(xss)
    ave_include_prob_cov1[i, ] = ave_include_prob1_sd
    ave_include_prob_cov2[i, ] = ave_include_prob601_sd
}
saveRDS(list(ave_include_prob_cov1=ave_include_prob_cov1, ave_include_prob_cov2=ave_include_prob_cov2), 'PIP_ds.RData')


#    ave_include_prob_cov1: mean of PIP for s_1
#    ave_include_prob_cov1: mean of PIP for s_1
#    saved to PIP.RData
ave_include_prob_cov1 = matrix(vector(length = 4*length(sample_interval)), 4, length(sample_interval))
ave_include_prob_cov2 = matrix(vector(length = 4*length(sample_interval)), 4, length(sample_interval))

for(i in 1:4){

    print(paste(method_names[i], ' good results'))
    
    repart1 = readRDS(paste('adaptive_tuning', method_names[i], '.RData', sep=''))
    repart2 = readRDS(paste('adaptive_tuning', method_names[i], 'cont.RData', sep=''))
    xss = abind(repart1$xss, repart2$xss, along = 3)

    rm(repart1)
    rm(repart2)
    repart3 = readRDS(paste('adaptive_tuning', method_names[i], 'fcont_good.RData', sep=''))
    xss = abind(xss, repart3[[1]][,,1:29000], along = 3)
    rm(repart3)
    xss = xss[,,8001:73000]
    registerDoParallel(core=no_cores)
    include_prob1 = simplify2array(foreach(i = 1:nrepeat) %dopar% inclusion_prob_single(xss[i,,], 1, sample_interval))
    ave_include_prob1 = rowMeans(include_prob1)

    registerDoParallel(core=no_cores)
    include_prob601 = simplify2array(foreach(i = 1:nrepeat) %dopar% inclusion_prob_single(xss[i,,], 601, sample_interval))
    ave_include_prob601 = rowMeans(include_prob601)

    ave_include_prob_cov1[i, ] = ave_include_prob1
    ave_include_prob_cov2[i, ] = ave_include_prob601
    ess2s = ESS_2(xss[,c(1, 601), ])
    print(paste('ESS for s_1', ess2s[1]))
    print(paste('ESS for s_601', ess2s[1]))

    rm(xss)
}
saveRDS(list(ave_include_prob_cov1=ave_include_prob_cov1, ave_include_prob_cov2=ave_include_prob_cov2), 'PIP.RData')
