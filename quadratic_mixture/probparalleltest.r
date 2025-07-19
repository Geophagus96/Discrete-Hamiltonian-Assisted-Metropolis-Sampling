require(parallel)
require(doParallel)
require(pracma)

sigma= 1.8
nsize = 10
ndim = 8

grid = -nsize:nsize 
ncomponent = 5
center = linspace(-7,7, ncomponent)

#--- Define the mean and variances for each mixture component
mus = list()
Ws = list()
for(i in 1:ncomponent){
  mus[[i]] = rep(center[i], ndim)
  Ws[[i]] = diag(rep(1/sigma^2, ndim))
}


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

no_cores <- detectCores() - 1
registerDoParallel(cores=no_cores)
mixture_weighting = function(mu, W, ndim){
  ans = 0
  for(j in 1:((2*nsize+1)^(ndim))){
    x = int_grid_transformer2(j-1, ndim) -nsize
    diff = x-mu
    mass = exp(-0.5 * t(diff) %*% W %*% diff)
    ans = ans + mass/10000
  }  
  return(ans)
}

sub_mus = list()
sub_Ws = list()
  
for(i in 1:ncomponent){
    sub_mus[[i]] = rep(center[i], 4)
    sub_Ws[[i]] = diag(rep(1/sigma^2, 4))
}

print('calculating mixture weights')
weights = simplify2array(foreach(i = 1:ncomponent) %dopar% {mixture_weighting(sub_mus[[i]], sub_Ws[[i]], 4)})
weights = weights/sum(weights)

weighted_mass <- function(x, mus, Ws, weights) {
  components <- mapply(function(mu, W, weight) {
    diff <- x - mu
    weight*exp(-0.5 * t(diff) %*% W %*% diff)
  }, mus, Ws, weights, SIMPLIFY = TRUE)
  
  return(sum(components))
}


#--- Functions to calculate probability tables
#    prob_table_one: probability table of one-dimensional marginal distribution
#    prob_table_two: probability table of two-dimensional marginal distribution
#    prob_table_four: probability table of four-dimensional marginal distribution

registerDoParallel(cores=no_cores)
start_time = Sys.time()
prob_table_four = foreach(i=1:((2*nsize+1)^4),.combine = 'c') %dopar% {
  x = int_grid_transformer2(i-1,4) -nsize
  weighted_mass(x, sub_mus, sub_Ws, weights)
}

prob_table_four = prob_table_four/sum(prob_table_four)
print(Sys.time()-start_time)

prob_table_two = 0*vector(length=(2*nsize+1)^2)
for(i in 1:(2*nsize+1)^4){
  prob_table_two[grid_int_transformer2(int_grid_transformer2(i-1, 4)[1:2]-nsize, 2)] = prob_table_two[grid_int_transformer2(int_grid_transformer2(i-1, 4)[1:2]-nsize, 2)]+prob_table_four[i]
}

prob_table_one = vector(length=(2*nsize+1))
for(i in 1:(2*nsize+1)){
  prob_table_one[i] = sum(prob_table_two[((i-1)*(2*nsize+1)+1):(i*(2*nsize+1))])
}

e12 = 0
for(j in 1:((2*nsize+1)^2)){
    
    x = int_grid_transformer2(j-1,2) -nsize
   
    e12 = e12 + x[1]*x[2]*prob_table_two[j]
}

e2 = 0
for(j in 1:(2*nsize+1)){
  e2 = e2+(j-nsize-1)^2*prob_table_one[j]
}

print(paste('moment estimations', e12, e2))


saveRDS(list(prob_table_one, prob_table_two, prob_table_four, mus, Ws), file='probtable.RData')