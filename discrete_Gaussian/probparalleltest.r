require(parallel)
require(doParallel)

rho = 0.9
sigma= 5
nsize = 10
ndim = 8
grid = -nsize:nsize 
cov_mat = sigma^2*(rho*matrix(rep(1, ndim^2), ndim, ndim)+(1-rho)*diag(rep(1, ndim)))
inv_cov_mat = solve(cov_mat)

print(paste('discrete variance', sigma))
print(paste('discrete correlation', rho))
grid_int_transformer = function(x){
  bins = nsize+x
  powers = (2*nsize+1)^(1:ndim-1)
  return(sum(bins*powers))
}

grid_int_transformer2 = function(x, dimsize){
  bins = nsize+x
  powers = (2*nsize+1)^(1:dimsize-1)
  return(sum(bins*powers)+1)
}

int_grid_transformer = function(intx){
  x = vector(length = ndim)
  for(i in 1:ndim){
    x[i] = intx%%(2*nsize+1)
    intx = intx%/%(2*nsize+1)
  }
  return(x)
}

int_grid_transformer2 = function(intx, dimsize){
  x = vector(length = dimsize)
  for(i in 1:dimsize){
    x[i] = intx%%(2*nsize+1)
    intx = intx%/%(2*nsize+1)
  }
  return(x)
}

gradient = function(x){
  return(-inv_cov_mat%*%x)
}

energy = function(x) {
  return(-0.5*t(x)%*%inv_cov_mat%*%x)
}

no_cores <- detectCores() - 1


#--- Functions to calculate probability tables
#    prob_table_one: probability table of one-dimensional marginal distribution
#    prob_table_two: probability table of two-dimensional marginal distribution
#    prob_table_four: probability table of four-dimensional marginal distribution

registerDoParallel(cores=no_cores)
start_time = Sys.time()
prob_table_four = foreach(i=1:((2*nsize+1)^4),.combine = 'c') %dopar% {
  ans = 0 
  for(j in 1:((2*nsize+1)^(ndim-4))){
    x = c(int_grid_transformer2(i-1,4), int_grid_transformer2(j-1, ndim-4)) -nsize
    ans = ans +exp(energy(x))
  }
  ans
  
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

saveRDS(list(prob_table_one/sum(prob_table_one), prob_table_two/sum(prob_table_two), prob_table_four), file='probtable.RData')