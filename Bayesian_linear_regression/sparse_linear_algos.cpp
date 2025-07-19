#include <RcppArmadillo.h>
#include <math.h>
#include <iostream>
#include <distou_revise.cpp>

using namespace Rcpp;
using namespace arma;
using namespace std;



// #########################################################
// // Codes for "Discrete Hamiltonian-Assisted Metropolis Sampling" for Bayesian Linear Regression
// #########################################################

// --- Implements the following MCMC algorithms:
// --- NCG, AVG, V-DHAMS, and O-DHAMS.
// --- log_quantity: function for querying gradient of f(s), with input
//        x: the X matrix
//        z: the response y
//        s: the binary mask
//        api: a_psi
//        bpi: b_psi
//        asigma: a_sigma
//        bsigma: b_sigma
//        g: g
//        l: l (our kappa is chosen as 0.999 and therefore is dropped)
//    the function return a list of f(s) and the gradient of f(s)

// --- In all samplers below, we have additional arguments
//         n: Number of draws
//         j: Random seed for reproducibility
//         param_list: List of sampler-specific parameters
//         s0: the initial value for s


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List log_quantity(const mat& x,  const vec& z, const vec& s, double api, double bpi, double asigma, double bsigma, double g, double l) {
    int n = x.n_rows; 
    int p = x.n_cols; 
    arma::ivec binmask = arma::conv_to<arma::ivec>::from(s > 0.0);
    arma::mat xs = x.cols(arma::find(binmask == 1));

    const int ds = xs.n_cols;    

    mat xxts = xs.t()*xs;
    double ps = R::lgammafn(sum(s) + api) + R::lgammafn(p - sum(s) + bpi);
    mat invmat1 = inv(l * eye(ds, ds) + xxts); 
    mat invmat2 = inv(l * eye(ds, ds) + (g + 1) * xxts); 
    mat gradxs_part1 = xs * invmat1 - (g + 1) * xs * invmat2;
    mat gradxs_part20 = z.t() * xs * invmat2;
    mat gradxs_part2 = (asigma + n / 2) * (2 * g * z * z.t() * xs * invmat2 - 2 * g * (1 + g) * xs*gradxs_part20.t()*gradxs_part20);
    double gradxs_part3 = 2 * bsigma + accu(square(z)) - g * as_scalar(z.t() * xs * invmat2 * xs.t() * z);
    mat gradxs = gradxs_part1 + (1 / gradxs_part3) * gradxs_part2;
    vec grads = arma::zeros<arma::vec>(s.n_elem);
    grads.elem(arma::find(binmask == 1)) = sum((gradxs % xs), 0);

    double add_up = (R::digamma(sum(s) + api) - R::digamma(p - sum(s) + bpi));
    grads += add_up;
    vec eigs = eig_sym(xxts); 
    double part1 = 0.5 * accu(log(l + eigs) - log(l + (g + 1) * eigs));
    double part3 = (-(2 * asigma + z.size()) / 2) * log(2 * bsigma + accu(square(z)) - as_scalar(g * z.t() * xs * invmat2 * xs.t() * z));
    double logmass = ps + part1 + part3;

    return List::create(Named("log_mass") = logmass, Named("grads") = grads);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double log_mass(const mat& x,  const vec& z, const vec& s, double api, double bpi, double asigma, double bsigma, double g, double l) {
  int n = x.n_rows; 
  int p = x.n_cols; 
  arma::ivec binmask = arma::conv_to<arma::ivec>::from(s > 0.0);
  arma::mat xs = x.cols(arma::find(binmask == 1));

  const int ds = xs.n_cols;    

  mat xxts = xs.t()*xs;
  double ps = R::lgammafn(sum(s) + api) + R::lgammafn(p - sum(s) + bpi);
  mat invmat1 = inv(l * eye(ds, ds) + xxts); 
  mat invmat2 = inv(l * eye(ds, ds) + (g + 1) * xxts); 
  vec eigs = eig_sym(xxts); 
  double part1 = 0.5 * accu(log(l + eigs) - log(l + (g + 1) * eigs));
  double part3 = (-(2 * asigma + z.size()) / 2) * log(2 * bsigma + accu(square(z)) - as_scalar(g * z.t() * xs * invmat2 * xs.t() * z));
  double logmass = ps + part1 + part3;

  return logmass;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List NCG(int n, List param_list, int j, mat& x,  const vec& z, const vec &s0, double api, double bpi, double asigma, double bsigma, double g, double l) {
  
  arma_rng::set_seed(j);
  double delta = param_list["a"];
  int ndim = x.n_cols; 
  int nstate = 2;
  mat ss(ndim, n, fill::zeros); 
  vec accs(n, fill::zeros); 
  vec s = s0; 
  mat ngrids(ndim, 2, fill::zeros); 
  ngrids.col(1).fill(1.0); 

  List old_quantity = log_quantity(x, z, s, api, bpi, asigma, bsigma, g, l);
  double oldmass = old_quantity["log_mass"];
  vec oldgrad = as<vec>(old_quantity["grads"]);
  mat olddelta =  pow(ngrids - repmat(s, 1, nstate), 2);
  // Main MCMC loop
  for (int i = 0; i < n; i++) {

    mat energys = 0.5 * (ngrids % repmat(oldgrad, 1, nstate)) - (1.0 / (2 * delta)) * olddelta;
    mat ps = exp(energys);
    ps.each_col() /= sum(ps, 1); 
    vec snew(ndim, fill::zeros); 
    vec w = randu<vec>(ndim); 
    snew = conv_to<vec>::from(w > ps.col(0));  
    List new_quantity = log_quantity(x, z, snew, api, bpi, asigma, bsigma, g, l);
    double newmass = new_quantity["log_mass"];
    vec newgrad = as<vec>(new_quantity["grads"]);

    mat newdelta = pow(ngrids-repmat(snew, 1, nstate), 2);
    double qdown = 1.0;
    double qup = 1.0;

    vec down_exp = exp(0.5 * oldgrad % snew - (1.0 / (2 * delta)) * pow(s - snew, 2));
    vec down_sums = sum(exp(0.5 * ngrids % repmat(oldgrad, 1, nstate) - (1.0 / (2 * delta)) * olddelta), 1);  
    vec up_exp = exp(0.5 * newgrad % s - (1.0 / (2 * delta)) * pow(s - snew, 2));
    vec up_sums = sum(exp(0.5 * ngrids %repmat(newgrad, 1, nstate) - (1.0 / (2 * delta)) * newdelta), 1);  
    double ratio = exp(newmass - oldmass);
    double Q = prod(up_exp%down_sums/(down_exp%up_sums)) * ratio;

    vec q = randu<vec>(1); 
    if (q[0] <= Q) {
      s = snew;
      accs[i] = 1;
      oldgrad = newgrad;
      oldmass = newmass;
      olddelta = newdelta;
    } else {
      accs[i] = 0;
    }

    ss.col(i) = s;
  }

  return List::create(Named("ss") = ss, Named("accs") = accs);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List AVG(int n, List param_list, int j, mat& x,  const vec& z, const vec &s0, double api, double bpi, double asigma, double bsigma, double g, double l) {
  
  arma_rng::set_seed(j);
  double delta = param_list["a"];
  int ndim = x.n_cols; 
  int nstate = 2;
  mat ss(ndim, n, fill::zeros); 
  vec accs(n, fill::zeros); 
  vec s = s0; 
  mat ngrids(ndim, 2, fill::zeros); 
  ngrids.col(1).fill(1.0); 

  List old_quantity = log_quantity(x, z, s, api, bpi, asigma, bsigma, g, l);
  double oldmass = old_quantity["log_mass"];
  vec oldgrad = as<vec>(old_quantity["grads"]);

  // Main MCMC loop
  for (int i = 0; i < n; i++) {
    vec y = sqrt(2/delta) * s + randn<vec>(ndim);
    mat deltamat = pow(sqrt(2/delta)*ngrids - repmat(y, 1, nstate), 2);
    mat energys = (ngrids % repmat(oldgrad, 1, nstate)) - 0.5*deltamat;
    mat ps = exp(energys);
    ps.each_col() /= sum(ps, 1); 
    vec snew(ndim, fill::zeros); 
    vec w = randu<vec>(ndim); 
    snew = conv_to<vec>::from(w > ps.col(0));  
    List new_quantity = log_quantity(x, z, snew, api, bpi, asigma, bsigma, g, l);
    double newmass = new_quantity["log_mass"];
    vec newgrad = as<vec>(new_quantity["grads"]);

    double qdown = 1.0;
    double qup = 1.0;

    vec down_exp = exp(oldgrad % snew - 0.5 * pow(y - sqrt(2/delta)*snew, 2));
    vec down_sums = sum(exp(ngrids % repmat(oldgrad, 1, nstate) - 0.5 * deltamat), 1);  
    vec up_exp = exp(newgrad % s - 0.5 * pow(y - sqrt(2/delta)*s, 2));
    vec up_sums =  sum(exp(ngrids % repmat(newgrad, 1, nstate) - 0.5 * deltamat), 1);   
    double ratio = newmass-oldmass-0.5*sum(pow(y-sqrt(2/delta)*snew,2))+0.5*sum(pow(y-sqrt(2/delta)*s,2));
    double Q = prod(up_exp%down_sums/(down_exp%up_sums)) * exp(ratio);

    vec q = randu<vec>(1); 
    if (q[0] <= Q) {
      s = snew;
      accs[i] = 1;
      oldgrad = newgrad;
      oldmass = newmass;
    } else {
      accs[i] = 0;
    }

    ss.col(i) = s;
  }

  return List::create(Named("ss") = ss, Named("accs") = accs);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List DHAMS(int n, List param_list, int j, mat& x,  const vec& z, const vec &s0, double api, double bpi, double asigma, double bsigma, double g, double l){
  
  arma_rng::set_seed(j);

  double delta = param_list["a"];
  double eps = param_list["b"];
  double ga = param_list["c"];
  double phi = param_list["d"];

  int ndim = x.n_cols; 
  int nstate = 2;
  mat ss(ndim, n, fill::zeros); 
  mat us(ndim, n, fill::zeros); 
  vec accs(n, fill::zeros); 
  vec s = s0;
  vec u = randn<vec>(ndim);
  mat ngrids(ndim, 2, fill::zeros); 
  ngrids.col(1).fill(1.0); 
  List old_quantity = log_quantity(x, z, s, api, bpi, asigma, bsigma, g, l);
  double oldmass = old_quantity["log_mass"];
  vec oldgrad = as<vec>(old_quantity["grads"]);

  for(int i = 0; i<n; i++){
    vec u12 = eps*u +sqrt(1-eps*eps)*randn<vec>(ndim);
    vec y = s-delta*u12;
    mat olddelta = pow(ngrids - repmat(y, 1, nstate), 2);
    mat energys = (ngrids % repmat(oldgrad, 1, nstate)) - (0.5/(delta*delta))*olddelta;
    mat ps = exp(energys);
    ps.each_col() /= sum(ps, 1); 
    vec embednoise = randu<vec>(ndim);
    vec oldembed = (ps.col(0))%embednoise%(1-2*s)+(ps.col(0)+embednoise)%s;
    vec temp = -oldembed + ga * randu<vec>(ndim);
    vec snew = conv_to<vec>::from((temp -arma::floor(temp))>ps.col(0));
    List new_quantity = log_quantity(x, z, snew, api, bpi, asigma, bsigma, g, l);
    double newmass = new_quantity["log_mass"];
    vec newgrad = as<vec>(new_quantity["grads"]);
    vec unew = (1/delta)*(y-snew)+phi*(newgrad-oldgrad);
    mat newdelta = pow(ngrids-repmat(snew+delta*unew, 1, nstate), 2);
    energys = (ngrids%repmat(newgrad, 1, nstate)) -(0.5/(delta*delta))*newdelta;
    mat pback = exp(energys);
    pback.each_col() /= sum(pback, 1);
    double qprod = 1.0;
    vec cump(2, fill::zeros);
    vec cumpback(2, fill::zeros);
    for(int j=0; j< ndim; j++){
      cump = cumsum(ps.row(j).t());
      cumpback = cumsum(pback.row(j).t());
      qprod = qprod*trans_prob(snew[j]+1, s[j]+1, cumpback, -1.0, ga)/trans_prob(s[j]+1, snew[j]+1, cump, -1.0, ga);
    }
    double ratio = newmass - oldmass - 0.5 * (sum(unew%unew) - sum(u12%u12));
    double Q = exp(ratio)*qprod;
    vec q = randu<vec>(1); 
    if (q[0] <= Q) {
      s = snew;
      u = unew;
      accs[i] = 1;
      oldgrad = newgrad;
      oldmass = newmass;
    } else {
      u = -u12;
      accs[i] = 0;
    }
    ss.col(i) = s;
    us.col(i) = u;
  }
  return List::create(Named("ss") = ss, Named("us") = us, Named("accs") = accs);
}