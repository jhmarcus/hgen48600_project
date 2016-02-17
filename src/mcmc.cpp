#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Eigen;
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]

//' Normal approximation to Wright-Fisher with selection transition probability 
//' @param x_0 ancestral allele frequency
//' @param x_t current allele frequency 
//' @param n_gen number of generations between ancestral and current populations
//' @param s selection coefficent
//' @param h dominence parameter
//' @param N effective p
//' @return probability of transitioning from state x_0 to state x_t
// [[Rcpp::export]]
double norm_trans_prob(double x_t, double x_0, int n_gen, double s, double h, int N) {
  double var = x_0 * (1 - x_0); // helper for computing mu and sigma
  double mu = x_0 + (2 * s) * var * (x_0 + (h * (1 - (2 * x_0))));
  double sigma = sqrt( var * ((double)n_gen / (2 * N)) ); 
  float eps = .001;
  float trans_prob = R::pnorm(x_t + eps, mu, sigma, 1, 0) - R::pnorm(x_t - eps, mu, sigma, 1, 0);
  return trans_prob; 
}

// [[Rcpp::export]]
MatrixXd foward(MatrixXd O, VectorXd states, int n_states, double s, double h, int N) {
  int i,j,k;
  int allele_count;
  int n_chr;
  int n_gen;
  VectorXd a(n_states);
  double emiss_prob;

  int n_obs = O.rows();
  MatrixXd F(n_states, n_obs);
  
  // intialization 
  for (i=0; i<n_states; i++){
    emiss_prob = R::dbinom(O(0, 0), O(0, 1), states[i], 0);
    F(i, 0) = (1.0 / n_states) * emiss_prob;
  }
  
  // induction 
  for (i=1; i<n_obs; i++){
    for (j=0; j<n_states; j++){
      allele_count = O(i, 0);
      n_chr = O(i, 1);
      n_gen = O(i, 2) - O(i-1, 2);
      for (k=0; k<n_states; k++){
        a[k] = norm_trans_prob(states[k], states[j], n_gen, s, h, N);
      }
      emiss_prob = R::dbinom(allele_count, n_chr, states[j], 0);
      F(j, i) = (F.col(i-1).dot(a)) * emiss_prob;
    }
  }
  return F;
}

/*
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
VectorXd mcmc(MatrixXd O, VectorXd states, int nStates, double s0, double h, int N, double propSd, int nIter){
  VectorXd posteriorSamples(nIter);
  int i;
  double currentS, newS, currentPrior, newPrior, currentLikelihood, newLikelihood, A;
  int nObs = O.rows();
  posteriorSamples[0] = s0;
  for (i=1; i<nIter; i++){
    currentS = posteriorSamples[i-1];
    newS = currentS + R::rnorm(0, propSd);
    currentPrior = R::dunif(currentS, -.3, .3, 0);
    
    currentLikelihood = foward(O, states, nStates, currentS, h, N).col(nObs).sum();
    newPrior = R::dunif(newS, -.3, .3, 0);
    
    newLikelihood = foward(O, states, nStates, newS, h, N).col(nObs).sum();
    A = (newPrior * newLikelihood) / (currentPrior * currentLikelihood) ;
    if (R::runif(0, 1) < A){
      posteriorSamples[i] = newS;
    }
    else {
      posteriorSamples[i] = currentS;
    }
  }
  return posteriorSamples;
}
*/