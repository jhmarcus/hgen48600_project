#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Eigen;
using namespace Rcpp;

// [[Rcpp::export]]
double normTransProb(double x0, double xt, int nGen, double s, double h, int N) {
  double var = x0 * (1 - x0);
  //if (var <= 0.0){
  //  var = .00001;
  //}
  double mu = x0 + (2 * s) * var * (x0 + (h * (1 - (2 * x0))));
  double sigma = sqrt( var * ((double)nGen / (2 * N)) ); 
  float eps = .001;
  float transProb = R::pnorm(xt + eps, mu, sigma, 1, 0) - R::pnorm(xt - eps, mu, sigma, 1, 0);
  return transProb; 
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
MatrixXd foward(MatrixXd O, VectorXd states, int nStates, double s, double h, int N) {
  int i,j,k;
  int alleleCount;
  int nChr;
  int nGen;
  VectorXd a(nStates);
  double emissProb;

  int nObs = O.rows();
  MatrixXd F(nStates, nObs);
  
  // intialization 
  for (i=0; i<nStates; i++){
    emissProb = R::dbinom(O(0, 0), O(0, 1), states[i], 0);
    F(i, 0) = (1.0 / nStates) * emissProb;
  }
  
  // induction 
  for (i=1; i<nObs; i++){
    for (j=0; j<nStates; j++){
      alleleCount = O(i, 0);
      nChr = O(i, 1);
      nGen = O(i, 2) - O(i-1, 2);
      for (k=0; k<nStates; k++){
        a[k] = normTransProb(states[k], states[j], nGen, s, h, N);
      }
      emissProb = R::dbinom(alleleCount, nChr, states[j], 0);
      F(j, i) = (F.col(i-1).dot(a)) * emissProb;
    }
  }
  return F;
}

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