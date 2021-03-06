#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

//' Get Wright-Fisher binomial success probability with selection
//' @param x frequency in the previous generation
//' @param s selection coefficent
//' @param h dominence parameter 
//' @return g(x)
// [[Rcpp::export]]
double g(double x, double s, double h){
  double a = s * x * (1 - x);
  double b = h + ((1 - 2 * h) * x);
  return(x + (a * b));
}

//' Second derivative of g for taylor series approx to moments
//' @param x frequency in the previous generation
//' @param s selection coefficent
//' @param h dominence parameter 
//' @return g''(x)
// [[Rcpp::export]]
double g_2(double x, double s, double h){
  double a = h * (6 * x - 3);
  double b = 3 * x + 1;
  return(2 * s * (a - b));
}


//' Get first two moments for the normal approximation using the delta method
//' @param x_0 ancestral allele frequency
//' @param n_gen number of generations between ancestral and current populations
//' @param s selection coefficent
//' @param h dominence parameter
//' @param N effective population size
// [[Rcpp::export]]
arma::vec moments_delta_method(double x_0, int n_gen, double s, double h, int N){
  int i;
  double a, b;
  arma::vec moments(2);
  double mu = x_0;
  double sigma2 = (1.0 / (2*N)) * x_0 * (1 - x_0);
  for (i=0; i<n_gen; i++){
    a = (1 / ( 2 * N)) * g(mu, s, h) * (1 - g(mu, s, h));
    b = pow(1 + s * (2 - (3 * mu)) * mu, 2) * sigma2;
    sigma2 = a + b;
    mu = g(mu, s, h);
    //mu = g(mu, s, h) + .5 * (1 / ( 2 * N)) * (mu * (1 - mu)) * g_2(mu, s, h);
  }
  moments[0] = mu;
  moments[1] = sigma2;
  return moments;
}

//' Normal approximation using the delta method approximation to the moments
//' @param x_0 ancestral allele frequency
//' @param x_t current allele frequency 
//' @param n_gen number of generations between ancestral and current populations
//' @param s selection coefficent
//' @param h dominence parameter
//' @param N effective population size
//' @return probability of transitioning from state x_0 to state x_t
// [[Rcpp::export]]
double norm_delta_trans_prob(double x_t, double x_0, int n_gen, double s, double h, int N){
    arma::vec moments = moments_delta_method(x_0, n_gen, s, h, N);
    double mu = moments[0];
    double sigma2 = moments[1];
    double eps = .0001;
    double trans_prob = R::pnorm(x_t + eps, mu, sqrt(sigma2), 1, 0) - R::pnorm(x_t - eps, mu, sqrt(sigma2), 1, 0);
    return trans_prob; 
}

//' Normal approximation from Lacerda and Seoighe 2014 under weak selection and mutation
//' @param x_0 ancestral allele frequency
//' @param x_t current allele frequency 
//' @param n_gen number of generations between ancestral and current populations
//' @param s selection coefficent
//' @param h dominence parameter
//' @param N effective population size
//' @return probability of transitioning from state x_0 to state x_t
// [[Rcpp::export]]
double norm_weak_s_trans_prob(double x_t, double x_0, int n_gen, double s, double h, int N) {
  double var = x_0 * (1 - x_0); // helper variable
  double exp_ns = exp(n_gen * s); // helper variable
  double mu = x_0 / (x_0 + ((1 - x_0) * exp(-s * n_gen)));
  double sigma2_num_a = var * exp_ns;
  double sigma2_num_b = pow(x_0, 2) * exp(2 * n_gen * s);
  double sigma2_num_c = (1 - (2 * x_0) - (2 * n_gen * s * x_0) * (x_0 - 1)) * exp_ns;
  double sigma2_num_d = pow(x_0 - 1, 2);
  double sigma2_denom = N * s * pow(1 + x_0 * (exp_ns - 1), 4);
  double sigma2_num = sigma2_num_a * (sigma2_num_b + (sigma2_num_c - sigma2_num_d));
  double sigma2 = sigma2_num / sigma2_denom;
  double eps = .0001;
  double trans_prob = R::pnorm(x_t + eps, mu, sqrt(sigma2), 1, 0) - R::pnorm(x_t - eps, mu, sqrt(sigma2), 1, 0);
  return trans_prob; 
}

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

//' Foward algorithim of the HMM using the WF normal approximation transition probabilites
//' @param O observation matrix first column is the allele_count, the second column is n_chromsomes, 
//' the third column is the generation
//' @param states is a vector of the discretized state space in this case discretized allele frequencies
//' @param s the selection coefficent
//' @param h the dominence parameter
//' @param N the effective population size
//' @return a n_states x n_obs matrix F that stores the forward variables 
// [[Rcpp::export]]
arma::mat foward(arma::mat O, arma::vec states, double s, double h, int N) {
  int i,j,k;
  int allele_count;
  int n_chr;
  int n_gen;
  int n_states = states.size();
  arma::vec a(n_states);
  double emiss_prob;

  int n_obs = O.n_rows; 
  arma::mat F(n_states, n_obs);
  F.zeros();
  //F.print("F = ");
  
  // intialization 
  for (i=0; i<n_states; i++){
    allele_count = O(0, 0);
    n_chr = O(0, 1);
    emiss_prob = R::dbinom(allele_count, n_chr, states[i], 0);
    //printf("%e \n", emiss_prob);
    F(i, 0) = (1.0 / n_states) * emiss_prob;
  }
  //F.print("F = ");
  
  // induction 
  for (i=1; i<n_obs; i++){
    for (j=0; j<n_states; j++){
      allele_count = O(i, 0);
      n_chr = O(i, 1);
      n_gen = O(i, 2) - O(i-1, 2);
      //printf("%i, %i, %i\n", allele_count, n_chr, n_gen);
      for (k=0; k<n_states; k++){
        //a[k] = norm_trans_prob(states[j], states[k], n_gen, s, h, N);
        //a[k] = norm_weak_s_trans_prob(states[j], states[k], n_gen, s, h, N);
        a[k] = norm_delta_trans_prob(states[j], states[k], n_gen, s, h, N);
        
      }
      //a.t().print("a = ");
      emiss_prob = R::dbinom(allele_count, n_chr, states[j], 0);
      F(j, i) = sum(F.col(i-1) % a) * emiss_prob;
    }
  }
  //F.print("F =");
  //printf("%e\n", sum(F.col(n_obs-1)));
  return F;
}

//' @param O observation matrix first column is the allele_count, the second column is n_chromsomes, 
//' the third column is the generation
//' @param states is a vector of the discretized state space in this case discretized allele frequencies
//' @param s the selection coefficent
//' @param h the dominence parameter
//' @param N the effective population size
//' @param prop_sd standerd deviation of the proposal distribution 
//' @param n_iter number of iterations of mcmc
//' @return a vector of samples from the posterior distribution of s
// [[Rcpp::export]]
arma::vec mcmc(arma::mat O, arma::vec states, double s_0, double h, int N, double prop_sd, int n_iter){
  arma::vec posterior_samples(n_iter);
  int i;
  double current_s, new_s, current_prior, new_prior, current_likelihood, new_likelihood, A;
  int n_obs = O.n_rows;
  int n_states = states.size();
  posterior_samples[0] = s_0;
  for (i=1; i<n_iter; i++){
    current_s = posterior_samples[i-1];
    new_s = current_s + R::rnorm(0, prop_sd);
    
    current_prior = R::dnorm(current_s, 0, .1, 0);
    current_likelihood = sum(foward(O, states, current_s, h, N).col(n_obs-1));
    
    new_prior = R::dnorm(new_s, 0, .1, 0);
    new_likelihood = sum(foward(O, states, new_s, h, N).col(n_obs-1));
    
    A = (new_prior * new_likelihood) / (current_prior * current_likelihood);
    if (R::runif(0, 1) < A){
      posterior_samples[i] = new_s;
    }
    else {
      posterior_samples[i] = current_s;
    }
    
    //if (i % 100 == 0){
      //printf("mcmc_iteration: %i \n", i);
    //}
  }
  return posterior_samples;
}