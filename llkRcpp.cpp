#include <Rcpp.h>
using namespace Rcpp;
// #include <cmath>
// #include <math.h>

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
double llkRcpp(List data, List par) {
  // ===========================================================================
  // Extract data objects
  // ---------------------------------------------------------------------------
  const NumericMatrix W = data["W"];            // detection histories
  const NumericMatrix R = data["R"];            // received levels
  const NumericMatrix Y_rec = data["Y_rec"];    // recorded bearings
  const NumericMatrix Y_grid = data["Y_grid"];       // grid bearings
  const NumericMatrix X = data["X"];            // distances
  const NumericVector A = data["A"];            // areas of grid cells of mesh
  const int trunc_level = data["trunc_level"];  // truncation level for received levels
  const NumericMatrix design_matrix = 
    data["design_matrix"];                // the design matrix from the GAM object
  const NumericVector S = data["S"];            // source level integration
  // ===========================================================================
  // ===========================================================================
  // Extract individual parameters
  // ---------------------------------------------------------------------------
  const double g0 = exp(double(par["logit_g0"])) / 
    (exp(double(par["logit_g0"])) + 1.0);
  const double beta_r = exp(double(par["log_beta_r"]));
  const double sd_r = exp(double(par["log_sd_r"]));
  const double mu_s = exp(double(par["log_mu_s"]));
  const double sd_s = exp(double(par["log_sd_s"]));
  
  const double kappa_low = exp(double(par["log_kappa_low"]));
  const double kappa_high = exp(double(par["log_kappa_high"])) + kappa_low;
  const double mix_bear = exp(double(par["logit_mix_bear"])) / 
    (exp(double(par["logit_mix_bear"])) + 1.0);
  
  // std::cout << "mix_bear: " << mix_bear << std::endl;
  // std::cout << "R::bessel_i(5.0, 0.0, 1.0): " << R::bessel_i(5.0, 0.0, 1.0) << std::endl;
  
  const double log_A_s = log(5.0);
  
  const NumericVector dens_pars = par["density_pars"];
  
  const double bessel_low = R::bessel_i(kappa_low, 0.0, 1.0);
  const double bessel_high = R::bessel_i(kappa_high, 0.0, 1.0);
  // ===========================================================================
  
  // ===========================================================================
  // Extract relevant R functions [NOT NEEDED RIGHT NOW]
  // ---------------------------------------------------------------------------
  // // Obtain environment containing function
  // Environment matrixStats = Environment::namespace_env("matrixStats");
  // // Environment matrixStats("package:matrixStats"); 
  // // Make function callable from C++
  // Function logSumExp = matrixStats["logSumExp"];    
  // ===========================================================================
  
  // Extract dimensions
  const int n_call = W.rows();
  const int n_det = W.cols();
  const int n_grid = X.rows();
  const int n_dens_pars = dens_pars.size();
  const int n_sl = S.size();
  
  // Predict the density, D, a vector of length = n_grid
  NumericVector D(n_grid);
  for (int m = 0; m < n_grid; m++) {
    double predictor = 0;
    for (int index = 0; index < n_dens_pars; index++) {
      predictor += dens_pars[index] * design_matrix(m, index);
    }
    D(m) = exp(predictor); // use a log link to have strictly positive density
  }
  
  // Get general probabilities
  NumericVector log_probs_s = dnorm(S, mu_s, sd_s, true);
  NumericVector probs_s = exp(log_probs_s);
  
  // std::cout << "Number of observations: " << n_call << std::endl;
  // std::cout << "Number of detectors: " << n_det << std::endl;
  // std::cout << "Grid size of the mesh: " << n_grid << std::endl;
  // std::cout << "Number of source level values tried: " << n_sl << std::endl;
  // std::cout << "kappa_high parameter: " << kappa_high << std::endl;
  
  // std::cout<< "1 + INFINITY: " << 1 + INFINITY << std::endl; // Is still inf!
  
  // std::cout << "test: " << - pow(10, 16) << std::endl;
  
  // std::cout << "X(0, 0): " << X(0, 0) << std::endl;
  
  // ===========================================================================
  // THE NON-CONDITIONAL PART OF THE LOG LIKELIHOOD
  // ===========================================================================
  
  // Define the likelihood of n_calls with constant density
  NumericVector lambda_vector(n_sl);
  
  for (int sl_index = 0; sl_index < n_sl; sl_index++) {
    NumericVector lambda_vector_sl(n_grid);
    
    for (int m = 0; m < n_grid; m++) { // loop through the mesh
      NumericVector probs_m(n_det);
      NumericVector non_probs_m(n_det);
      for (int j = 0; j < n_det; j++) { // loop through the detectors
        // std::cout << j << std::endl; // works
        
        double expected_rl = S[sl_index] - beta_r * log10(X(m, j));
        // std::cout << expected_rl << std::endl; // works
        
        double xx = (trunc_level - expected_rl) / sd_r;
        double p = g0 * (1.0 - R::pnorm(xx, 0.0, 1.0, true, false));
        probs_m[j] = p;
        // std::cout << p << std::endl; // works
        
        non_probs_m[j] = 1.0 - probs_m[j];
      }
      
      // Get probability of at least two detections
      NumericVector v1(n_det + 1);
      for (int a = 0; a < n_det; a++) { // derive all probabilities of one detection
        // Vector for configurations of a single detection
        NumericVector config_vector(n_det);
        // Vector to store the probabilities associated with v2 configurations
        NumericVector probability_vector(n_det); 
        for (int iii = 0; iii < n_det; iii++) {
          config_vector[iii] = 0; // fill the vector v2 with zeros
        }
        config_vector(a) = 1; // set one position equal to
        // std::cout << "config_vector: " << config_vector << std::endl;
        for (int ii = 0; ii < n_det; ii++) {
          if (config_vector[ii] == 1) probability_vector[ii] = probs_m[ii];
          if (config_vector[ii] == 0) probability_vector[ii] = non_probs_m[ii];
        }
        // std::cout << "probability_vector: " << probability_vector << std::endl;
        // std::cout << "probs_m: " << probs_m << std::endl;
        // std::cout << "non_probs_m: " << non_probs_m << std::endl;
        NumericVector temp = cumprod(probability_vector);
        v1(a) = rev(temp)[0];
      }
      NumericVector temp = cumprod(non_probs_m);
      v1(n_det) = rev(temp)[0]; // probability of no detections
      double p_twice = 1.0 - sum(v1); // p_twice is complement of p_none and p_once
      // std::cout << "p_twice: " << p_twice << " and v1: " << v1 << std::endl;
      
      // store lambda given x in lambda vector
      lambda_vector_sl[m] = D[m] * p_twice * A[m];
    }
    lambda_vector[sl_index] = sum(lambda_vector_sl) * exp(log_A_s) * 
      probs_s[sl_index];
  }
  
  double lambda = sum(lambda_vector);
  
  // std::cout << "The value for lambda is: " << lambda << std::endl;
  
  // Type llk_n = dpois(Type(n_call), lambda, true); 
  double llk_n = -lambda; // Only need part due to logarithm in llk
  
  if (llk_n  == -INFINITY) llk_n = -pow(10, 16);
  if (llk_n  == INFINITY) llk_n = pow(10, 16);
  
  // ===========================================================================
  // THE CONDITIONAL PART OF THE LIKELIHOOD
  // ===========================================================================
  
  // Define the conditional negative log likelihood
  double llk_cond = 0.0; // define vector to store likelihoods per individual call
  
  for (int i = 0; i < n_call; i++) { // loop over the calls
    // std::cout << "call number: " << i << std::endl; 
    NumericVector log_probs_i(n_sl);
    for (int sl_index = 0; sl_index < n_sl; sl_index++) {
      double log_prob_s = log_probs_s[sl_index];
      NumericVector log_probs_i_sl(n_grid);
      for (int m = 0; m < n_grid; m++) { // loop through the mesh
        // =====================================================================
        // Likelihood of capture histories
        // ---------------------------------------------------------------------
        double log_prob_w_im = 0.0;
        
        for (int j = 0; j < n_det; j++) {
          double expected_rl = S[sl_index] - beta_r * log10(X(m, j));
          double xx = (trunc_level - expected_rl) / sd_r;
          double prob = g0 * (1.0 - R::pnorm(xx, 0.0, 1.0, true, false));
          
          int test = 1 * W(i, j); // Multiply W(i, j)  by one so to not use DATA directly
          
          if (test == 1) log_prob_w_im += log(prob); 
          else log_prob_w_im += log(1.0 - prob);
        }
        if (log_prob_w_im == -INFINITY) log_prob_w_im = -pow(10, 16);
        if (log_prob_w_im == INFINITY) log_prob_w_im = pow(10, 16);
        // =====================================================================
        
        // =====================================================================
        // Likelihood of bearings
        // ---------------------------------------------------------------------
        double log_prob_y_im = 0.0;
        
        for (int j = 0; j < n_det; j++) {
          if (W(i, j) == 1) {
            double cos_obs_minus_exp = cos(Y_rec(i, j) - Y_grid(m, j));
            if (false) { // REPLACE LATER FOR MIXTURE CHECK
              // log_prob_y_im += kappa * cos(obs_minus_exp) - 
              //   log(2 * M_PI * besselI(kappa, Type(0)));// log version
            } else {
              log_prob_y_im += log(mix_bear * exp(kappa_low * cos_obs_minus_exp) / 
                (2.0 * M_PI * bessel_low) +
                (1.0 - mix_bear) * exp(kappa_high * cos_obs_minus_exp) / 
                (2.0 * M_PI * bessel_high));
            }
          }
        }
        if (log_prob_y_im == -INFINITY) log_prob_y_im = -pow(10, 16);
        if (log_prob_y_im == INFINITY) log_prob_y_im = pow(10, 16);
        // =====================================================================
        
        // =====================================================================
        // Likelihood of received levels
        // ---------------------------------------------------------------------
        double log_prob_r_im = 0.0; // sum the log probabilities
        
        // Loop through all detectors and sum the log probabilities
        for (int j = 0; j < n_det; j++) {
          if (W(i, j) == 1) {
            double expected_rl = S[sl_index] - beta_r * log10(X(m, j));
            double recorded_rl = 1.0 * R(i, j);
            
            // Probability of exp_rl - rec_rl
            double log_prob = R::dnorm(recorded_rl, expected_rl, sd_r, true); // log probabilities
            
            log_prob -= log(sd_r * (1.0 - 
              R::pnorm(trunc_level, expected_rl, sd_r, true, false))); // Different from R code,
                                                                       // but maths is the same.
            
            log_prob_r_im += log_prob;
          }
        }
        // Set log_prob_r_im to -pow(10, 16) if it is -inf
        if (log_prob_r_im == -INFINITY) log_prob_r_im = -pow(10, 16);
        if (log_prob_r_im == INFINITY) log_prob_r_im = pow(10, 15); // 16 - 1 = 15, just to create a difference
        // =====================================================================
        
        if (log_prob_s == -INFINITY) log_prob_s = -pow(10, 16);
        
        // std::cout << "log(D): " << log(D) << std::endl; // works
        // std::cout << "prob_w_im: " << prob_w_im << std::endl; // works for now
        // std::cout << "prob_y_im: " << prob_y_im << std::endl; // works
        // std::cout << "prob_r_im: " << prob_r_im << std::endl; // works
        // log likelihood of a single call given a single location
        double cond_llk_i_sl_m = 
          log(A[m]) + 
          log(D[m]) + 
          log_prob_w_im + 
          log_prob_y_im + 
          log_prob_r_im + 
          log_prob_s; 
        
        bool test = std::isnan(cond_llk_i_sl_m); // isnan from cmath requires a double
        if (test)  {
          std::cout << "log_prob_w_im: " << log_prob_w_im <<
            "\nlog_prob_y_im: " << log_prob_y_im <<
              "\nlog_prob_r_im: " << log_prob_r_im << 
                "\nlog_prob_s: " << R::dnorm(S[sl_index], mu_s, sd_s, true) << std::endl;
        }
        
        log_probs_i_sl[m] = cond_llk_i_sl_m;
      }
      // Sum the values in log_probs_i_sl for every grid point, in a similar way to logsumexp()
      // Can we use the function logspace_add() ?
      double log_sum_i_sl = log_probs_i_sl[0];
      for (int ii = 1; ii < n_grid; ii++) {
        log_sum_i_sl = R::logspace_add(log_sum_i_sl, log_probs_i_sl[ii]);
      }
      log_probs_i[sl_index] = log_sum_i_sl + log_A_s;
      // // REPLACE THE ABOVE WITH 
      // NumericVector log_sum_i_sl = logSumExp(log_probs_i_sl); // output as NumericVector, as this works
      // log_probs_i[sl_index] = log_sum_i_sl[0] + log_A_s;
    }
    
    // Sum the values in log_probs_i for every source level, in a similar way to logsumexp()
    double log_sum_i = log_probs_i[0];
    for (int ii = 1; ii < n_sl; ii++) {
      log_sum_i = R::logspace_add(log_sum_i, log_probs_i[ii]);
    }
    llk_cond += log_sum_i;
    // // REPLACE THE ABOVE WITH
    // NumericVector log_sum_i = logSumExp(log_probs_i); // output as NumericVector, as this works
    // llk_cond += log_sum_i[0];
    // std::cout << "llk_cond after " << i << " additions: " << llk_cond << std::endl;
  }
  // Combine the two log likelihood elements and make it negative
  double nll = -(llk_n + llk_cond);

  std::cout << "The total n part of the log likelihood: " << llk_n << std::endl; // nan
  std::cout << "The conditional part of the log likelihood: " << llk_cond << std::endl; // nan
  std::cout << "The negative log likelihood: " << nll << std::endl;
  
  // Return the nll
  return nll;
}
