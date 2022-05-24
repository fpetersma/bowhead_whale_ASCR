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

// =============================================================================
// VARIABLE SOURCE LEVEL SCRIPTS
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
double uncondNllRcpp(List dat, List par) {
  // ===========================================================================
  // 1. EXTRACT DATA OBJECTS
  // ---------------------------------------------------------------------------
  const NumericMatrix X = dat["X"];            // distances
  const NumericVector A_x = dat["A_x"];        // areas of grid cells of mesh
  const int trunc_level = dat["trunc_level"];  // truncation level for received levels
  const NumericVector S = dat["S"];            // source level integration
  const NumericVector S_probs = dat["S_probs"];// source level probabilities
  const double A_s = dat["A_s"];               // spacing of S integration
  const NumericVector D = dat["D"];            // density for every grid point
  
  // Extract dimensions
  const int n_det = X.cols();
  const int n_grid = X.rows();
  const int n_sl = S.size();
  // ===========================================================================
  
  // ===========================================================================
  // 2. EXTRACT INDIVIDUAL PARAMETERS
  // ---------------------------------------------------------------------------
  const double g0 = exp(double(par["logit_g0"])) / 
    (exp(double(par["logit_g0"])) + 1.0);
  const double beta_r = exp(double(par["log_beta_r"]));
  const double sd_r = exp(double(par["log_sd_r"]));
  // const double mu_s = exp(double(par["log_mu_s"]));
  // const double sd_s = exp(double(par["log_sd_s"]));
  // ===========================================================================
  
  // ===========================================================================
  // 3. DERIVE LAMBDA, THE EXPECTED NUMBER OF DETECTED CALLS
  // ---------------------------------------------------------------------------
  // Define the likelihood of n_calls with constant density
  double lambda = 0.0;
  
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
      lambda_vector_sl[m] = D[m] * p_twice * A_x[m];
    }
    lambda += sum(lambda_vector_sl) * A_s * S_probs[sl_index];
  }
  
  // std::cout << "The value for lambda is: " << lambda << std::endl;
  
  
  // ===========================================================================
  // 4. RETURN THE UNCOND. NEGATIVE LOG LIKELIHOOD
  // ---------------------------------------------------------------------------
  double nll_n = lambda; // Only need part due to logarithm in nll
  
  // Check if nll_cond is (-)Inf
  if (nll_n  == -INFINITY) nll_n = -pow(10, 16);
  if (nll_n  == INFINITY) nll_n = pow(10, 16);
  
  // Return the nll_cond
  return nll_n;
  // ===========================================================================
}

// [[Rcpp::export]]
double singleNllRcpp(List dat, List par) {
  // ===========================================================================
  // 1. EXTRACT DATA OBJECTS
  // ---------------------------------------------------------------------------
  const NumericVector W = dat["W"];            // detection histories
  const NumericVector R = dat["R"];            // received levels
  const NumericVector Y_rec = dat["Y_rec"];    // recorded bearings
  const NumericMatrix Y_grid = dat["Y_grid"];  // grid bearings
  const NumericMatrix X = dat["X"];            // distances
  const NumericVector A_x = dat["A_x"];        // areas of grid cells of mesh
  const int trunc_level = dat["trunc_level"];  // truncation level for received levels
  const NumericVector S = dat["S"];            // source level integration
  const NumericVector S_probs_log = 
    dat["S_probs_log"];                        // source level probabilities
  const double A_s = dat["A_s"];               // spacing of S integration
  const double log_A_s = log(A_s);
  const NumericVector D = dat["D"];            // density for every grid point
  
  const int USE_BEARINGS = dat["USE_BEARINGS"];// 0: no bearings, 1: single
                                               // kappa; 2: mix kappa
  // Declare all to allow compiling
  double bessel, bessel_low, bessel_high; 
  
  if (USE_BEARINGS == 1) {
    bessel = dat["bessel"]; // the BesselI() values
  }
  if (USE_BEARINGS == 2) {
    bessel_low = dat["bessel_low"]; // the BesselI() values
    bessel_high = dat["bessel_high"];

  }

  // ===========================================================================
  
  // ===========================================================================
  // 2. EXTRACT INDIVIDUAL PARAMETERS
  // ---------------------------------------------------------------------------
  const double g0 = exp(double(par["logit_g0"])) / 
    (exp(double(par["logit_g0"])) + 1.0);
  const double beta_r = exp(double(par["log_beta_r"]));
  const double sd_r = exp(double(par["log_sd_r"]));
  const double mu_s = exp(double(par["log_mu_s"]));
  const double sd_s = exp(double(par["log_sd_s"]));
  
  // Declare all at 0.0 to enable compiling
  double kappa, kappa_low, kappa_high, mix_bear;
  
  if (USE_BEARINGS == 1) {
    kappa = exp(double(par["log_kappa"]));
  }
  if (USE_BEARINGS == 2) {
    kappa_low = exp(double(par["log_kappa_low"]));
    kappa_high = exp(double(par["log_kappa_high"])) + kappa_low;
    mix_bear = exp(double(par["logit_mix_bear"])) / 
      (exp(double(par["logit_mix_bear"])) + 1.0);
  }
  // const NumericVector dens_pars = par["density_pars"];
  
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
  const int n_det =  dat["n_det"];
  const int n_grid =  dat["n_grid"];
  const int n_sl =  dat["n_sl"];
  
  // std::cout << "Number of observations: " << n_call << std::endl;
  // std::cout << "Number of detectors: " << n_det << std::endl;
  // std::cout << "Grid size of the mesh: " << n_grid << std::endl;
  // std::cout << "Number of source level values tried: " << n_sl << std::endl;
  // std::cout << "kappa_high parameter: " << kappa_high << std::endl;
  
  // std::cout<< "1 + INFINITY: " << 1 + INFINITY << std::endl; // Is still inf!
  
  // std::cout << "test: " << - pow(10, 16) << std::endl;
  
  // std::cout << "X(0, 0): " << X(0, 0) << std::endl;
  
  // ===========================================================================
  // THE CONDITIONAL PART OF THE LIKELIHOOD
  // ===========================================================================
  
  // // Define the conditional negative log likelihood
  // double llk_cond = 0.0; // define vector to store likelihoods per individual call
  
  NumericVector log_probs_i(n_sl);
  for (int sl_index = 0; sl_index < n_sl; sl_index++) {
    double log_prob_s = S_probs_log[sl_index];
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
        
        int test = 1 * W[j]; // Multiply W[j]  by one so to not use DATA directly
        
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
      
      if (USE_BEARINGS != 0) {
        for (int j = 0; j < n_det; j++) {
          if (W[j] == 1) {
            double cos_obs_minus_exp = cos(Y_rec[j] - Y_grid(m, j));
            if (USE_BEARINGS == 1) { 
              log_prob_y_im += kappa * cos_obs_minus_exp -
                log(2 * M_PI * bessel);// log version
            } 
            if (USE_BEARINGS == 2) {
              log_prob_y_im += log(mix_bear * exp(kappa_low * cos_obs_minus_exp) / 
                (2.0 * M_PI * bessel_low) +
                (1.0 - mix_bear) * exp(kappa_high * cos_obs_minus_exp) / 
                (2.0 * M_PI * bessel_high));
            }
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
        if (W[j] == 1) {
          double expected_rl = S[sl_index] - beta_r * log10(X(m, j));
          double recorded_rl = 1.0 * R[j];
          
          // Probability of exp_rl - rec_rl
          double log_prob = R::dnorm(recorded_rl, expected_rl, sd_r, true); // log probabilities
          
          log_prob -= log((1.0 - 
            R::pnorm(trunc_level, expected_rl, sd_r, true, false))); // Different from R code,
          // but maths is the same.
          
          log_prob_r_im += log_prob;
        }
      }
      // Set log_prob_r_im to -pow(10, 16) if it is -inf
      if (log_prob_r_im == -INFINITY) log_prob_r_im = -pow(10, 16);
      if (log_prob_r_im == INFINITY) log_prob_r_im = pow(10, 15); // 16 - 1 = 15, just to create a difference
      // =====================================================================
      
      // Check if log_prob_s is -Inf; if so, set to -1e-16
      if (log_prob_s == -INFINITY) log_prob_s = -pow(10, 16);
      
      // std::cout << "log(D): " << log(D) << std::endl; // works
      // std::cout << "prob_w_im: " << prob_w_im << std::endl; // works for now
      // std::cout << "prob_y_im: " << prob_y_im << std::endl; // works
      // std::cout << "prob_r_im: " << prob_r_im << std::endl; // works
      // log likelihood of a single call given a single location
      double cond_llk_i_sl_m = 
        log(A_x[m]) + 
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
  // Make it negative
  double nll_cond = -log_sum_i;
  
  // std::cout << "The conditional part of the log likelihood: " << llk_cond << std::endl; // nan
  // std::cout << "The negative log likelihood: " << nll << std::endl;
  
  // Return the nll
  return nll_cond;
}
// =============================================================================

// =============================================================================
// SINGLE SOURCE LEVEL SCRIPTS
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
double uncondNllRcppFixedSL(List dat, List par) {
  // ===========================================================================
  // 1. EXTRACT DATA OBJECTS
  // ---------------------------------------------------------------------------
  const NumericMatrix X = dat["X"];            // distances
  const NumericVector A_x = dat["A_x"];        // areas of grid cells of mesh
  const int trunc_level = dat["trunc_level"];  // truncation level for received levels
  const NumericVector D = dat["D"];            // density for every grid point
  
  // Extract dimensions
  const int n_det = X.cols();
  const int n_grid = X.rows();
  // ===========================================================================
  
  // ===========================================================================
  // 2. EXTRACT INDIVIDUAL PARAMETERS
  // ---------------------------------------------------------------------------
  const double g0 = exp(double(par["logit_g0"])) / 
    (exp(double(par["logit_g0"])) + 1.0);
  const double beta_r = exp(double(par["log_beta_r"]));
  const double sd_r = exp(double(par["log_sd_r"]));
  const double mu_s = exp(double(par["log_mu_s"]));

  // ===========================================================================
  
  // ===========================================================================
  // 3. DERIVE LAMBDA, THE EXPECTED NUMBER OF DETECTED CALLS
  // ---------------------------------------------------------------------------
  // Define the likelihood of n_calls with constant density
  double lambda = 0.0;
    
  for (int m = 0; m < n_grid; m++) { // loop through the mesh
    NumericVector probs_m(n_det);
    NumericVector non_probs_m(n_det);
    for (int j = 0; j < n_det; j++) { // loop through the detectors
      // std::cout << j << std::endl; // works
      
      double expected_rl = mu_s - beta_r * log10(X(m, j));
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
    lambda += D[m] * p_twice * A_x[m];
  }
  // std::cout << "The value for lambda is: " << lambda << std::endl;
  
  
  // ===========================================================================
  // 4. RETURN THE UNCOND. NEGATIVE LOG LIKELIHOOD
  // ---------------------------------------------------------------------------
  double nll_n = lambda; // Only need part due to logarithm in nll
  
  // Check if nll_cond is (-)Inf
  if (nll_n  == -INFINITY) nll_n = -pow(10, 16);
  if (nll_n  == INFINITY) nll_n = pow(10, 16);
  
  // Return the nll_cond
  return nll_n;
  // ===========================================================================
}

// [[Rcpp::export]]
double singleNllRcppFixedSL(List dat, List par) {
  // ===========================================================================
  // 1. EXTRACT DATA OBJECTS
  // ---------------------------------------------------------------------------
  const NumericVector W = dat["W"];            // detection histories
  const NumericVector R = dat["R"];            // received levels
  const NumericVector Y_rec = dat["Y_rec"];    // recorded bearings
  const NumericMatrix Y_grid = dat["Y_grid"];  // grid bearings
  const NumericMatrix X = dat["X"];            // distances
  const NumericVector A_x = dat["A_x"];        // areas of grid cells of mesh
  const int trunc_level = dat["trunc_level"];  // truncation level for received levels
  const NumericVector D = dat["D"];            // density for every grid point
  
  const int USE_BEARINGS = dat["USE_BEARINGS"];// 0: no bearings, 1: single
                                               // kappa; 2: mix kappa
  // Declare all to allow compiling
  double bessel, bessel_low, bessel_high;   
  
  if (USE_BEARINGS == 1) {
   bessel = dat["bessel"]; // the BesselI() values
  }
  if (USE_BEARINGS == 2) {
   bessel_low = dat["bessel_low"]; // the BesselI() values
   bessel_high = dat["bessel_high"];
   
  }
  // ===========================================================================
  
  // ===========================================================================
  // 2. EXTRACT INDIVIDUAL PARAMETERS
  // ---------------------------------------------------------------------------
  const double g0 = exp(double(par["logit_g0"])) / 
    (exp(double(par["logit_g0"])) + 1.0);
  const double beta_r = exp(double(par["log_beta_r"]));
  const double sd_r = exp(double(par["log_sd_r"]));
  const double mu_s = exp(double(par["log_mu_s"]));

  // Declare all at 0.0 to enable compiling
  double kappa, kappa_low, kappa_high, mix_bear;

  if (USE_BEARINGS == 1) {
    kappa = exp(double(par["log_kappa"]));
  }
  if (USE_BEARINGS == 2) {
    kappa_low = exp(double(par["log_kappa_low"]));
    kappa_high = exp(double(par["log_kappa_high"])) + kappa_low;
    mix_bear = exp(double(par["logit_mix_bear"])) / 
      (exp(double(par["logit_mix_bear"])) + 1.0);
  }
  
  // const NumericVector dens_pars = par["density_pars"];
  
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
  const int n_det =  dat["n_det"];
  const int n_grid =  dat["n_grid"];

  // std::cout << "Number of observations: " << n_call << std::endl;
  // std::cout << "Number of detectors: " << n_det << std::endl;
  // std::cout << "Grid size of the mesh: " << n_grid << std::endl;
  // std::cout << "Number of source level values tried: " << n_sl << std::endl;
  // std::cout << "kappa_high parameter: " << kappa_high << std::endl;
  
  // std::cout<< "1 + INFINITY: " << 1 + INFINITY << std::endl; // Is still inf!
  
  // std::cout << "test: " << - pow(10, 16) << std::endl;
  
  // std::cout << "X(0, 0): " << X(0, 0) << std::endl;
  
  // ===========================================================================
  // THE CONDITIONAL PART OF THE LIKELIHOOD
  // ===========================================================================
  
  // // Define the conditional negative log likelihood
  // double llk_cond = 0.0; // define vector to store likelihoods per individual call
  
  NumericVector log_probs_i(n_grid);
  for (int m = 0; m < n_grid; m++) { // loop through the mesh
    // =====================================================================
    // Likelihood of capture histories
    // ---------------------------------------------------------------------
    double log_prob_w_im = 0.0;
    
    for (int j = 0; j < n_det; j++) {
      double expected_rl = mu_s - beta_r * log10(X(m, j));
      double xx = (trunc_level - expected_rl) / sd_r;
      double prob = g0 * (1.0 - R::pnorm(xx, 0.0, 1.0, true, false));
      
      int test = 1 * W[j]; // Multiply W[j]  by one so to not use DATA directly
      
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
    
    if (USE_BEARINGS != 0) {
      for (int j = 0; j < n_det; j++) {
        if (W[j] == 1) {
          double cos_obs_minus_exp = cos(Y_rec[j] - Y_grid(m, j));
          if (USE_BEARINGS == 1) { // REPLACE LATER FOR MIXTURE CHECK
            log_prob_y_im += kappa * cos_obs_minus_exp -
              log(2 * M_PI * bessel);// log version
          } 
          if (USE_BEARINGS == 2) {
            log_prob_y_im += log(mix_bear * exp(kappa_low * cos_obs_minus_exp) / 
              (2.0 * M_PI * bessel_low) +
              (1.0 - mix_bear) * exp(kappa_high * cos_obs_minus_exp) / 
              (2.0 * M_PI * bessel_high));
          }
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
      if (W[j] == 1) {
        double expected_rl = mu_s - beta_r * log10(X(m, j));
        double recorded_rl = 1.0 * R[j];
        
        // Probability of exp_rl - rec_rl
        double log_prob = R::dnorm(recorded_rl, expected_rl, sd_r, true); // log probabilities
        
        log_prob -= log((1.0 - // -= log(sd_r * 1.0 -
          R::pnorm(trunc_level, expected_rl, sd_r, true, false))); // Different from R code,
        // but maths is the same <- NO, this is not true.
        // 09/09/2021: where is the detection probability in the denominator (see notebook)?
        
        log_prob_r_im += log_prob;
      }
    }
    // Set log_prob_r_im to -pow(10, 16) if it is -inf
    if (log_prob_r_im == -INFINITY) log_prob_r_im = -pow(10, 16);
    if (log_prob_r_im == INFINITY) log_prob_r_im = pow(10, 15); // 16 - 1 = 15, just to create a difference
    // =====================================================================
    
    // if (m == 0) {
    //   std::cout << "log(D): " << log(D[m]) << std::endl; // works
    //   std::cout << "log_prob_w_im: " << log_prob_w_im << std::endl; // works for now
    //   std::cout << "log_prob_y_im: " << log_prob_y_im << std::endl; // works
    //   std::cout << "log_prob_r_im: " << log_prob_r_im << std::endl; // works
    // }

    // log likelihood of a single call given a single location
    double cond_llk_i_m = 
      log(A_x[m]) + 
      log(D[m]) + 
      log_prob_w_im + 
      log_prob_y_im + 
      log_prob_r_im;
    
    bool test = std::isnan(cond_llk_i_m); // isnan from cmath requires a double
    if (test)  {
      std::cout << "log_prob_w_im: " << log_prob_w_im <<
        "\nlog_prob_y_im: " << log_prob_y_im <<
          "\nlog_prob_r_im: " << log_prob_r_im << std::endl;
    }
    
    log_probs_i[m] = cond_llk_i_m;
  }
  // Sum the values in log_probs_i_sl for every grid point, in a similar way to logsumexp()
  // Can we use the function logspace_add() ?
  // Sum the values in log_probs_i for every source level, in a similar way to logsumexp()
  double log_sum_i = log_probs_i[0];
  for (int ii = 1; ii < n_grid; ii++) {
    log_sum_i = R::logspace_add(log_sum_i, log_probs_i[ii]);
  }
  // Make it negative
  double nll_cond = -log_sum_i;
  
  // std::cout << "The conditional part of the log likelihood: " << llk_cond << std::endl; // nan
  // std::cout << "The negative log likelihood: " << nll << std::endl;
  
  // Return the nll
  return nll_cond;
}

// =============================================================================
