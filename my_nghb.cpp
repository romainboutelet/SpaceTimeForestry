#include <Rcpp.h>
#include <cmath>
#include <limits>

// [[Rcpp::export]]
Rcpp::IntegerVector my_nghb_cpp(Rcpp::NumericVector x, Rcpp::NumericMatrix X, Rcpp::NumericVector D) {
  int num_points = X.nrow();
  int dim_x = X.ncol();
  Rcpp::NumericVector D_tmp(Rcpp::clone(D));
  
  // Replace D == 0 with Inf
  for (int i = 0; i < num_points; ++i) {
    if (D_tmp(i) == 0) {
      D_tmp(i) = std::numeric_limits<double>::infinity();
    }
  }

  
  // Find the index of the minimum distance
  int nghb = 0;
  double min_dist = D_tmp(0);
  for (int i = 1; i < num_points; ++i) {
    if (D_tmp(i) < min_dist) {
      min_dist = D_tmp(i);
      nghb = i;
    }
  }
  int nghb_r = nghb + 1;
  
  // Initialize the neighbors list
  Rcpp::IntegerVector nghb_indices = Rcpp::IntegerVector::create(nghb_r);
  Rcpp::NumericVector x0 = X.row(nghb);
  // Mask initialization
  Rcpp::LogicalVector mask(num_points, true);
  for (int i = 0; i < num_points; ++i) {
    mask[i] = (X(i, 0) * (x[0] - x0[0]) > (x[0] * x[0] + x[1] * x[1] - x0[0] * x0[0] - x0[1] * x0[1]) / 2 - X(i, 1) * (x[1] - x0[1]));
  }
  D_tmp(nghb) = std::numeric_limits<double>::infinity();  // set the distance of the chosen nghb to Inf
  
  int count = 1;
  
  // Main loop
  while (Rcpp::any(mask).is_true() && count < 5) {
    // Update D based on mask
    for (int i = 0; i < num_points; ++i) {
      if (!mask[i]) {
        D_tmp(i) = std::numeric_limits<double>::infinity();
      }
    }
    
    // Find the next minimum distance
    nghb = 0;
    min_dist = D_tmp(0);
    for (int i = 1; i < num_points; ++i) {
      if (D_tmp(i) < min_dist) {
        min_dist = D_tmp(i);
        nghb = i;
      }
    }
    nghb_r = nghb + 1;
    
    // Add the new index to nghb_indices
    nghb_indices.push_back(nghb_r);
    x0 = X.row(nghb);
    
    // Update the mask based on the new chosen nghb
    for (int i = 0; i < num_points; ++i) {
      mask[i] = (mask[i] && (X(i, 0) * (x[0] - x0[0]) > (x[0] * x[0] + x[1] * x[1] - x0[0] * x0[0] - x0[1] * x0[1]) / 2 - X(i, 1) * (x[1] - x0[1])));
    }
    
    count++;
  }
  
  return nghb_indices;
}
