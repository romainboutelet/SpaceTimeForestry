#include <Rcpp.h>
#include <stdio.h>

// [[Rcpp::export]]
Rcpp::IntegerVector my_nghb_cpp(Rcpp::NumericVector x, Rcpp::NumericMatrix X, Rcpp::NumericVector D) {
  int num_points = X.nrow();
  Rcpp::LogicalVector mask(num_points, true);
  
  // Replace D == 0 with Inf
  for (int i = 0; i < num_points; ++i) {
    if (D(i) == 0) {
      mask(i) = 0;
    }
  }
  // Find the index of the minimum distance
  int nghb = 0;
  while (!mask(nghb)) {
    nghb++;
  }
  double min_dist = D(nghb);
  for (int i = 1; i < num_points; ++i) {
    if ((D(i) < min_dist)&& mask[i]) {
      min_dist = D(i);
      nghb = i;
    }
  }
  int nghb_r = nghb + 1;
  
  // Initialize the neighbors list
  Rcpp::IntegerVector nghb_indices (5);
  nghb_indices(0) = nghb_r;
  // Mask initialization
  bool cond = 0;
  for (int i = 0; i < num_points; ++i) {
    mask[i] = mask(i) && (X(i, 0) * (x[0] - X.row(nghb)[0]) > (x[0] * x[0] + x[1] *
      x[1] - X.row(nghb)[0] * X.row(nghb)[0] - X.row(nghb)[1] *
      X.row(nghb)[1]) / 2 - X(i, 1) * (x[1] - X.row(nghb)[1]));
    cond = mask[i] || cond;
  }
  int count = 1;
  // Main loop
  while (cond == 1) {
    // Find the next minimum distance
    nghb = 0;
    while (!mask(nghb)) {
      nghb++;
    }
    min_dist = D(nghb);
    for (int i = 1; i < num_points; ++i) {
      if ((D(i) < min_dist) && mask[i]) {
        min_dist = D(i);
        nghb = i;
      }
    }
    nghb_r = nghb + 1;
    // Add the new index to nghb_indices
    nghb_indices(count) = nghb_r;

    // Update the mask based on the new chosen nghb
    cond = 0;
    for (int i = 0; i < num_points; ++i) {
      mask[i] = mask[i] && ((X(i, 0) * (x[0] - X.row(nghb)[0]) >
                               (x[0] * x[0] + x[1] * x[1] -
                               X.row(nghb)[0] * X.row(nghb)[0] - 
                               X.row(nghb)[1] * X.row(nghb)[1]) / 2 -
                               X(i, 1) * (x[1] - X.row(nghb)[1])));
      cond = mask[i] || cond;
    }
    count++;
    cond = cond && (count < 5);
    if (cond == 0){
      break;
    }
  }
  int n_nghb = 0;
  for (int i = 0; i < 5; ++i) {
    if (nghb_indices(i) != 0) {
      n_nghb++;
    }
  }
  Rcpp::IntegerVector nghb_idx(n_nghb);
  for (int i = 0; i < n_nghb; ++i) {
    nghb_idx(i) = nghb_indices(i);
  }
  return nghb_idx;
}

