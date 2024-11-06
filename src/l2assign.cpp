#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
IntegerVector l2assign(NumericMatrix spec, NumericMatrix p) {
  int nsubpop = p.nrow();
  int nbands = (p.ncol() + 1) / 2;
  int nfreq = spec.nrow();
  int nrep = spec.ncol();

  // assignment vector
  IntegerVector result(nrep);

  // distance matrix
  NumericMatrix dist(nrep, nsubpop);

  // loop over subpopulations
  for (int j = 0; j < nsubpop; j++) {
    // boundaries and collapsed measures for the jth subpopulation encoded
    NumericVector endpoints = p(j, _);
    NumericVector collapsed = p(j, _);
    endpoints = endpoints[Range(0, nbands - 2)];
    collapsed = collapsed[Range(nbands - 1, p.ncol() - 1)];

    // Compute widths of each band and expand collapsed measures
    std::vector<int> widths(nbands);
    for (int i = 0; i < nbands; i++) {
      if (i == 0) {
        widths[i] = endpoints[i];
      } else if (i == nbands - 1) {
        widths[i] = (nfreq + 1) - endpoints[i - 1];
      } else {
        widths[i] = endpoints[i] - endpoints[i - 1];
      }
    }

    // Create center vector by repeating collapsed values according to band widths
    std::vector<double> center_vec(nfreq);
    for (int i = 0; i < nbands; i++) {
      for (int w = 0; w < widths[i]; w++) {
        center_vec[w] = collapsed[i];
      }
    }

    // Convert center vector back to a NumericVector
    // NumericVector center(center_vec.begin(), center_vec.end());

    // Compute L2 distance between each spectrum and the j-th cluster
    for (int i = 0; i < nrep; i++) {
      double dist_val = 0;
      for (int k = 0; k < nfreq; k++) {
        dist_val += pow(spec(k, i) - center_vec[k], 2);
      }
      dist(i, j) = dist_val / (2 * (nfreq + 1));
    }
  }

  // Assign the minimum distance index for each spectrum
  for (int i = 0; i < nrep; i++) {
    double min_dist = dist(i, 0);
    int min_idx = 0;
    for (int j = 1; j < nsubpop; j++) {
      if (dist(i, j) < min_dist) {
        min_dist = dist(i, j);
        min_idx = j;
      }
    }
    result[i] = min_idx + 1; // R uses 1-based indexing
  }

  return result;
}
