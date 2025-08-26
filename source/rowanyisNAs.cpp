#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
LogicalVector rowanyisNAs(NumericMatrix mat) {
  int nrow = mat.nrow();
  int ncol = mat.ncol();
  LogicalVector out(nrow);
  
  for (int i = 0; i < nrow; ++i) {
    bool has_na = false;
    for (int j = 0; j < ncol; ++j) {
      if (NumericMatrix::is_na(mat(i, j))) {
        has_na = true;
        break;
      }
    }
    out[i] = has_na;
  }
  
  return out;
}