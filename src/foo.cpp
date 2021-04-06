#include <Rcpp.h>
using namespace Rcpp;

#define IS_GFF  (1<<0)

// [[Rcpp::export]]
double one(double x) {
  printf("%f", NA_REAL);
  if (NumericVector::is_na(x))
    return 12.28;
  else
    return x + 2.2;
}



