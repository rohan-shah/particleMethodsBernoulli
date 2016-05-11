#ifndef IMPORTANCE_RESAMPLING_HEADER_GUARD
#define IMPORTANCE_RESAMPLING_HEADER_GUARD
#include "Rcpp.h"
namespace particleMethodsBernoulli
{
	SEXP importanceResampling(SEXP lowerBound, SEXP trueProbabilities, SEXP n, SEXP seed);
}
#endif
