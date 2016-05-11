#ifndef IMPORTANCE_SAMPLING_HEADER_GUARD
#define IMPORTANCE_SAMPLING_HEADER_GUARD
#include "Rcpp.h"
namespace particleMethodsBernoulli
{
	SEXP importanceSampling(SEXP lowerBound, SEXP trueProbabilities, SEXP n, SEXP seed);
}
#endif
