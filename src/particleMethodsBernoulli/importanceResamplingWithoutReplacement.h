#ifndef IMPORTANCE_RESAMPLING_WITHOUT_REPLACEMENT_HEADER_GUARD
#define IMPORTANCE_RESAMPLING_WITHOUT_REPLACEMENT_HEADER_GUARD
#include "Rcpp.h"
namespace particleMethodsBernoulli
{
	SEXP importanceResamplingWithoutReplacement(SEXP lowerBound, SEXP trueProbabilities, SEXP n, SEXP seed);
}
#endif
