#ifndef IMPORTANCE_RESAMPLING_WITHOUT_REPLACEMENT_S2_HEADER_GUARD
#define IMPORTANCE_RESAMPLING_WITHOUT_REPLACEMENT_S2_HEADER_GUARD
#include "Rcpp.h"
namespace particleMethodsBernoulli
{
	SEXP importanceResamplingWithoutReplacementS2(SEXP lowerBound, SEXP trueProbabilities, SEXP n, SEXP seed, SEXP k);
}
#endif
