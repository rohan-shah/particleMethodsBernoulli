#ifndef IMPORTANCE_RESAMPLING_WITHOUT_REPLACEMENT_HEADER_GUARD
#define IMPORTANCE_RESAMPLING_WITHOUT_REPLACEMENT_HEADER_GUARD
#include "Rcpp.h"
namespace particleMethodsBernoulli
{
	SEXP importanceResamplingWithoutReplacement(SEXP nBernoullis, SEXP lowerBound, SEXP trueProbability, SEXP n, SEXP seed);
}
#endif
