#ifndef IMPORTANCE_RESAMPLING_WITHOUT_REPLACEMENT_AUXILIARY_HEADER_GUARD
#define IMPORTANCE_RESAMPLING_WITHOUT_REPLACEMENT_AUXILIARY_HEADER_GUARD
#include "Rcpp.h"
namespace particleMethodsBernoulli
{
	SEXP importanceResamplingWithoutReplacementAuxiliary(SEXP lowerBound, SEXP trueProbabilities, SEXP n, SEXP seed);
}
#endif
