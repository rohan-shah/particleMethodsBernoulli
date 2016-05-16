#ifndef IMPORTANCE_RESAMPLING_AUXILIARY_HEADER_GUARD
#define IMPORTANCE_RESAMPLING_AUXILIARY_HEADER_GUARD
#include "Rcpp.h"
namespace particleMethodsBernoulli
{
	SEXP importanceResamplingAuxiliary(SEXP lowerBound, SEXP trueProbabilities, SEXP n, SEXP seed);
}
#endif
