#ifndef IMPORTANCE_SAMPLING_HEADER_GUARD
#define IMPORTANCE_SAMPLING_HEADER_GUARD
#include "Rcpp.h"
namespace particleMethodsBernoulli
{
	SEXP importanceSampling(SEXP nBernoullis, SEXP lowerBound, SEXP trueProbability, SEXP n, SEXP seed);
}
#endif
