#ifndef IMPORTANCE_RESAMPLING_HEADER_GUARD
#define IMPORTANCE_RESAMPLING_HEADER_GUARD
#include "Rcpp.h"
namespace particleMethodsBernoulli
{
	SEXP importanceResampling(SEXP nBernoullis, SEXP lowerBound, SEXP trueProbability, SEXP n, SEXP seed);
}
#endif
