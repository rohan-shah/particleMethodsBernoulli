#include "importanceSampling.h"
#include "includeMPFR.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/bernoulli_distribution.hpp>
namespace particleMethodsBernoulli
{
	SEXP importanceSampling(SEXP nBernoullis_sexp, SEXP lowerBound_sexp, SEXP trueProbability_sexp, SEXP n_sexp, SEXP seed_sexp)
	{
	BEGIN_RCPP
		int nBernoullis;
		try
		{
			nBernoullis = Rcpp::as<int>(nBernoullis_sexp);
		}
		catch(...)
		{
			throw std::runtime_error("Input nBernoullis must be an integer");
		}
		if(nBernoullis < 0)
		{
			throw std::runtime_error("Input nBernoullis must be positive");
		}

		int lowerBound;
		try
		{
			lowerBound = Rcpp::as<int>(lowerBound_sexp);
		}
		catch(...)
		{
			throw std::runtime_error("Input lowerBound must be an integer");
		}
		if(lowerBound <= nBernoullis/2)
		{
			throw std::runtime_error("Input lowerBound must be bigger than nBernoullis/2");
		}
		if(lowerBound >= nBernoullis)
		{
			throw std::runtime_error("Input lowerBound must be smaller than nBernoullis");
		}

		double trueProbability;
		try
		{
			trueProbability = Rcpp::as<double>(trueProbability_sexp);
		}
		catch(...)
		{
			throw std::runtime_error("Input trueProbability must be a number");
		}
		if(trueProbability <= 0 || trueProbability >= 1)
		{
			throw std::runtime_error("Input trueProbability must be in (0, 1)");
		}

		int n;
		try
		{
			n = Rcpp::as<int>(n_sexp);
		}
		catch(...)
		{
			throw std::runtime_error("Input n must be an integer");
		}
		if(n < 1)
		{
			throw std::runtime_error("Input n must be positive");
		}

		int seed;
		try
		{
			seed = Rcpp::as<int>(seed_sexp);
		}
		catch(...)
		{
			throw std::runtime_error("Input seed must be an integer");
		}

		boost::mt19937 randomSource;
		randomSource.seed(seed);
		double newProbability = (double)lowerBound / (double)nBernoullis;
		double ratio1 = trueProbability / newProbability, ratio2 = (1 - trueProbability) / (1 - newProbability);
		boost::random::bernoulli_distribution<> bernoulli(newProbability);
		mpfr_class sum = 0, sumSquared = 0;
		for(int i = 0; i < n; i++)
		{
			double likelihoodRatio = 1;
			int valuesOfOne = 0;
			for(int bernoulliCounter = 0; bernoulliCounter < nBernoullis; bernoulliCounter++)
			{
				int value = bernoulli(randomSource);
				if(value)
				{
					valuesOfOne++;
					likelihoodRatio *= ratio1;
				}
				else
				{
					likelihoodRatio *= ratio2;
				}
				//Break if we can never be above the lower bound
				if(valuesOfOne + nBernoullis - 1 - bernoulliCounter <= lowerBound) break;
			}
			if(valuesOfOne > lowerBound)
			{
				sum += likelihoodRatio;
				sumSquared += likelihoodRatio*likelihoodRatio;
			}
		}
		mpfr_class estimate = sum/n;
		mpfr_class estimateOfVariance = (sumSquared / n - estimate * estimate)/n;
		return Rcpp::List::create(Rcpp::Named("estimate") = estimate.convert_to<double>(), Rcpp::Named("estimateOfVariance") = estimateOfVariance.convert_to<double>());
	END_RCPP
	}
}
