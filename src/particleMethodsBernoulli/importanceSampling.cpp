#include "importanceSampling.h"
#include "includeMPFR.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/bernoulli_distribution.hpp>
namespace particleMethodsBernoulli
{
	SEXP importanceSampling(SEXP lowerBound_sexp, SEXP trueProbabilities_sexp, SEXP n_sexp, SEXP seed_sexp)
	{
	BEGIN_RCPP
		std::vector<double> trueProbabilities;
		try
		{
			trueProbabilities = Rcpp::as<std::vector<double> >(trueProbabilities_sexp);
		}
		catch(...)
		{
			throw std::runtime_error("Input trueProbabilities must be a numeric vector");
		}
		for(std::vector<double>::iterator trueProbability = trueProbabilities.begin(); trueProbability != trueProbabilities.end(); trueProbability++)
		{
			if(*trueProbability <= 0 || *trueProbability >= 1)
			{
				throw std::runtime_error("Input trueProbability must be in (0, 1)");
			}
		}
		int nBernoullis = trueProbabilities.size();

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
		std::vector<double> ratio1, ratio2;
		for(std::vector<double>::iterator trueProbability = trueProbabilities.begin(); trueProbability != trueProbabilities.end(); trueProbability++)
		{
			ratio1.push_back(*trueProbability / newProbability);
			ratio2.push_back((1 - *trueProbability) / (1 - newProbability));
		}
		boost::random::bernoulli_distribution<> bernoulli(newProbability);
		mpfr_class sum = 0, sumSquared = 0;
		for(int i = 0; i < n; i++)
		{
			mpfr_class likelihoodRatio = 1;
			int valuesOfOne = 0;
			for(int bernoulliCounter = 0; bernoulliCounter < nBernoullis; bernoulliCounter++)
			{
				int value = bernoulli(randomSource);
				if(value)
				{
					valuesOfOne++;
					likelihoodRatio *= ratio1[bernoulliCounter];
				}
				else
				{
					likelihoodRatio *= ratio2[bernoulliCounter];
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
