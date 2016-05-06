#include "importanceResamplingWithoutReplacement.h"
#include "includeMPFR.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
namespace particleMethodsBernoulli
{
	SEXP importanceResamplingWithoutReplacement(SEXP nBernoullis_sexp, SEXP lowerBound_sexp, SEXP trueProbability_sexp, SEXP n_sexp, SEXP seed_sexp)
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
		double ratio1 = newProbability / trueProbability, ratio2 = (1 - newProbability) / (1 - trueProbability);
		boost::random::bernoulli_distribution<> bernoulli(newProbability);
		//Initially we have two samples, corresponding to the first bernoulli being 0 or 1. Note that nBernoullis == 1 gives an error above, so we can assume that there are at least 2 bernoullis
		std::vector<int> samples, newSamples;
		std::vector<std::pair<int, int> > sampleWeights, newSampleWeights;
		samples.push_back(0);
		samples.push_back(1);
		sampleWeights.push_back(std::make_pair(0, 0));
		sampleWeights.push_back(std::make_pair(0, 0));
		std::vector<int> choicesUp, choicesDown;

		mpfr_class product = 1;
		for(int bernoulliCounter = 1; bernoulliCounter < nBernoullis; bernoulliCounter++)
		{
			choicesUp.clear();
			choicesDown.clear();
			//Sample and update the weights. Everything in weightRatio1 has a weight of averageWeight * ratio1. Everything in weightRatio2 has a weight of averageWeight * ratio2. Anything in neither has a weight of 0. 
			for(std::size_t i = 0; i < samples.size(); i++)
			{
				int maxPossible = samples[i] + nBernoullis - bernoulliCounter;
				if(maxPossible == lowerBound + 1)
				{
					choicesUp.push_back((int)i);
				}
				else if(maxPossible > lowerBound + 1)
				{
					choicesUp.push_back((int)i);
					choicesDown.push_back((int)i);
				}
			}
			newSamples.clear();
			newSampleWeights.clear();
			//We take every successor
			if(choicesUp.size() + choicesDown.size() <= (std::size_t)n)
			{
				for(std::size_t i = 0; i < choicesDown.size(); i++)
				{
					newSamples.push_back(samples[choicesDown[i]]);
					newSampleWeights.push_back(sampleWeights[choicesDown[i]]);
				}
				for(std::size_t i = 0; i < choicesUp.size(); i++)
				{
					newSamples.push_back(samples[choicesUp[i]]+1);
					newSampleWeights.push_back(sampleWeights[choicesUp[i]]);
				}
			}
			else
			{
				throw std::runtime_error("Not Implemented");
			}
			samples.swap(newSamples);
			sampleWeights.swap(newSampleWeights);
		}
		std::vector<int> table((nBernoullis+1)*(nBernoullis+1), -1);
		std::vector<mpfr_class> powerPairs, densityValues(nBernoullis+1);
		mpfr_class power1 = trueProbability, power2 = 1 - trueProbability;
		for(int i = 0; i <= nBernoullis; i++)
		{
			densityValues[i] = boost::multiprecision::pow(power1, i) * boost::multiprecision::pow(power2, nBernoullis-i);
		}

		mpfr_class estimate = 0;
		mpfr_class ratio1_mpfr = ratio1;
		mpfr_class ratio2_mpfr = ratio2;
		for(std::size_t i = 0; i < samples.size(); i++)
		{
			std::pair<int, int> weight = sampleWeights[i];
			if(table[weight.first + weight.second * (nBernoullis+1)] == -1)
			{
				mpfr_class power = boost::multiprecision::pow(ratio1_mpfr, weight.first) * boost::multiprecision::pow(ratio2_mpfr, weight.second);
				table[weight.first + weight.second * (nBernoullis+1)] = powerPairs.size();
				powerPairs.push_back(power);
			}
			estimate += powerPairs[table[weight.first + weight.second*(nBernoullis+1)]] * densityValues[samples[i]];
		}
		return Rcpp::List::create(Rcpp::Named("estimate") = estimate.convert_to<double>());
	END_RCPP
	}
}
