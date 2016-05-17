#include "importanceResamplingAuxiliary.h"
#include "includeMPFR.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
namespace particleMethodsBernoulli
{
	SEXP importanceResamplingAuxiliary(SEXP lowerBound_sexp, SEXP trueProbabilities_sexp, SEXP n_sexp, SEXP seed_sexp)
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
		std::vector<int> valuesOfOne(n, 0);
		std::vector<int> newValuesOfOne(n, 0);
		std::vector<int> toSample;
		std::vector<mpfr_class> weights(n, 1), newWeights(n, 1);

		mpfr_class product = 1;
		for(int bernoulliCounter = 0; bernoulliCounter < nBernoullis; bernoulliCounter++)
		{
			toSample.clear();
			//Sample and update the weights. Everything in weightRatio1 has a weight of averageWeight * ratio1. Everything in weightRatio2 has a weight of averageWeight * ratio2. Anything in neither has a weight of 0. 
			for(int i = 0; i < n; i++)
			{
				int value = bernoulli(randomSource);
				if(value)
				{
					valuesOfOne[i]++;
					weights[i] *= ratio1[bernoulliCounter];
					toSample.push_back(i);
				}
				else
				{
					if(valuesOfOne[i] + nBernoullis - 1 - bernoulliCounter > lowerBound)
					{
						toSample.push_back(i);
						weights[i] *= ratio2[bernoulliCounter];
					}
				}
			}
			product *= (double) toSample.size() / (double)n;
			//Resampling step
			if(toSample.size() == 0)
			{
				return Rcpp::List::create(Rcpp::Named("estimate") = 0.0);
			}
			else
			{
				//With this probability we select something with weight ratio1
				boost::random::uniform_int_distribution<> simpleRandom(0, toSample.size()-1);
				for(int i = 0; i < n; i++)
				{
					int sampled = toSample[simpleRandom(randomSource)];
					newValuesOfOne[i] = valuesOfOne[sampled];
					newWeights[i] = weights[sampled];
				}
			}
			valuesOfOne.swap(newValuesOfOne);
			weights.swap(newWeights);
		}
		mpfr_class sum = 0;
		for(int i = 0; i < n; i++)
		{
			sum += weights[i];
		}
		sum /= n;
		sum *= product;
		return Rcpp::List::create(Rcpp::Named("estimate") = sum.convert_to<double>());
	END_RCPP
	}
}
