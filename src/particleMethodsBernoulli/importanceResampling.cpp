#include "importanceResampling.h"
#include "includeMPFR.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
namespace particleMethodsBernoulli
{
	SEXP importanceResampling(SEXP lowerBound_sexp, SEXP trueProbabilities_sexp, SEXP n_sexp, SEXP seed_sexp)
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
		std::vector<int> weightRatio1, weightRatio2;
		mpfr_class averageWeight = 1;
		for(int bernoulliCounter = 0; bernoulliCounter < nBernoullis; bernoulliCounter++)
		{
			weightRatio1.clear();
			weightRatio2.clear();
			//Sample and update the weights. Everything in weightRatio1 has a weight of averageWeight * ratio1. Everything in weightRatio2 has a weight of averageWeight * ratio2. Anything in neither has a weight of 0. 
			for(int i = 0; i < n; i++)
			{
				int value = bernoulli(randomSource);
				if(value)
				{
					valuesOfOne[i]++;
					weightRatio1.push_back(i);
				}
				else
				{
					if(valuesOfOne[i] + nBernoullis - 1 - bernoulliCounter > lowerBound)
					{
						weightRatio2.push_back(i);
					}
				}
			}
			//Resampling step
			if(weightRatio1.size() == 0 && weightRatio2.size() == 0)
			{
				averageWeight = 0;
				goto returnAnswer;
			}
			else if(weightRatio1.size() == 0)
			{
				boost::random::uniform_int_distribution<> ratio2Sampler(0, weightRatio2.size()-1);
				for(int i = 0; i < n; i++)
				{
					newValuesOfOne[i] = valuesOfOne[weightRatio2[ratio2Sampler(randomSource)]];
				}
			}
			else if(weightRatio2.size() == 0)
			{
				boost::random::uniform_int_distribution<> ratio1Sampler(0, weightRatio1.size()-1);
				for(int i = 0; i < n; i++)
				{
					newValuesOfOne[i] = valuesOfOne[weightRatio1[ratio1Sampler(randomSource)]];
				}
			}
			else
			{
				//With this probability we select something with weight ratio1
				boost::random::bernoulli_distribution<> bernoulliForResampling(weightRatio1.size() * ratio1[bernoulliCounter] / (weightRatio1.size() * ratio1[bernoulliCounter] + weightRatio2.size() * ratio2[bernoulliCounter]));
				boost::random::uniform_int_distribution<> ratio2Sampler(0, weightRatio2.size()-1);
				boost::random::uniform_int_distribution<> ratio1Sampler(0, weightRatio1.size()-1);
				for(int i = 0; i < n; i++)
				{
					if(bernoulliForResampling(randomSource))
					{
						newValuesOfOne[i] = valuesOfOne[weightRatio1[ratio1Sampler(randomSource)]];
					}
					else
					{
						newValuesOfOne[i] = valuesOfOne[weightRatio2[ratio2Sampler(randomSource)]];
					}
				}
			}
			valuesOfOne.swap(newValuesOfOne);
			//Update average weight
			averageWeight *= (weightRatio1.size() * ratio1[bernoulliCounter] + weightRatio2.size() * ratio2[bernoulliCounter]) / n;
		}
returnAnswer:
		return Rcpp::List::create(Rcpp::Named("estimate") = averageWeight.convert_to<double>());
	END_RCPP
	}
}
