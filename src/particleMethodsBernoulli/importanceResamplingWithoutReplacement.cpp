#include "importanceResamplingWithoutReplacement.h"
#include "includeMPFR.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/random_number_generator.hpp>
#include <boost/range/algorithm/random_shuffle.hpp>
namespace particleMethodsBernoulli
{
	SEXP importanceResamplingWithoutReplacement(SEXP lowerBound_sexp, SEXP trueProbabilities_sexp, SEXP n_sexp, SEXP seed_sexp)
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
		boost::random_number_generator<boost::mt19937> generator(randomSource);

		double newProbability = (double)lowerBound / (double)nBernoullis;
		std::vector<double> ratio1, ratio2, totals;
		for(std::vector<double>::iterator trueProbability = trueProbabilities.begin(); trueProbability != trueProbabilities.end(); trueProbability++)
		{
			ratio1.push_back(*trueProbability / newProbability);
			ratio2.push_back((1 - *trueProbability) / (1 - newProbability));
			totals.push_back(ratio1.back() + ratio2.back());
		}
		//Initially we have two samples, corresponding to the first bernoulli being 0 or 1. Note that nBernoullis == 1 gives an error above, so we can assume that there are at least 2 bernoullis
		std::vector<int> samples, newSamples;
		std::vector<mpfr_class> sampleDensityOnWeights, newSampleDensityOnWeights;
		std::vector<mpfr_class> remainingDensityValues, newRemainingDensityValues;

		samples.push_back(0);
		samples.push_back(1);
		sampleDensityOnWeights.push_back(1);
		sampleDensityOnWeights.push_back(1);
		remainingDensityValues.push_back(1-trueProbabilities[0]);
		remainingDensityValues.push_back(trueProbabilities[0]);
		std::vector<int> choicesUp, choicesDown;

		mpfr_class product = 1;
		for(int bernoulliCounter = 1; bernoulliCounter < nBernoullis; bernoulliCounter++)
		{
			choicesUp.clear();
			choicesDown.clear();
			//Sample and update the weights. Everything in weightRatio1 has a weight of averageWeight * ratio1[bernoulliCounter]. Everything in weightRatio2 has a weight of averageWeight * ratio2[bernouliCounter]. Anything in neither has a weight of 0. 
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
			newSampleDensityOnWeights.clear();
			newRemainingDensityValues.clear();
			double complementaryTrueProb = 1 - trueProbabilities[bernoulliCounter];
			//We take every successor
			if(choicesUp.size() + choicesDown.size() <= (std::size_t)n)
			{
				for(std::size_t i = 0; i < choicesDown.size(); i++)
				{
					newSamples.push_back(samples[choicesDown[i]]);
					newSampleDensityOnWeights.push_back(sampleDensityOnWeights[choicesDown[i]]);
					newRemainingDensityValues.push_back(remainingDensityValues[choicesDown[i]] * complementaryTrueProb);
				}
				for(std::size_t i = 0; i < choicesUp.size(); i++)
				{
					newSamples.push_back(samples[choicesUp[i]]+1);
					newSampleDensityOnWeights.push_back(sampleDensityOnWeights[choicesUp[i]]);
					newRemainingDensityValues.push_back(remainingDensityValues[choicesUp[i]] * trueProbabilities[bernoulliCounter]);
				}
			}
			else
			{
				boost::random_shuffle(choicesUp, generator);
				boost::random_shuffle(choicesDown, generator);
				int sampled = 0;
				int nPairs = choicesUp.size() + choicesDown.size() - n;
				for(; sampled < nPairs; sampled++)
				{
					mpfr_class upSize_mpfr = remainingDensityValues[choicesUp[sampled]] * trueProbabilities[bernoulliCounter];
					mpfr_class downSize_mpfr = remainingDensityValues[choicesDown[sampled]] * complementaryTrueProb;
					double upSize = upSize_mpfr.convert_to<double>(), downSize = downSize_mpfr.convert_to<double>();
					double total = upSize + downSize;
					double upProbability = upSize / total;
					boost::random::bernoulli_distribution<> raoHartleyCochranBernoulli(upProbability);
					if(raoHartleyCochranBernoulli(randomSource))
					{
						newSamples.push_back(samples[choicesUp[sampled]]+1);
						newSampleDensityOnWeights.push_back(sampleDensityOnWeights[choicesUp[sampled]] * upSize_mpfr / upProbability);
						newRemainingDensityValues.push_back(1);
					}
					else
					{
						newSamples.push_back(samples[choicesDown[sampled]]);
						newSampleDensityOnWeights.push_back(sampleDensityOnWeights[choicesDown[sampled]] * downSize_mpfr / (1 - upProbability));
						newRemainingDensityValues.push_back(1);
					}
				}
				for(; sampled < (int)choicesUp.size(); sampled++)
				{
					newSamples.push_back(samples[choicesUp[sampled]]+1);
					newSampleDensityOnWeights.push_back(sampleDensityOnWeights[choicesUp[sampled]]);
					newRemainingDensityValues.push_back(remainingDensityValues[choicesUp[sampled]] * trueProbabilities[bernoulliCounter]);
					if(sampled < (int)choicesDown.size()) 
					{
						newSamples.push_back(samples[choicesDown[sampled]]);
						newSampleDensityOnWeights.push_back(sampleDensityOnWeights[choicesDown[sampled]]);
						newRemainingDensityValues.push_back(remainingDensityValues[choicesDown[sampled]] * complementaryTrueProb);
					}
				}
				if((int)newSamples.size() != n) throw std::runtime_error("Internal error");
			}
			samples.swap(newSamples);
			sampleDensityOnWeights.swap(newSampleDensityOnWeights);
			remainingDensityValues.swap(newRemainingDensityValues);
		}
		mpfr_class estimate = 0;
		for(std::size_t i = 0; i < samples.size(); i++)
		{
			estimate += sampleDensityOnWeights[i] * remainingDensityValues[i];
		}
		return Rcpp::List::create(Rcpp::Named("estimate") = estimate.convert_to<double>());
	END_RCPP
	}
}
