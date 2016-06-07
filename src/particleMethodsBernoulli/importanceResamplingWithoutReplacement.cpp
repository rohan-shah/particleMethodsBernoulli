#include "importanceResamplingWithoutReplacement.h"
#include "includeMPFR.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/random_number_generator.hpp>
#include <boost/range/algorithm/random_shuffle.hpp>
#include "sampford.h"
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
		mpfr_class newProbability_mpfr = newProbability, compNewProbability_mpfr = 1 - newProbability;
		//Initially we have two samples, corresponding to the first bernoulli being 0 or 1. Note that nBernoullis == 1 gives an error above, so we can assume that there are at least 2 bernoullis
		std::vector<int> samples, newSamples;
		std::vector<mpfr_class> sampleDensityOnWeight, newSampleDensityOnWeight, sampfordWeights, newSampfordWeights, rescaledWeights;

		samples.push_back(0);
		samples.push_back(1);
		sampleDensityOnWeight.push_back(1-trueProbabilities[0]);
		sampleDensityOnWeight.push_back(trueProbabilities[0]);
		sampfordWeights.push_back(compNewProbability_mpfr);
		sampfordWeights.push_back(newProbability_mpfr);
		std::vector<int> choicesUp, choicesDown;

		sampling::sampfordFromParetoNaiveArgs sampfordArgs;
		sampfordArgs.n = n;
		std::vector<int> sampfordSampleIndices;
		std::vector<mpfr_class> sampfordSampleInclusionProbabilities;
		sampfordArgs.indices = &sampfordSampleIndices;
		sampfordArgs.inclusionProbabilities = &sampfordSampleInclusionProbabilities;
		sampfordArgs.weights = &newSampfordWeights;
		sampfordArgs.rescaledWeights = &rescaledWeights;

		for(int bernoulliCounter = 1; bernoulliCounter < nBernoullis; bernoulliCounter++)
		{
			choicesUp.clear();
			choicesDown.clear();
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
			newSampleDensityOnWeight.clear();
			newSampfordWeights.clear();
			double complementaryTrueProb = 1 - trueProbabilities[bernoulliCounter];
			//We take every successor
			if(choicesUp.size() + choicesDown.size() <= (std::size_t)n)
			{
				for(std::size_t i = 0; i < choicesDown.size(); i++)
				{
					newSamples.push_back(samples[choicesDown[i]]);
					newSampleDensityOnWeight.push_back(sampleDensityOnWeight[choicesDown[i]] * complementaryTrueProb);
					newSampfordWeights.push_back(sampfordWeights[choicesDown[i]]*compNewProbability_mpfr);
				}
				for(std::size_t i = 0; i < choicesUp.size(); i++)
				{
					newSamples.push_back(samples[choicesUp[i]]+1);
					newSampleDensityOnWeight.push_back(sampleDensityOnWeight[choicesUp[i]] * trueProbabilities[bernoulliCounter]);
					newSampfordWeights.push_back(sampfordWeights[choicesUp[i]]*newProbability_mpfr);
				}
				newSampfordWeights.swap(sampfordWeights);
			}
			else
			{
				for(std::size_t i = 0; i < choicesDown.size(); i++)
				{
					newSampfordWeights.push_back(sampfordWeights[choicesDown[i]] * compNewProbability_mpfr);
				}
				for(std::size_t i = 0; i < choicesUp.size(); i++)
				{
					newSampfordWeights.push_back(sampfordWeights[choicesUp[i]] * newProbability_mpfr);
				}
				sampling::sampfordFromParetoNaive(sampfordArgs, randomSource);
				sampfordWeights.clear();
				for(std::vector<int>::iterator j = sampfordSampleIndices.begin(); j != sampfordSampleIndices.end(); j++)
				{
					if(*j < (int)choicesDown.size())
					{
						newSamples.push_back(samples[choicesDown[*j]]);
						newSampleDensityOnWeight.push_back(sampleDensityOnWeight[choicesDown[*j]] * complementaryTrueProb / sampfordSampleInclusionProbabilities[*j]);
					}
					else
					{
						newSamples.push_back(samples[choicesUp[*j - choicesDown.size()]]+1);
						newSampleDensityOnWeight.push_back(sampleDensityOnWeight[choicesUp[*j - choicesDown.size()]] * trueProbabilities[bernoulliCounter]/ sampfordSampleInclusionProbabilities[*j]);
					}
					sampfordWeights.push_back(newSampfordWeights[*j] / sampfordSampleInclusionProbabilities[*j]);
				}
			}
			samples.swap(newSamples);
			sampleDensityOnWeight.swap(newSampleDensityOnWeight);
		}
		mpfr_class estimate = 0;
		for(std::size_t i = 0; i < samples.size(); i++)
		{
			estimate += sampleDensityOnWeight[i];
		}
		return Rcpp::List::create(Rcpp::Named("estimate") = estimate.convert_to<double>());
	END_RCPP
	}
}
