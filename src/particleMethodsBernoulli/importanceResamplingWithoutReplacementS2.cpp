#include "importanceResamplingWithoutReplacementS2.h"
#include "includeMPFR.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/random_number_generator.hpp>
#include <boost/range/algorithm/random_shuffle.hpp>
#include "sampford.h"
namespace particleMethodsBernoulli
{

	int countBits(int i)
	{
		i = i - ((i >> 1) & 0x55555555);
		i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
		return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
	}
	SEXP importanceResamplingWithoutReplacementS2(SEXP lowerBound_sexp, SEXP trueProbabilities_sexp, SEXP n_sexp, SEXP seed_sexp, SEXP k_sexp)
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
		int k;
		try
		{
			k = Rcpp::as<int>(k_sexp);
		}
		catch(...)
		{
			throw std::runtime_error("Input k muts be an integer");
		}
		if(k > 32 || k < 1)
		{
			throw std::runtime_error("Input k must be between 1 and 32 (inclusive)");
		}
		
		boost::mt19937 randomSource;
		randomSource.seed(seed);
		boost::random_number_generator<boost::mt19937> generator(randomSource);

		double newProbability = (double)lowerBound / (double)nBernoullis;
		mpfr_class newProbability_mpfr = newProbability, compNewProbability_mpfr = 1 - newProbability;
		std::vector<mpfr_class> powers(nBernoullis+1);
		for(int i = 0; i < nBernoullis+1; i++)
		{
			powers[i] = boost::multiprecision::pow(newProbability_mpfr, i) * boost::multiprecision::pow(compNewProbability_mpfr, nBernoullis - i);
		}
		//Initially we have two samples, corresponding to the first bernoulli being 0 or 1. Note that nBernoullis == 1 gives an error above, so we can assume that there are at least 2 bernoullis
		std::vector<int> samples, newSamples;
		std::vector<mpfr_class> sampleDensityOnWeight, newSampleDensityOnWeight, sampfordWeights, sampfordWeights2, newSampfordWeights;
		std::vector<unsigned int> bits, newBits;

		samples.push_back(0);
		samples.push_back(1);
		sampleDensityOnWeight.push_back(1-trueProbabilities[0]);
		sampleDensityOnWeight.push_back(trueProbabilities[0]);
		sampfordWeights.push_back(1);
		sampfordWeights.push_back(1);
		std::vector<int> choicesUp, choicesDown;

		bits.push_back(0U);
		bits.push_back(1U << (k-1));

		sampling::sampfordFromParetoNaiveArgs sampfordArgs;
		sampfordArgs.n = n;

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
			newBits.clear();
			newSamples.clear();
			newSampleDensityOnWeight.clear();
			newSampfordWeights.clear();
			double complementaryTrueProbPrevious = std::numeric_limits<double>::quiet_NaN(), trueProbPrevious = std::numeric_limits<double>::quiet_NaN();
			mpfr_class trueProb = trueProbabilities[bernoulliCounter];
			mpfr_class complementaryTrueProb = 1 - trueProb;

			if(bernoulliCounter > k-1)
			{
				trueProbPrevious = trueProbabilities[bernoulliCounter - k];
				complementaryTrueProbPrevious = 1 - trueProbPrevious;
			}
			//We take every successor
			if(choicesUp.size() + choicesDown.size() <= (std::size_t)n)
			{
				for(std::size_t i = 0; i < choicesDown.size(); i++)
				{
					newSamples.push_back(samples[choicesDown[i]]);
					newSampleDensityOnWeight.push_back(sampleDensityOnWeight[choicesDown[i]] * complementaryTrueProb);
					newBits.push_back(bits[choicesDown[i]] >> 1);
					if(bernoulliCounter <= k-1)
					{
						newSampfordWeights.push_back(sampfordWeights[choicesDown[i]]);
					}
					else
					{
						if(bits[choicesDown[i]] & 1U)
						{
							newSampfordWeights.push_back(sampfordWeights[choicesDown[i]]*trueProbPrevious);
						}
						else
						{
							newSampfordWeights.push_back(sampfordWeights[choicesDown[i]]*complementaryTrueProbPrevious);
						}
					}
				}
				for(std::size_t i = 0; i < choicesUp.size(); i++)
				{
					newSamples.push_back(samples[choicesUp[i]]+1);
					newSampleDensityOnWeight.push_back(sampleDensityOnWeight[choicesUp[i]] * trueProb);
					newBits.push_back((bits[choicesUp[i]] >> 1) + (1U << (k-1)));
					if(bernoulliCounter <= k-1)
					{
						newSampfordWeights.push_back(sampfordWeights[choicesUp[i]]);
					}
					else
					{
						if(bits[choicesUp[i]] & 1U)
						{
							newSampfordWeights.push_back(sampfordWeights[choicesUp[i]]*trueProb);
						}
						else
						{
							newSampfordWeights.push_back(sampfordWeights[choicesUp[i]]*complementaryTrueProb);
						}
					}
				}
				newSampfordWeights.swap(sampfordWeights);
			}
			else
			{
				for(std::size_t i = 0; i < choicesDown.size(); i++)
				{
					newSampfordWeights.push_back(sampfordWeights[choicesDown[i]] * powers[countBits(bits[choicesDown[i]])]);
				}
				for(std::size_t i = 0; i < choicesUp.size(); i++)
				{
					newSampfordWeights.push_back(sampfordWeights[choicesUp[i]] * powers[countBits(bits[choicesUp[i]])+1]);
				}
				sampfordArgs.weights.swap(newSampfordWeights);
				sampling::sampfordFromParetoNaive(sampfordArgs, randomSource);
				sampfordArgs.weights.swap(newSampfordWeights);

				sampfordWeights2.clear();
				for(std::vector<int>::iterator j = sampfordArgs.indices.begin(); j != sampfordArgs.indices.end(); j++)
				{
					if(*j < (int)choicesDown.size())
					{
						newSamples.push_back(samples[choicesDown[*j]]);
						newSampleDensityOnWeight.push_back(sampleDensityOnWeight[choicesDown[*j]] * complementaryTrueProb / sampfordArgs.inclusionProbabilities[*j]);
						if(bernoulliCounter <= k-1)
						{
							sampfordWeights2.push_back(sampfordWeights[choicesDown[*j]] / sampfordArgs.inclusionProbabilities[*j]);
						}
						else
						{
							if(bits[choicesDown[*j]] & 1U)
							{
								sampfordWeights2.push_back(sampfordWeights[choicesDown[*j]]*trueProbPrevious / sampfordArgs.inclusionProbabilities[*j]);
							}
							else
							{
								sampfordWeights2.push_back(sampfordWeights[choicesDown[*j]]*complementaryTrueProbPrevious / sampfordArgs.inclusionProbabilities[*j]);
							}
						}
						newBits.push_back((bits[choicesDown[*j]] >> 1));
					}
					else
					{
						int choiceUpIndex = choicesUp[*j - choicesDown.size()];
						newSamples.push_back(samples[choiceUpIndex]+1);
						newSampleDensityOnWeight.push_back(sampleDensityOnWeight[choiceUpIndex] * trueProb / sampfordArgs.inclusionProbabilities[*j]);
						if(bernoulliCounter <= k-1)
						{
							sampfordWeights2.push_back(sampfordWeights[choiceUpIndex] / sampfordArgs.inclusionProbabilities[*j]);
						}
						else
						{
							if(bits[choiceUpIndex] & 1U)
							{
								sampfordWeights2.push_back(sampfordWeights[choiceUpIndex]*trueProbPrevious / sampfordArgs.inclusionProbabilities[*j]);
							}
							else
							{
								sampfordWeights2.push_back(sampfordWeights[choiceUpIndex]*complementaryTrueProbPrevious / sampfordArgs.inclusionProbabilities[*j]);
							}
						}
						newBits.push_back((bits[choiceUpIndex] >> 1) + (1U << (k-1)));
					}
				}
				sampfordWeights.swap(sampfordWeights2);
			}
			samples.swap(newSamples);
			sampleDensityOnWeight.swap(newSampleDensityOnWeight);
			bits.swap(newBits);
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
