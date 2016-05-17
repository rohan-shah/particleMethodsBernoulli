#include "importanceResamplingWithoutReplacementAuxiliary.h"
#include "includeMPFR.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/random_number_generator.hpp>
#include <boost/range/algorithm/random_shuffle.hpp>
namespace particleMethodsBernoulli
{
	SEXP importanceResamplingWithoutReplacementAuxiliary(SEXP lowerBound_sexp, SEXP trueProbabilities_sexp, SEXP n_sexp, SEXP seed_sexp)
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
		//Initially we have two samples, corresponding to the first bernoulli being 0 or 1. Note that nBernoullis == 1 gives an error above, so we can assume that there are at least 2 bernoullis
		std::vector<int> samples, newSamples;
		std::vector<mpfr_class> sampleDensity, newSampleDensity;

		samples.push_back(0);
		samples.push_back(1);
		sampleDensity.push_back(1 - trueProbabilities[0]);
		sampleDensity.push_back(trueProbabilities[0]);
		std::vector<int> choicesUp, choicesDown;

		mpfr_class product = 1;
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
			newSampleDensity.clear();
			double complementaryTrueProb = 1 - trueProbabilities[bernoulliCounter];
			//We take every successor
			if(choicesUp.size() + choicesDown.size() <= (std::size_t)n)
			{
				for(std::size_t i = 0; i < choicesDown.size(); i++)
				{
					newSamples.push_back(samples[choicesDown[i]]);
					newSampleDensity.push_back(sampleDensity[choicesDown[i]] * complementaryTrueProb);
				}
				for(std::size_t i = 0; i < choicesUp.size(); i++)
				{
					newSamples.push_back(samples[choicesUp[i]]+1);
					newSampleDensity.push_back(sampleDensity[choicesUp[i]] * trueProbabilities[bernoulliCounter]);
				}
			}
			else
			{
				product *= (double) (choicesUp.size() + choicesDown.size()) /  (double) n;
				for(int i = 0; i < n; i++)
				{
					boost::random::bernoulli_distribution<> upOrDown((double)choicesUp.size() / ((double)choicesUp.size() + (double)choicesDown.size()));
					if(upOrDown(randomSource))
					{
						boost::random::uniform_int_distribution<> randomUp(0, choicesUp.size() - 1);
						int choice = randomUp(randomSource);
						newSamples.push_back(samples[choicesUp[choice]]+1);
						newSampleDensity.push_back(sampleDensity[choicesUp[choice]] * trueProbabilities[bernoulliCounter]);
						if(choice != (int)choicesUp.size() - 1) std::swap(choicesUp.back(), choicesUp[choice]);
						choicesUp.pop_back();
					}
					else
					{
						boost::random::uniform_int_distribution<> randomDown(0, choicesDown.size() - 1);
						int choice = randomDown(randomSource);
						newSamples.push_back(samples[choicesDown[choice]]);
						newSampleDensity.push_back(sampleDensity[choicesDown[choice]] * complementaryTrueProb);
						if(choice != (int)choicesDown.size() - 1) std::swap(choicesDown.back(), choicesDown[choice]);
						choicesDown.pop_back();
					}
				}
			}
			samples.swap(newSamples);
			sampleDensity.swap(newSampleDensity);
		}
		mpfr_class estimate = 0;
		for(std::size_t i = 0; i < samples.size(); i++)
		{
			estimate += sampleDensity[i];
		}
		estimate *= product;
		return Rcpp::List::create(Rcpp::Named("estimate") = estimate.convert_to<double>());
	END_RCPP
	}
}
