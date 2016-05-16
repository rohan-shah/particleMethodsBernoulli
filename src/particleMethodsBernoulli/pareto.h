#ifndef PARETO_SAMPLING_HEADER_GUARD
#define PARETO_SAMPLING_HEADER_GUARD
#include <vector>
#include <boost/random/mersenne_twister.hpp>
#include "includeMPFR.h"
namespace particleMethodsBernoulli
{
	struct paretoSamplingArgs
	{
	public:
		paretoSamplingArgs ()
		{}
		std::size_t n;
		std::vector<bool> deterministicInclusion;
		struct paretoStatistic
		{
			double statistic;
			int order;
			bool operator<(const paretoStatistic& other)
			{
				return statistic < other.statistic;
			}
		};
		std::vector<paretoStatistic> paretoStatistics;
		bool calculateInclusionProbabilities;
	};
	void pareto(paretoSamplingArgs& args, std::vector<int>& indices, std::vector<mpfr_class>& inclusionProbabilities, std::vector<mpfr_class>& weights, boost::mt19937& randomSource);
}
#endif