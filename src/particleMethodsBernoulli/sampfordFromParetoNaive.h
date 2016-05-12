#ifndef SAMPFORD_FROM_PARETO_HEADER_GUARD
#define SAMPFORD_FROM_PARETO_HEADER_GUARD
#include "pareto.h"
namespace particleMethodsBernoulli
{
	struct sampfordFromParetoNaiveArgs
	{
	public:
		paretoSamplingArgs paretoArgs;
		sampfordFromParetoNaiveArgs()
		{}
		std::size_t n;
	};
	void sampfordFromParetoNaive(sampfordFromParetoNaiveArgs& args, std::vector<int>& indices, std::vector<mpfr_class>& inclusionProbabilities, std::vector<mpfr_class>& weights, boost::mt19937& randomSource);
}
#endif
