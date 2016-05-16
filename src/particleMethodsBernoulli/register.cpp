#include <Rcpp.h>
#include <internal.h>
#ifdef _MSC_VER
	#undef RcppExport
	#define RcppExport extern "C" __declspec(dllexport)
#endif
#include "importanceSampling.h"
#include "importanceResampling.h"
#include "importanceResamplingWithoutReplacement.h"
#include "importanceResamplingWithoutReplacementAuxiliary.h"
#include "importanceResamplingAuxiliary.h"
extern "C" const char* package_name = "particleMethodsBernoulli";
R_CallMethodDef callMethods[] = 
{
	{"importanceSampling", (DL_FUNC)&particleMethodsBernoulli::importanceSampling, 4},
	{"importanceResampling", (DL_FUNC)&particleMethodsBernoulli::importanceResampling, 4},
	{"importanceResamplingAuxiliary", (DL_FUNC)&particleMethodsBernoulli::importanceResamplingAuxiliary, 4},
	{"importanceResamplingWithoutReplacement", (DL_FUNC)&particleMethodsBernoulli::importanceResamplingWithoutReplacement, 4},
	{"importanceResamplingWithoutReplacementAuxiliary", (DL_FUNC)&particleMethodsBernoulli::importanceResamplingWithoutReplacementAuxiliary, 4},
	{NULL, NULL, 0}
};
RcppExport void R_init_particleMethodsBernoulli(DllInfo *info)
{
	std::vector<R_CallMethodDef> callMethodsVector;
	R_CallMethodDef* packageCallMethods = callMethods;
	while(packageCallMethods->name != NULL) packageCallMethods++;
	callMethodsVector.insert(callMethodsVector.begin(), callMethods, packageCallMethods);

	R_CallMethodDef* RcppStartCallMethods = Rcpp_get_call();
	R_CallMethodDef* RcppEndCallMethods = RcppStartCallMethods;
	while(RcppEndCallMethods->name != NULL) RcppEndCallMethods++;
	callMethodsVector.insert(callMethodsVector.end(), RcppStartCallMethods, RcppEndCallMethods);
	R_CallMethodDef blank = {NULL, NULL, 0};
	callMethodsVector.push_back(blank);

	R_registerRoutines(info, NULL, &(callMethodsVector[0]), NULL, NULL);
	init_Rcpp_cache();
}
