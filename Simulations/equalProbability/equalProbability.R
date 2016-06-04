library(particleMethodsBernoulli)
source("./generateScenarios.R")
exact <- 1.35138e-10

SCENARIO_INDEX <- as.integer(Sys.getenv("SCENARIO_INDEX"))
replications <- scenarios[SCENARIO_INDEX, "replications"]
sampleSize <- scenarios[SCENARIO_INDEX, "sampleSize"]
method <- scenarios[SCENARIO_INDEX, "method"]
file <- scenarios[SCENARIO_INDEX, "file"]
if(method == "IS")
{
	#Importance sampling simulation
	if(file.exists(file))
	{
		load(file)
	}
	else results <- list()
	if(length(results) != replications)
	{
		for(i in (length(results)+1):replications)
		{
			results[[i]] <- importanceSampling(nBernoullis = 100, lowerBound = 80, probability = 0.5, n = sampleSize, seed = i)
			if(i %% 1000 == 0) save(results, file = file)
		}
	}
	estimate <- mean(unlist(lapply(results, function(x) x@estimate)))
	var <- var(unlist(lapply(results, function(x) x@estimate)))
	times <- unlist(lapply(results, function(x) difftime(x@end, x@start, units = "sec")))
	save(results, estimate, var, times, file = file)
} else if(method == "WO-Replacement")
{
	#Without replacement simulation
	if(file.exists(file))
	{
		load(file)
	}
	else results <- list()
	if(length(results) != replications)
	{
		for(i in (length(results)+1):replications)
		{
			results[[i]] <- importanceResamplingWithoutReplacement(nBernoullis = 100, lowerBound = 80, probability = 0.5, n = sampleSize, seed = i)
			if(i %% 1000 == 0) save(results, file = file)
		}
	}
	estimate <- mean(unlist(lapply(results, function(x) x@estimate)))
	var <- var(unlist(lapply(results, function(x) x@estimate)))
	times <- unlist(lapply(results, function(x) difftime(x@end, x@start, units = "sec")))
	save(results, estimate, var, times, file = file)
} else if(method == "WO-Replacement-S2")
{
	#Without replacement simulation
	if(file.exists(file))
	{
		load(file)
	}
	else results <- list()
	if(length(results) != replications)
	{
		for(i in (length(results)+1):replications)
		{
			results[[i]] <- importanceResamplingWithoutReplacementS2(nBernoullis = 100, lowerBound = 80, probability = 0.5, n = sampleSize, seed = i, k = 1)
			if(i %% 1000 == 0) save(results, file = file)
		}
	}
	estimate <- mean(unlist(lapply(results, function(x) x@estimate)))
	var <- var(unlist(lapply(results, function(x) x@estimate)))
	times <- unlist(lapply(results, function(x) difftime(x@end, x@start, units = "sec")))
	save(results, estimate, var, times, file = file)

} else if(method == "Bootstrap")
{
	#With replacement simulation
	if(file.exists(file))
	{
		load(file)
	}
	else results <- list()
	if(length(results) != replications)
	{
		for(i in (length(results)+1):replications)
		{
			results[[i]] <- importanceResampling(nBernoullis = 100, lowerBound = 80, probability = 0.5, n = sampleSize, seed = i)
			if(i %% 1000 == 0) save(results, file = file)
		}
	}
	estimate <- mean(unlist(lapply(results, function(x) x@estimate)))
	var <- var(unlist(lapply(results, function(x) x@estimate)))
	times <- unlist(lapply(results, function(x) difftime(x@end, x@start, units = "sec")))
	save(results, estimate, var, times, file = file)
} else if(method == "Bootstrap-Auxiliary")
{
	if(file.exists(file))
	{
		load(file)
	}
	else results <- list()
	if(length(results) != replications)
	{
		for(i in (length(results)+1):replications)
		{
			results[[i]] <- importanceResamplingAuxiliary(nBernoullis = 100, lowerBound = 80, probability = 0.5, n = sampleSize, seed = i)
			if(i %% 1000 == 0) save(results, file = file)
		}
	}
	estimate <- mean(unlist(lapply(results, function(x) x@estimate)))
	var <- var(unlist(lapply(results, function(x) x@estimate)))
	times <- unlist(lapply(results, function(x) difftime(x@end, x@start, units = "sec")))
	save(results, estimate, var, times, file = file)
} else
{
	stop("Unrecognized method")
}
