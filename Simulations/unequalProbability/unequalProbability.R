library(particleMethodsBernoulli)
source("./generateScenarios.R")
exact <- 2.29224e-7

SCENARIO_INDEX <- as.integer(Sys.getenv("SCENARIO_INDEX"))
replications <- scenarios[SCENARIO_INDEX, "replications"]
sampleSize <- scenarios[SCENARIO_INDEX, "sampleSize"]
method <- scenarios[SCENARIO_INDEX, "method"]
file <- scenarios[SCENARIO_INDEX, "file"]
tmpFile <- paste0(file, ".tmp")
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
			results[[i]] <- importanceSampling(lowerBound = 95, probabilities = c(rep(0.001, 5), rep(0.9, 95)), n = sampleSize, seed = i)
			if(i %% 1000 == 0)
			{
				save(results, file = tmpFile)
				file.rename(tmpFile, file)
			}
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
			results[[i]] <- importanceResamplingWithoutReplacement(lowerBound = 95, probabilities = c(rep(0.001, 5), rep(0.9, 95)), n = sampleSize, seed = i)
			if(i %% 1000 == 0) 
			{
				save(results, file = tmpFile)
				file.rename(tmpFile, file)
			}
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
			results[[i]] <- importanceResamplingWithoutReplacementS2(lowerBound = 95, probabilities = c(rep(0.001, 5), rep(0.9, 95)), n = sampleSize, seed = i, k = 1)
			if(i %% 1000 == 0)
			{
				save(results, file = tmpFile)
				file.rename(tmpFile, file)
			}
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
			results[[i]] <- importanceResampling(lowerBound = 95, probabilities = c(rep(0.001, 5), rep(0.9, 95)), n = sampleSize, seed = i)
			if(i %% 1000 == 0)
			{
				save(results, file = tmpFile)
				file.rename(tmpFile, file)
			}
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
			results <- importanceResamplingAuxiliary(lowerBound = 95, probabilities = c(rep(0.001, 5), rep(0.9, 95)), n = sampleSize, seed = i)
			if(i %% 1000 == 0)
			{
				save(results, file = tmpFile)
				file.rename(tmpFile, file)
			}
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
