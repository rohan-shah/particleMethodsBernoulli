library(particleMethodsBernoulli)
exact <- 2.29224e-7
summarised <- data.frame()
source("./generateScenarios.R")
for(i in 1:nrow(scenarios))
{
	if(file.exists(scenarios[i, "file"]))
	{
		a <- load(scenarios[i, "file"])
		if(!("estimate" %in% a))
		{
			estimate <- mean(unlist(lapply(results, function(x) x@estimate)))
			var <- var(unlist(lapply(results, function(x) x@estimate)))
			times <- unlist(lapply(results, function(x) difftime(x@end, x@start, units = "sec")))
		}
		allEstimates <- unlist(lapply(results, function(x) x@estimate))
		summarised <- rbind(summarised, data.frame(estimate = estimate, variance = var, wnv = var*mean(times), wnrv = var*mean(times)/(exact^2), re = sqrt(var)/exact, meanTime = mean(times), mse = mean((allEstimates - exact)^2), stringsAsFactors = FALSE))
	}
	else
	{
		summarised <- rbind(summarised, data.frame(estimate = NA, variance = NA, wnv = NA, wnrv = NA, re = NA, meanTime = NA, mse = NA, stringsAsFactors = FALSE))
	}
}
summarised <- cbind(scenarios, summarised)
library(xtable)

summarised <- summarised[,c("sampleSize", "method", "mse", "re")]
rewriteMethods <- function(x)
{
	if(x == "WO-Replacement") return("WO-Replacement-$S_1$")
	else if(x == "WO-Replacement-S2") return("WO-Replacement-$S_{2, 1}$")
	return(x)
}
summarised$method <- sapply(summarised$method, rewriteMethods)

summarised100 <- subset(summarised, sampleSize == 100L)
summarised1000 <- subset(summarised, sampleSize == 1000L)
summarised100 <- summarised100[,-1]
summarised1000 <- summarised1000[,-1]

#Standardise ordering
summarised100 <- summarised100[match(c("IS", "WO-Replacement-$S_1$", "WO-Replacement-$S_{2, 1}$", "Bootstrap"), summarised100$method),]
summarised1000 <- summarised1000[match(c("IS", "WO-Replacement-$S_1$", "WO-Replacement-$S_{2, 1}$", "Bootstrap"), summarised1000$method),]

colnames(summarised100) <- colnames(summarised1000) <- c("Method", "MSE", "RE")
print(xtable(summarised100, display = c("s", "s", "e", "e"), label = "table:particleBernoulliUnequal_100", caption = "Simulation results for Bernoullis with unequal probabilities, with $\\sampleSize = 100$"), include.rownames=FALSE, math.style.exponents = TRUE, sanitize.text.function = identity)
print(xtable(summarised1000, display = c("s", "s", "e", "e"), label = "table:particleBernoulliUnequal_1000", caption = "Simulation results for Bernoullis with unequal probabilities, with $\\sampleSize = 1000$"), include.rownames=FALSE, math.style.exponents = TRUE, sanitize.text.function = identity)

