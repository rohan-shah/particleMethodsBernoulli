source("./generateScenarios.R")
if(!file.exists("./results")) dir.create("results")
if(length(list.files(path="./results", pattern=".RData.lock")) > 0)
{
	stop("Please remove lock directories")
}
for(i in 1:nrow(scenarios))
{
	submit <- FALSE
	resultFile <- scenarios[i, "file"]
	if(!file.exists(resultFile))
	{
		submit <- TRUE
	}
	else
	{
		a <- load(resultFile)
		if(any(!(c("results", "estimate", "var", "times") %in% a))) submit <- TRUE
		else if(length(results) != scenarios[i, "replications"]) submit <- TRUE
	}
	if(submit)
	{
		system2(command = "qsub", args = "equalProbability", env = paste0("SCENARIO_INDEX=", i), wait=TRUE)
	}
}
