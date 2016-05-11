importanceSampling <- function(nBernoullis, lowerBound, probability, probabilities, n, seed)
{
	if(missing(probabilities))
	{
		if(missing(probability) || missing(nBernoullis))
		{
			stop("If input probabilities is missing, then inputs nBernoullis and probability must be given")
		}
		probabilities <- rep(probability, nBernoullis)
	}
	else
	{
		if(!missing(probability) || !missing(nBernoullis))
		{
			stop("If input probabilities is used, then inputs nBernoullis and probability must be missing")
		}
	}
	start <- Sys.time()
	result <- .Call("importanceSampling", lowerBound, probabilities, n, seed, PACKAGE="particleMethodsBernoulli")
	end <- Sys.time()
	return(new("importanceSamplingResult", estimate = result$estimate, estimateOfVariance = result$estimateOfVariance, call = match.call(), start = start, end = end))
}
