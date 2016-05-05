importanceSampling <- function(nBernoullis, lowerBound, probability, n, seed)
{
	start <- Sys.time()
	result <- .Call("importanceSampling", nBernoullis, lowerBound, probability, n, seed, PACKAGE="particleMethodsBernoulli")
	end <- Sys.time()
	return(new("importanceSamplingResult", estimate = result$estimate, estimateOfVariance = result$estimateOfVariance, call = match.call(), start = start, end = end))
}
