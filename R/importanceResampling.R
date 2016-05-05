importanceResampling <- function(nBernoullis, lowerBound, probability, n, seed)
{
	start <- Sys.time()
	result <- .Call("importanceResampling", nBernoullis, lowerBound, probability, n, seed, PACKAGE="particleMethodsBernoulli")
	end <- Sys.time()
	return(new("importanceResamplingResult", estimate = result$estimate, call = match.call(), start = start, end = end))
}
