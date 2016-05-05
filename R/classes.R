setClass("importanceSamplingResult", slots = list(estimate = "numeric", call = "call", estimateOfVariance = "numeric", start = "POSIXct", end = "POSIXct"))
setClass("importanceResamplingResult", slots = list(estimate = "numeric", call = "call", start = "POSIXct", end = "POSIXct"))
