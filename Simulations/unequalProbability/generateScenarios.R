methods <- c("IS", "WO-Replacement", "Bootstrap", "WO-Replacement-S2")
scenarios <- expand.grid(method = methods, sampleSize = c(100L, 1000L), replications = 50000L, stringsAsFactors=FALSE)
scenarios$file <- sapply(1:nrow(scenarios), function(x) paste0("./results/", scenarios[x, "method"], "-", scenarios[x, "sampleSize"], ".RData")) 
