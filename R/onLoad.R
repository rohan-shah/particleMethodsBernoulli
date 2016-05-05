.onLoad <- function(libname, pkgname)
{
	library.dynam(package="particleMethodsBernoulli", chname="particleMethodsBernoulli", lib.loc = .libPaths())
}
