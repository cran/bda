.onUnload <- function(libpath)
  library.dynam.unload("bda",  libpath)

.onAttach <- function(libname, pkgname)
    packageStartupMessage(paste("bda -", as.character(packageVersion("bda"))))

.bdaConnect <- NULL

