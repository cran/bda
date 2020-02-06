.onUnload <- function(libpath)
  library.dynam.unload("bda",  libpath)

.onAttach <- function(libname, pkgname)
    packageStartupMessage("bda v14 (Bin Wang, 2020)")

.bdaConnect <- NULL

