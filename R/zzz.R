.onUnload <- function(libpath)
  library.dynam.unload("bda",  libpath)

.onAttach <- function(libname, pkgname)
    packageStartupMessage("bda v15 (Bin Wang, 2020)")

.bdaConnect <- NULL

