.onUnload <- function(libpath)
  library.dynam.unload("bda",  libpath)

.onAttach <- function(libname, pkgname)
    packageStartupMessage("bda v17 (Bin Wang, 2023)")

.bdaConnect <- NULL

