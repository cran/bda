.onLoad <- function(lib, pkg){
  packageStartupMessage("BDA 1.0 Copyright B. Wang 2011.\n")
  assign('.bdaConnect',NULL,pos=.GlobalEnv) 
}

.onUnload <- function(libpath)
    library.dynam.unload("bda",  libpath)

.bdaConnect <- NULL

