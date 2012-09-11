.onLoad <- function(lib, pkg){
  packageStartupMessage("bda 2.0.11-11 Copyright B. Wang 2011-2012.")
  assign('.bdaConnect',NULL,pos=.GlobalEnv) 
}

.onUnload <- function(libpath)
  library.dynam.unload("bda",  libpath)

.bdaConnect <- NULL

