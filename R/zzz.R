.onLoad <- function(lib, pkg){
  packageStartupMessage("bda 1.2 Copyright B. Wang 2011-2012.\n")
  assign('.bdaConnect',NULL,pos=.GlobalEnv) 
}

.onUnload <- function(libpath)
    library.dynam.unload("bda",  libpath)

.bdaConnect <- NULL

