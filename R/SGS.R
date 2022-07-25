#' @useDynLib SGS, .registration = TRUE
NULL

.onLoad = function(lib, pkg) 
{
  # library.dynam("SGS", package = pkg, lib.loc = lib)
	if( R.version$arch == "x86_64" )
		MAX_NODES <- 64
	else
		MAX_NODES <- 32

  assign(".SGS.env", new.env(), envir=parent.env(environment()))
  assign("SGS.log.indent.tracker", 0, envir = .SGS.env)
}#.ONLOAD

.onUnload = function(lib) 
{
  # library.dynam.unload("SGS", libpath = lib)
}#.ONUNLOAD
