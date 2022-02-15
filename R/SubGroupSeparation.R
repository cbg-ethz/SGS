#' @useDynLib SubGroupSeparation, .registration = TRUE
NULL

.onLoad = function(lib, pkg) 
{
  # library.dynam("SubGroupSeparation", package = pkg, lib.loc = lib)
	if( R.version$arch == "x86_64" )
		MAX_NODES <- 64
	else
		MAX_NODES <- 32

  assign(".SubGroupSeparation.env", new.env(), envir=parent.env(environment()))
  assign("SubGroupSeparation.log.indent.tracker", 0, envir = .SubGroupSeparation.env)
}#.ONLOAD

.onUnload = function(lib) 
{
  # library.dynam.unload("SubGroupSeparation", libpath = lib)
}#.ONUNLOAD
