## Attach short description of HEMDAG library when load it
.onAttach <- function(libname=.libPaths(), pkgname="HEMDAG"){
	packageStartupMessage("HEMDAG: Hierarchical Ensemble Methods for DAG-structured taxonomies\nPlease cite HEMDAG if you use it: see citation('HEMDAG') for details\n");
}
