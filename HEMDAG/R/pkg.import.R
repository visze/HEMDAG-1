#' @import graph
#' @import RBGL
#' @import  precrec				   
#' @import  PerfMeas			   
#' @import  preprocessCore
#' @import methods
#' @importFrom utils read.table write.table

## Quiet concerns of R CMD check. 
## Avoid this warning: no visible binding for global variable
if(getRversion() >= "2.15.1"){
	utils::globalVariables("curvetypes");
}
