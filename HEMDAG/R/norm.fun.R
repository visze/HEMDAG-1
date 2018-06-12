##**********##
## MAX-NORM ##
##**********##

#' @title Max normalization
#' @description Function to normalize the scores of a flat scores matrix per class
#' @details The scores of each class are normalized by dividing the score values for the maximum score of that class.
#' If the max score of a class is zero, no normalization is needed, otherwise \code{NaN} value will be printed as results of 0 out of 0 division.
#' @param S matrix with the raw non normalized scores. Rows are examples and columns are classes
#' @return A score matrix with the same dimensions of \code{S}, but with scores max/normalized separately for each class
#' @export
#' @examples
#' data(scores);
#' maxnorm <- normalize.max(S);
normalize.max <- function(S){
	classes <- colnames(S);
	maximum <- apply(S,2,max);
	for(class in classes){
		if(maximum[class] != 0){
			S[,class] <- S[,class]/maximum[class];
		}
	}
	return(S);
}

##********************##
## FLAT NORMALIZATION ##
##********************##

#' @title Flat scores normalization
#' @description High level functions to normalize a flat scores matrix w.r.t. max normalization (MaxNorm) or quantile normalization (Qnorm) 
#' @details To apply the quantile normalization the \pkg{preprocessCore} library is uded.
#' @param norm.type can be one of the following two values:
#' \itemize{
#' \item MaxNorm: each score is divided w.r.t. the max of each class;
#' \item Qnorm: a quantile normalization is applied. Library preprocessCore is used.
#' }
#' @param flat.file name of the flat scores matrix (without rda extension)
#' @param dag.file name of the graph that represents the hierarchy of the classes
#' @param flat.dir relative path to folder where flat normalized scores matrix is stored 
#' @param dag.dir relative path to folder where graph is stored
#' @param flat.norm.dir the directory where the normalized flat scores matrix must be stored
#' @return the matrix of the scores flat normalized w.r.t. MaxNorm or Qnorm
#' @export
#' @examples
#' data(scores);
#' data(graph);
#' if (!dir.exists("data")){
#' 	dir.create("data");
#' }
#' if (!dir.exists("results")){
#' 	dir.create("results");
#' }
#' save(S,file="data/scores.rda");
#' save(g,file="data/graph.rda");
#' flat.dir <- dag.dir <- "data/";
#' flat.norm.dir <- "results/";
#' flat.file <- "scores";
#' dag.file <- "graph";
#' norm.types <- c("MaxNorm","Qnorm");
#' for(norm.type in norm.types){
#' 	Do.FLAT.scores.normalization(norm.type=norm.type, flat.file=flat.file, 
#' 	dag.file=dag.file, flat.dir=flat.dir, dag.dir=dag.dir, 
#' 	flat.norm.dir=flat.norm.dir);
#' }
Do.FLAT.scores.normalization <- function(norm.type="MaxNorm", flat.file=flat.file, dag.file=dag.file, 
	flat.dir=flat.dir, dag.dir=dag.dir, flat.norm.dir=flat.norm.dir){
	
	## loading dag
	dag.path <- paste0(dag.dir, dag.file,".rda");
	g <- get(load(dag.path));
	
	## root node
	root <- root.node(g);

	## loading flat scores matrix 
	flat.path <- paste0(flat.dir, flat.file,".rda");
	S <- get(load(flat.path));
		
	## removing root node from flat matrix if it exists
	if(root %in% colnames(S)){
		S <- S[,-which(colnames(S)==root)];
	}

	## normalization
	if(norm.type=="MaxNorm"){
		## Max Normalization
		S <- normalize.max(S);		
	}else if(norm.type=="Qnorm"){
		## Quantile Normalization 
		## NOTE: normalize.quantiles function returns a unnamed matrix. colnames are essential for hier.corr..
		S.norm <- normalize.quantiles(S);
		dimnames(S.norm) <- list(rownames(S),colnames(S));
		S <- S.norm;
		rm(S.norm);	
	}else{
		stop("The chosen normalization method is not among those available or it was misspelled");
	}

	## Storing results
	save(S, file=paste0(flat.norm.dir, norm.type, ".", flat.file, ".rda"));
}
