##*********##
## HTD-DAG ##
##*********##
#' @name HTD-DAG
#' @title HTD-DAG
#' @description Implementation of a top-down procedure to correct the scores of the hierarchy according to the 
#' constraints that the score of a node cannot be greater than a score of its parents.
#' @details The HTD-DAG algorithm modifies the flat scores according to the hierarchy of a DAG through a unique run across
#' the nodes of the graph. For a given example \eqn{x \in X}, the flat predictions \eqn{f(x) = \hat{y}} are hierarchically corrected to
#' \eqn{\bar{y}}, by per-level visiting the nodes of the DAG from top to bottom according to the following simple rule:
#' \deqn{
#' \bar{y}_i := \left\{
#'    \begin{array}{lll}
#'      \hat{y}_i  & {\rm if} \quad i \in root(G) \\
#'      \min_{j \in par(i)} \bar{y}_j & {\rm if} \quad \min_{j \in par(i)} \bar{y}_j < \hat{y}_i \\
#'      \hat{y}_i & {\rm otherwise}
#'    \end{array}
#'   \right.
#' }
#' The node levels correspond to their maximum path length from the root.
#' @seealso \code{\link{graph.levels}}, \code{\link{hierarchical.checkers}}
#' @param S a named flat scores matrix with examples on rows and classes on columns
#' @param g a graph of class \code{graphNEL}. It represents the hierarchy of the classes
#' @param root name of the class that it is the top-level (root) of the hierarchy (\code{def:00})
#' @return a matrix with the scores of the classes corrected according to the HTD-DAG algorithm.
#' @export
#' @examples
#' data(graph);
#' data(scores);
#' root <- root.node(g);
#' S.htd <- htd(S,g,root);
htd <- function(S, g, root="00"){
	levels <- graph.levels(g,root);

	# a dummy root is added
	if(!(root %in% colnames(S))){
		max.score <- max(S);
		z <- rep(max.score,nrow(S));
		S <- cbind(z,S);
		colnames(S)[1] <- root;
	}
	## check consistency between nodes of g and classes of S
	class.check <- ncol(S)!=numNodes(g);
	if(class.check)
		stop("HTD: The number of nodes of the graph and the number of classes of the scores matrix does not match", call.=FALSE);
		
	# nodes are scanned from top to bottom: a list par.tod with the parents for each node (ordered from top to bottom) is obtained	
	par.tod <- get.parents.top.down(g,levels,root)
	for(i in 1:length(par.tod)){  								    
		child <- S[,names(par.tod[i])]; 				    
		parents <- as.matrix(S[,par.tod[[i]]]);	 			    
		# colnames(parents) <- par.tod[[i]];
		# Note: the version with an apply and an ifelse statement is slower ...						    
		for(j in 1:length(child)){ 							    
			x <- min(parents[j,]);								    
			if(x < child[j])									    
				child[j] <- x;									    
		}  													    
		S[,names(par.tod[i])] <- child;						    
	}													    
	# the dummy root is removed
	S <- S[,-which(colnames(S)==root)];
	return(S);
}

##************##
## DO HTD-DAG ##
##************##
#' @title HTD-DAG vanilla
#' @description High level function to correct the computed scores in a hierarchy according to the HTD-DAG algorithm
#' @details The function check if the number of classes between the flat scores matrix and the annotations matrix mismatched
#' (e.g. some terms might be removed from the flat scores matrix before applying an hierarchy-unaware algorithm). 
#' If so, the number of classes of the annotations matrix is narrowed to the number of terms of the flat scores matrix and
#' the corresponding subgraph is computed as well.
#' @param norm boolean value: 
#' \itemize{
#' \item \code{TRUE} (def.): the flat scores matrix has been already normalized in according to a normalization method;	
#' \item \code{FALSE}: the flat scores matrix has not been normalized yet. See the parameter \code{norm.type} for which normalization can be applied.
#' }
#' @param norm.type can be one of the following three values:
#'  \enumerate{
#'  \item \code{NULL} (def.): set \code{norm.type} to \code{NULL} if and only if the parameter \code{norm} is set to \code{TRUE};
#'  \item \code{MaxNorm}: each score is divided for the maximum of each class;
#'  \item \code{Qnorm}: quantile normalization. \pkg{preprocessCore} package is used. 
#'  }
#' @param folds number of folds of the cross validation on which computing the performance metrics averaged across folds (\code{def. 5}).
#' If \code{folds=NULL}, the performance metrics are computed one-shot, otherwise the performance metrics are averaged across folds.
#' @param seed initialization seed for the random generator to create folds (\code{def. 23}). If \code{NULL} folds are generated without seed 
#' initialization. 
#' @param n.round number of rounding digits to be applied to the hierarchical scores matrix (\code{def. 3}). It is used for choosing 
#' the best threshold on the basis of the best F-measure
#' @param f.criterion character. Type of F-measure to be used to select the best F-measure. Two possibilities:
#' \enumerate{
#' \item \code{F} (def.): corresponds to the harmonic mean between the average precision and recall
#' \item \code{avF}: corresponds to the per-example \code{F-score} averaged across all the examples
#' }
#' @param rec.levels a vector with the desired recall levels (\code{def:} \code{from:0.1}, \code{to:0.9}, \code{by:0.1}) to compute the 
#' the Precision at fixed Recall level (PXR)
#' @param flat.file name of the file containing the flat scores matrix to be normalized or already normalized (without rda extension)
#' @param ann.file name of the file containing the the label matrix of the examples (without rda extension)
#' @param dag.file name of the file containing the graph that represents the hierarchy of the classes (without rda extension)
#' @param flat.dir relative path where flat scores matrix is stored
#' @param ann.dir relative path where annotation matrix is stored
#' @param dag.dir relative path where graph is stored
#' @param hierScore.dir relative path where the hierarchical scores matrix must be stored
#' @param perf.dir relative path where all the performance measures must be stored
#' @return Two \code{rda} files stored in the respective output directories:
#' \enumerate{
#' 	\item \code{Hierarchical Scores Results}: a matrix with examples on rows and classes on columns representing the computed hierarchical scores 
#' 	for each example and for each considered class. It is stored in the \code{hierScore.dir} directory.
#' 	\item \code{Performance Measures}: \emph{flat} and \emph{hierarchical} performace results:
#' 	\enumerate{
#' 		\item AUPRC results computed though \code{AUPRC.single.over.classes} (\code{\link{AUPRC}});
#'		\item AUROC results computed through \code{AUROC.single.over.classes} (\code{\link{AUROC}}); 
#' 		\item PXR results computed though \code{PXR.at.multiple.recall.levels.over.classes} (\code{\link{PXR}});
#' 		\item FMM results computed though \code{compute.Fmeasure.multilabel} (\code{\link{FMM}}); 
#' }}
#' It is stored in the \code{perf.dir} directory.
#' @seealso \code{\link{HTD-DAG}}
#' @export
#' @examples
#' data(graph);
#' data(scores);
#' data(labels);
#' if (!dir.exists("data")){
#' 	dir.create("data");
#' }
#' if (!dir.exists("results")){
#' 	dir.create("results");
#' }
#' save(g,file="data/graph.rda");
#' save(L,file="data/labels.rda");
#' save(S,file="data/scores.rda");
#' dag.dir <- flat.dir <- ann.dir <- "data/";
#' hierScore.dir <- perf.dir <- "results/";
#' rec.levels <- seq(from=0.1, to=1, by=0.1);
#' dag.file <- "graph";
#' flat.file <- "scores";
#' ann.file <- "labels";
#' Do.HTD(norm=FALSE, norm.type="MaxNorm", folds=5, seed=23, n.round=3, 
#' f.criterion="F", rec.levels=rec.levels, flat.file=flat.file, ann.file=ann.file, 
#' dag.file=dag.file, flat.dir=flat.dir, ann.dir=ann.dir, dag.dir=dag.dir, 
#' hierScore.dir=hierScore.dir, perf.dir=perf.dir);
Do.HTD <- function(norm=TRUE, norm.type=NULL, folds=5, seed=23, n.round=3, f.criterion ="F", 
	rec.levels=seq(from=0.1, to=1, by=0.1), flat.file=flat.file, ann.file=ann.file, dag.file=dag.file, 
	flat.dir=flat.dir, ann.dir=ann.dir, dag.dir=dag.dir, hierScore.dir=hierScore.dir, perf.dir=perf.dir){
	
	## Setting Check
	if(norm==FALSE && is.null(norm.type))
		stop("HTD: If norm is set to FALSE, you need also to specify a normalization method among those available", call.=FALSE)
	if(norm==TRUE && !is.null(norm.type))
		warning("HTD: If norm is set to TRUE, the input flat matrix is already normalized.", 
			paste0(" Set norm.type to NULL and not to '", norm.type, "' to avoid this warning message"), call.=FALSE);
				
	## loading dag
	dag.path <- paste0(dag.dir, dag.file,".rda");
	g <- get(load(dag.path));
	## compute root node 
	root <- root.node(g);

	## loading annotation matrix
	ann.path <- paste0(ann.dir, ann.file,".rda");
	ann <- get(load(ann.path));
	## removing root node from annotation table if it exists
	if(root %in% colnames(ann))
		ann <- ann[,-which(colnames(ann)==root)];

	## loading flat matrix
	flat.path <- paste0(flat.dir, flat.file,".rda");
	if(norm){
		S <- get(load(flat.path));
		## removing root node from flat norm matrix if it exists
		if(root %in% colnames(S))
			S <- S[,-which(colnames(S)==root)];
	}else{
		S <- get(load(flat.path));
		S <- scores.normalization(norm.type=norm.type, S);
		cat(norm.type, "NORMALIZATION: DONE", "\n");
		## removing root node from flat norm matrix if it exists
		if(root %in% colnames(S))
			S <- S[,-which(colnames(S)==root)];
	}

	## check if |flat matrix classes| = |annotation matrix classes| 
	## if not the classes of annotation matrix are shrunk to those of flat matrix
	class.check <- ncol(S)!=ncol(ann);
	if(class.check){
		ann <- ann[,colnames(S)];
		nd <- colnames(S);
		g <- do.subgraph(nd, g, edgemode="directed");
	}

	## Compute FLAT PRC, AUC, PXR (average and per class) and FMM (average and per-example) one-shoot or cross-validated 
	PRC.flat <- AUPRC.single.over.classes(ann, S, folds=folds, seed=seed);
	AUC.flat <- AUROC.single.over.classes(ann, S, folds=folds, seed=seed);
	PXR.flat <- PXR.at.multiple.recall.levels.over.classes(ann, S, rec.levels=rec.levels, folds=folds, seed=seed);
	FMM.flat <- compute.Fmeasure.multilabel(ann, S, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, 
		b.per.example=TRUE, folds=folds, seed=seed);

	## Hierarchical Top Down Correction
	S <- htd(S, g, root);
	
	## Compute HIER PRC, AUC, PXR (average and per class) and FMM (average and per-example) one-shoot or cross-validated 
	PRC.hier <- AUPRC.single.over.classes(ann, S, folds=folds, seed=seed);
	AUC.hier <- AUROC.single.over.classes(ann, S, folds=folds, seed=seed);
	PXR.hier <- PXR.at.multiple.recall.levels.over.classes(ann, S, rec.levels=rec.levels, folds=folds, seed=seed);
	FMM.hier <- compute.Fmeasure.multilabel(ann, S, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, 
		b.per.example=TRUE, folds=folds, seed=seed);
	S.hier <- S;
	rm(S);

	## Storing Results 
	if(norm){
		save(S.hier, file=paste0(hierScore.dir, flat.file, ".hierScores.htd.rda"), compress=TRUE);
		save(PRC.flat, PRC.hier, AUC.flat, AUC.hier, PXR.flat, PXR.hier, FMM.flat, FMM.hier,
			file=paste0(perf.dir, "PerfMeas.", flat.file, ".hierScores.htd.rda"), compress=TRUE);
	}else{
		save(S.hier, file=paste0(hierScore.dir, norm.type,".", flat.file, ".hierScores.htd.rda"), compress=TRUE);	
		save(PRC.flat, PRC.hier, AUC.flat, AUC.hier, PXR.flat, PXR.hier, FMM.flat, FMM.hier,
			file=paste0(perf.dir, "PerfMeas.", norm.type,".", flat.file, ".hierScores.htd.rda"), compress=TRUE);
	}
}

#' @title HTD-DAG holdout
#' @description High level function to correct the computed scores in a hierarchy according to the HTD-DAG algorithm applying a 
#' classical holdout procedure
#' @details The function check if the number of classes between the flat scores matrix and the annotations matrix mismatched
#' (e.g. some terms might be removed from the flat scores matrix before applying an hierarchy-unaware algorithm). 
#' If so, the number of classes of the annotations matrix is narrowed to the number of terms of the flat scores matrix and
#' the corresponding subgraph is computed as well.
#' @param norm boolean value: 
#' \itemize{
#' \item \code{TRUE} (def.): the flat scores matrix has been already normalized in according to a normalization method;	
#' \item \code{FALSE}: the flat scores matrix has not been normalized yet. See the parameter \code{norm} for which normalization can be applied.
#' }
#' @param norm.type can be one of the following three values:
#'  \enumerate{
#'  \item \code{NULL} (def.): set \code{norm.type} to \code{NULL} if and only if the parameter \code{norm} is set to \code{TRUE};
#'  \item \code{MaxNorm}: each score is divided for the maximum of each class;
#'  \item \code{Qnorm}: quantile normalization. \pkg{preprocessCore} package is used. 
#'  }
#' @param folds number of folds of the cross validation on which computing the performance metrics averaged across folds (\code{def. 5}).
#' If \code{folds=NULL}, the performance metrics are computed one-shot, otherwise the performance metrics are averaged across folds.
#' @param seed initialization seed for the random generator to create folds (\code{def. 23}). If \code{NULL} folds are generated without seed 
#' initialization. 
#' @param n.round number of rounding digits to be applied to the hierarchical scores matrix (\code{def. 3}). It is used for choosing 
#' the best threshold on the basis of the best F-measure
#' @param f.criterion character. Type of F-measure to be used to select the best F-measure. Two possibilities:
#' \enumerate{
#' \item \code{F} (def.): corresponds to the harmonic mean between the average precision and recall
#' \item \code{avF}: corresponds to the per-example \code{F-score} averaged across all the examples
#' }
#' @param rec.levels a vector with the desired recall levels (\code{def:} \code{from:0.1}, \code{to:0.9}, \code{by:0.1}) to compute the 
#' the Precision at fixed Recall level (PXR)
#' @param flat.file name of the file containing the flat scores matrix to be normalized or already normalized (without rda extension)
#' @param ann.file name of the file containing the the label matrix of the examples (without rda extension)
#' @param dag.file name of the file containing the graph that represents the hierarchy of the classes (without rda extension)
#' @param ind.test.set name of the file containing a vector of integer numbers corresponding to the indices of the elements (rows) of scores 
#' matrix to be used in the	test set 
#' @param ind.dir relative path to folder where \code{ind.test.set} is stored
#' @param flat.dir relative path where flat scores matrix is stored
#' @param ann.dir relative path where annotation matrix is stored
#' @param dag.dir relative path where graph is stored
#' @param hierScore.dir relative path where the hierarchical scores matrix must be stored
#' @param perf.dir relative path where all the performance measures must be stored
#' @return Two \code{rda} files stored in the respective output directories:
#' \enumerate{
#' 	\item \code{Hierarchical Scores Results}: a matrix with examples on rows and classes on columns representing the computed hierarchical scores 
#' 	for each example and for each considered class. It is stored in the \code{hierScore.dir} directory.
#' 	\item \code{Performance Measures}: \emph{flat} and \emph{hierarchical} performace results:
#' 	\enumerate{
#' 		\item AUPRC results computed though \code{AUPRC.single.over.classes} (\code{\link{AUPRC}});
#'		\item AUROC results computed through \code{AUROC.single.over.classes} (\code{\link{AUROC}}); 
#' 		\item PXR results computed though \code{PXR.at.multiple.recall.levels.over.classes} (\code{\link{PXR}});
#' 		\item FMM results computed though \code{compute.Fmeasure.multilabel} (\code{\link{FMM}}); 
#' }}
#' It is stored in the \code{perf.dir} directory.
#' @seealso \code{\link{HTD-DAG}}
#' @export
#' @examples
#' data(graph);
#' data(scores);
#' data(labels);
#' data(test.index);
#' if (!dir.exists("data")){
#' 	dir.create("data");
#' }
#' if (!dir.exists("results")){
#' 	dir.create("results");
#' }
#' save(g,file="data/graph.rda");
#' save(L,file="data/labels.rda");
#' save(S,file="data/scores.rda");
#' save(test.index, file="data/test.index.rda");
#' ind.dir <- dag.dir <- flat.dir <- ann.dir <- "data/";
#' hierScore.dir <- perf.dir <- "results/";
#' rec.levels <- seq(from=0.1, to=1, by=0.1);
#' ind.test.set <- "test.index";
#' dag.file <- "graph";
#' flat.file <- "scores";
#' ann.file <- "labels";
#' Do.HTD.holdout(norm=FALSE, norm.type="MaxNorm", n.round=3, f.criterion ="F", folds=NULL, seed=23,
#' rec.levels=rec.levels, flat.file=flat.file, ann.file=ann.file, dag.file=dag.file, 
#' ind.test.set=ind.test.set, ind.dir=ind.dir, flat.dir=flat.dir, ann.dir=ann.dir, dag.dir=dag.dir, 
#' hierScore.dir=hierScore.dir, perf.dir=perf.dir);
Do.HTD.holdout <- function(norm=TRUE, norm.type=NULL, folds=5, seed=23, n.round=3, f.criterion ="F", 
	rec.levels=seq(from=0.1, to=1, by=0.1), flat.file=flat.file, ann.file=ann.file, dag.file=dag.file, 
	ind.test.set=ind.test.set, ind.dir=ind.dir, flat.dir=flat.dir, ann.dir=ann.dir, dag.dir=dag.dir, 
	hierScore.dir=hierScore.dir, perf.dir=perf.dir){
	
	## Setting Check
	if(norm==FALSE && is.null(norm.type))
		stop("HTD: If norm is set to FALSE, you need also to specify a normalization method among those available", call.=FALSE)
	if(norm==TRUE && !is.null(norm.type))
		warning("HTD: If norm is set to TRUE, the input flat matrix is already normalized.", 
			paste0(" Set norm.type to NULL and not to '", norm.type, "' to avoid this warning message"), call.=FALSE);
	
	## Loading Data
	## loading examples indices of the test set
	ind.set <- paste0(ind.dir, ind.test.set, ".rda");
	ind.test <- get(load(ind.set));

	## loading dag
	dag.path <- paste0(dag.dir, dag.file,".rda");
	g <- get(load(dag.path));
	## compute root node 
	root <- root.node(g);

	## loading annotation matrix
	ann.path <- paste0(ann.dir, ann.file,".rda");
	ann <- get(load(ann.path));
	## removing root node from annotation table if it exists
	if(root %in% colnames(ann))
		ann <- ann[,-which(colnames(ann)==root)];

	## loading flat matrix
	flat.path <- paste0(flat.dir, flat.file,".rda");
	if(norm){
		S <- get(load(flat.path));
		## removing root node from flat norm matrix if it exists
		if(root %in% colnames(S))
			S <- S[,-which(colnames(S)==root)];
	}else{
		S <- get(load(flat.path));
		S <- scores.normalization(norm.type=norm.type, S);
		cat(norm.type, "NORMALIZATION: DONE", "\n");
		## removing root node from flat norm matrix if it exists
		if(root %in% colnames(S))
			S <- S[,-which(colnames(S)==root)];
	}

	## check if |flat matrix classes| = |annotation matrix classes| 
	## if not the classes of annotation matrix are shrinked to those of flat matrix
	class.check <- ncol(S)!=ncol(ann);
	if(class.check){
		ann <- ann[,colnames(S)];
		nd <- colnames(S);
		g <- do.subgraph(nd, g, edgemode="directed");
	}

	## scores flat matrix shrinked to test test
	S <- S[ind.test,];
	
	## annotation table  shrinked to test test 
	ann <- ann[ind.test,];

	## Compute FLAT PRC, AUC, PXR (average and per class) and FMM (average and per-example) one-shoot or cross-validated 
	PRC.flat <- AUPRC.single.over.classes(ann, S, folds=folds, seed=seed);
	AUC.flat <- AUROC.single.over.classes(ann, S, folds=folds, seed=seed);
	PXR.flat <- PXR.at.multiple.recall.levels.over.classes(ann, S, rec.levels=rec.levels, folds=folds, seed=seed);
	FMM.flat <- compute.Fmeasure.multilabel(ann, S, n.round=n.round, f.criterion=f.criterion, verbose=FALSE,
		b.per.example=TRUE, folds=folds, seed=seed);

	## Hierarchical Top Down Correction 
	S <- htd(S,g,root);
	
	## Compute HIER PRC, AUC, PXR (average and per class) and FMM (average and per-example) one-shoot or cross-validated 
	PRC.hier <- AUPRC.single.over.classes(ann, S, folds=folds, seed=seed);
	AUC.hier <- AUROC.single.over.classes(ann, S, folds=folds, seed=seed);
	PXR.hier <- PXR.at.multiple.recall.levels.over.classes(ann, S, rec.levels=rec.levels, folds=folds, seed=seed);
	FMM.hier <- compute.Fmeasure.multilabel(ann, S, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, 
		b.per.example=TRUE, folds=folds, seed=seed);
	S.hier <- S;
	rm(S);
	
	## Storing Results 
	if(norm){
		save(S.hier, file=paste0(hierScore.dir, flat.file, ".hierScores.htd.holdout.rda"), compress=TRUE);
		save(PRC.flat, PRC.hier, AUC.flat, AUC.hier, PXR.flat, PXR.hier, FMM.flat, FMM.hier,
			file=paste0(perf.dir, "PerfMeas.", flat.file, ".hierScores.htd.holdout.rda"), compress=TRUE);
	}else{
		save(S.hier, file=paste0(hierScore.dir, norm.type,".", flat.file, ".hierScores.htd.holdout.rda"), compress=TRUE);	
		save(PRC.flat, PRC.hier, AUC.flat, AUC.hier, PXR.flat, PXR.hier, FMM.flat, FMM.hier,
		 file=paste0(perf.dir, "PerfMeas.", norm.type,".", flat.file, ".hierScores.htd.holdout.rda"), compress=TRUE);
	}
}
