##*********##
## HTD-DAG ##
##*********##

#' @name HTD-DAG
#' @title HTD-DAG
#' @description Implementetion of a top-down procedure to correct the scores of the hierarchy according to the 
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
htd <- function(S,g, root="00"){
	levels <- graph.levels(g,root);
	# a dummy root is added
	if(!(root %in% colnames(S))){
		max.score <- max(S);
		z <- rep(max.score,nrow(S));
		S <- cbind(z,S);
		colnames(S)[1] <- root;
	}
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
#' @param norm boolean value: 
#' \itemize{
#' \item \code{TRUE} (def.): the flat scores matrix has been already normalized in according to a normalization method;	
#' \item \code{FALSE}: the flat scores matrix has not been normalized yet. See the parameter \code{norm.type} for which normalization can be applied.
#' }
#' @param norm.type can be one of the following three values:
#'  \enumerate{
#'  \item \code{NONE} (def.): set \code{norm.type} to \code{NONE} if and only if the parameter \code{norm} is set to \code{TRUE};
#'  \item \code{MaxNorm}: each score is divided w.r.t. the max of each class;
#'  \item \code{Qnorm}: quantile normalization. \pkg{preprocessCore} package is used. 
#'  }
#' @param flat.file name of the file containing the flat scores matrix to be normalized or already normalized (without rda extension)
#' @param ann.file name of the file containing the the label matrix of the examples (without rda extension)
#' @param dag.file name of the file containing the graph that represents the hierarchy of the classes (without rda extension)
#' @param flat.dir relative path where flat scores matrix is stored
#' @param ann.dir relative path where annotation matrix is stored
#' @param dag.dir relative path where graph is stored
#' @param flat.norm.dir relative path where flat normalized scores matrix must be stored. Use this parameter if and only if \code{norm} is
#' set to \code{FALSE}, otherwise set \code{flat.norm.dir} to \code{NULL} (def.)
#' @param n.round number of rounding digits to be applied to the hierarchical scores matrix (\code{def. 3}). It is used for choosing 
#' the best threshold on the basis of the best F-measure
#' @param f.criterion character. Type of F-measure to be used to select the best F-measure. Two possibilities:
#' \enumerate{
#' \item \code{F} (def.): corresponds to the harmonic mean between the average precision and recall
#' \item \code{avF}: corresponds to the per-example \code{F-score} averaged across all the examples
#' }
#' @param hierScore.dir relative path where the hierarchical scores matrix must be stored
#' @param perf.dir relative path where the term-centric and protein-centric measures must be stored
#' @return Five \code{rda} files stored in the respective output directories:
#' \enumerate{
#' \item \code{hierarchical scores matrix}: a matrix with examples on rows and classes on columns representing the computed hierarchical scores 
#' for each example and for each considered class. This file is stored in \code{hierScore.dir} directory.
#' \item \code{FMM} (F-Measure Multilabel) \code{results}: \code{F-score} computed by \code{find.best.f} function. 
#' Both \emph{flat} and \emph{hierarchical} results are reported. This file is stored in \code{perf.dir} directory.
#' \item \code{PRC} (area under Precision-Recall Curve) \code{results}: \code{PRC} computed by \pkg{precrec} package. 
#' Both \emph{flat} and \emph{hierarchical} results are reported. This file is stored in \code{perf.dir} directory.
#' \item \code{AUC} (Area Under ROC Curve) \code{results}: \code{AUC} computed by \pkg{precrec} package. 
#' Both \emph{flat} and \emph{hierarchical} results are reported. This file is stored in \code{perf.dir} directory.
#' \item \code{PXR} (Precision at fixed Recall levels) average and per classes: \code{PXR} computed by \pkg{PerfMeas} package. 
#' It is stored in \code{perf.dir} directory.
#' }
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
#' dag.dir <- flat.dir <- flat.norm.dir <- ann.dir <- "data/";
#' hierScore.dir <- perf.dir <- "results/";
#' dag.file <- "graph";
#' flat.file <- "scores";
#' ann.file <- "labels";
#' Do.HTD(norm=FALSE, norm.type= "MaxNorm", flat.file=flat.file, ann.file=ann.file, 
#' dag.file=dag.file, flat.dir=flat.dir, ann.dir=ann.dir, dag.dir=dag.dir, 
#' flat.norm.dir=flat.norm.dir, n.round=3, f.criterion ="F", hierScore.dir=hierScore.dir, 
#' perf.dir=perf.dir);
Do.HTD <- function(norm=TRUE, norm.type= "NONE", flat.file=flat.file, ann.file=ann.file, dag.file=dag.file, flat.dir=flat.dir, 
	ann.dir=ann.dir, dag.dir=dag.dir, flat.norm.dir=NULL, n.round=3, f.criterion ="F", hierScore.dir=hierScore.dir, perf.dir=perf.dir){
	
	## Loading Data ############
	## loading dag
	dag.path <- paste0(dag.dir, dag.file,".rda");
	g <- get(load(dag.path));
	
	##root node
	root <- root.node(g);

	## loading flat scores matrix relative to a specific subontology
	flat.path <- paste0(flat.dir, flat.file,".rda");
	if(norm){
		S <- get(load(flat.path));
		
		## removing root node from flat norm matrix if it exists
		if(root %in% colnames(S)){
			S <- S[,-which(colnames(S)==root)];
		}
	}else{
		Do.FLAT.scores.normalization(norm.type=norm.type, flat.file=flat.file, dag.file=dag.file, flat.dir=flat.dir, 
			dag.dir=dag.dir, flat.norm.dir=flat.norm.dir);
		flat.path <- paste0(flat.norm.dir, norm.type,".",flat.file,".rda");
		S <- get(load(flat.path));
	}

	## loading annotation matrix
	ann.path <- paste0(ann.dir, ann.file,".rda");
	ann <- get(load(ann.path));
	gc();

	## removing root node from annotation table 
	ann <- ann[,-which(colnames(ann)==root)];

	## Computing FLAT Performances
	## FLAT AUC computed by precrec package
	AUC.flat <- AUROC.single.over.classes(ann, S); gc();
	
	## FLAT PxRs computed by PerfMeas pacakge
	PXR.flat <- precision.at.multiple.recall.level.over.classes(ann, S); 

	## F.measure: Computing Flat Examples-Measures 
	FMM.flat <- find.best.f(ann, S, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=TRUE); 

	## FLAT PRC computed by precrec package (more precise and accurate than PerfMeas)
	PRC.flat <- AUPRC.single.over.classes(ann, S); gc();

	## Hierarchical Top Down Correction ####################
	## in this way we fill memory because we store two double-float matrix. Solution overwrite!! we have already calculated the flat performances..
	# S.htd <- htd(S,g,root);
	S <- htd(S, g, root);
	
	## Computing Hier Performances
	## Hierarchical AUC (average and per.class) computed by precrec package
	AUC.hier <- AUROC.single.over.classes(ann, S); gc();

	## Hierarchical PxR at fixed recall levels 
	PXR.hier <- precision.at.multiple.recall.level.over.classes(ann, S); gc();

	## Computing Hierarchical Examples-Measures 
	FMM.hier <- find.best.f(ann, S, n.round=n.round, f.criterion =f.criterion, verbose=FALSE, b.per.example=TRUE);

	## Hierarchical PRC (average and per.class) computed by precrec package
	PRC.hier <- AUPRC.single.over.classes(ann, S); 

	## storing the hierarchical matrix
	S.hier <- S;
	rm(S); gc();

	## Storing Results #########
	if(norm){
		save(S.hier, file=paste0(hierScore.dir, flat.file, ".hierScores.htd.rda"), compress=TRUE);
		save(AUC.flat, AUC.hier, file=paste0(perf.dir, "AUC.", flat.file, ".hierScores.htd.rda"), compress=TRUE);
		save(PXR.flat, PXR.hier, file=paste0(perf.dir, "PXR.", flat.file, ".hierScores.htd.rda"), compress=TRUE);
		save(FMM.flat, FMM.hier, file=paste0(perf.dir, "FMM.", flat.file, ".hierScores.htd.rda"), compress=TRUE);
		save(PRC.flat, PRC.hier, file=paste0(perf.dir, "PRC.", flat.file, ".hierScores.htd.rda"), compress=TRUE);
	}else{
		save(S.hier, file=paste0(hierScore.dir, norm.type,".", flat.file, ".hierScores.htd.rda"), compress=TRUE);	
		save(AUC.flat, AUC.hier, file=paste0(perf.dir, "AUC.", norm.type,".", flat.file, ".hierScores.htd.rda"), compress=TRUE);	
		save(PXR.flat, PXR.hier, file=paste0(perf.dir, "PXR.", norm.type,".", flat.file, ".hierScores.htd.rda"), compress=TRUE);	
		save(FMM.flat, FMM.hier, file=paste0(perf.dir, "FMM.", norm.type,".", flat.file, ".hierScores.htd.rda"), compress=TRUE);
		save(PRC.flat, PRC.hier, file=paste0(perf.dir, "PRC.", norm.type,".", flat.file, ".hierScores.htd.rda"), compress=TRUE);
	}
}

#' @title HTD-DAG holdout
#' @description High level function to correct the computed scores in a hierarchy according to the HTD-DAG algorithm applying a 
#' classical holdout procedure
#' @param norm boolean value: 
#' \itemize{
#' \item \code{TRUE} (def.): the flat scores matrix has been already normalized in according to a normalization method;	
#' \item \code{FALSE}: the flat scores matrix has not been normalized yet. See the parameter \code{norm} for which normalization can be applied.
#' }
#' @param norm.type can be one of the following three values:
#'  \enumerate{
#'  \item \code{NONE} (def.): set \code{norm.type} to \code{NONE} if and only if the parameter \code{norm} is set to \code{TRUE};
#'  \item \code{MaxNorm}: each score is divided w.r.t. the max of each class;
#'  \item \code{Qnorm}: quantile normalization. \pkg{preprocessCore} package is used. 
#'  }
#' @param flat.file name of the file containing the flat scores matrix to be normalized or already normalized (without rda extension)
#' @param ann.file name of the file containing the the label matrix of the examples (without rda extension)
#' @param dag.file name of the file containing the graph that represents the hierarchy of the classes (without rda extension)
#' @param ind.test.set name of the file containing a vector of integer numbers corresponding to the indices of the elements (rows) of scores 
#' matrix to be used in the	test set 
#' @param ind.dir relative path to folder where \code{ind.test.set} is stored
#' @param flat.dir relative path where flat scores matrix is stored
#' @param ann.dir relative path where annotation matrix is stored
#' @param dag.dir relative path where graph is stored
#' @param flat.norm.dir relative path where flat normalized scores matrix must be stored. Use this parameter if and only if \code{norm} is
#' set to \code{FALSE}, otherwise set \code{flat.norm.dir} to \code{NULL} (def.)
#' @param n.round number of rounding digits to be applied to the hierarchical scores matrix (\code{def. 3}). It is used for choosing 
#' the best threshold on the basis of the best F-measure
#' @param f.criterion character. Type of F-measure to be used to select the best F-measure. Two possibilities:
#' \enumerate{
#' \item \code{F} (def.): corresponds to the harmonic mean between the average precision and recall
#' \item \code{avF}: corresponds to the per-example \code{F-score} averaged across all the examples
#' }
#' @param hierScore.dir relative path where the hierarchical scores matrix must be stored
#' @param perf.dir relative path where the term-centric and protein-centric measures must be stored
#' @return Five \code{rda} files stored in the respective output directories:
#' \enumerate{
#' \item \code{hierarchical scores matrix}: a matrix with examples on rows and classes on columns representing the computed hierarchical scores 
#' for each example and for each considered class. This file is stored in \code{hierScore.dir} directory.
#' \item \code{FMM} (F-Measure Multilabel) \code{results}: \code{F-score} computed by \code{find.best.f} function. 
#' Both \emph{flat} and \emph{hierarchical} results are reported. This file is stored in \code{perf.dir} directory.
#' \item \code{PRC} (area under Precision-Recall Curve) \code{results}: \code{PRC} computed by \pkg{precrec} package. 
#' Both \emph{flat} and \emph{hierarchical} results are reported. This file is stored in \code{perf.dir} directory.
#' \item \code{AUC} (Area Under ROC Curve) \code{results}: \code{AUC} computed by \pkg{precrec} package. 
#' Both \emph{flat} and \emph{hierarchical} results are reported. This file is stored in \code{perf.dir} directory.
#' \item \code{PXR} (Precision at fixed Recall levels) average and per classes: \code{PXR} computed by \pkg{PerfMeas} package. 
#' It is stored in \code{perf.dir} directory.
#' }
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
#' ind.dir <- dag.dir <- flat.dir <- flat.norm.dir <- ann.dir <- "data/";
#' hierScore.dir <- perf.dir <- "results/";
#' ind.test.set <- "test.index";
#' dag.file <- "graph";
#' flat.file <- "scores";
#' ann.file <- "labels";
#' Do.HTD.holdout(norm=FALSE, norm.type= "MaxNorm", flat.file=flat.file, ann.file=ann.file, 
#' dag.file=dag.file, ind.test.set=ind.test.set, ind.dir=ind.dir, flat.dir=flat.dir, 
#' ann.dir=ann.dir, dag.dir=dag.dir, flat.norm.dir=flat.norm.dir, n.round=3, f.criterion ="F", 
#' hierScore.dir=hierScore.dir, perf.dir=perf.dir);
Do.HTD.holdout <- function(norm=TRUE, norm.type= "NONE", flat.file=flat.file, ann.file=ann.file, dag.file=dag.file, 
	ind.test.set=ind.test.set, ind.dir=ind.dir, flat.dir=flat.dir, ann.dir=ann.dir, dag.dir=dag.dir, flat.norm.dir=NULL, 
	n.round=3, f.criterion ="F", hierScore.dir=hierScore.dir, perf.dir=perf.dir){

	## Loading Data ############
	# loading examples indices of the test set
	ind.set <- paste0(ind.dir, ind.test.set, ".rda");
	ind.test <- get(load(ind.set));

	## loading dag
	dag.path <- paste0(dag.dir, dag.file,".rda");
	g <- get(load(dag.path));
	
	##root node
	root <- root.node(g);

	## loading flat scores matrix relative to a specific subontology
	flat.path <- paste0(flat.dir, flat.file,".rda");
	if(norm){
		S <- get(load(flat.path));
		gc();	##in order to save ram memory..

		## removing root node from flat norm matrix if it exists
		if(root %in% colnames(S)){
			S <- S[,-which(colnames(S)==root)];
		}
	}else{
		Do.FLAT.scores.normalization(norm.type=norm.type, flat.file=flat.file, dag.file=dag.file, flat.dir=flat.dir, 
			dag.dir=dag.dir, flat.norm.dir=flat.norm.dir);
		flat.path <- paste0(flat.norm.dir, norm.type,".",flat.file,".rda");
		S <- get(load(flat.path));
	}

	## shrinking the size of S to the examples of test set
	S <- S[ind.test,];

	## loading annotation matrix
	ann.path <- paste0(ann.dir, ann.file,".rda");
	ann <- get(load(ann.path));

	## removing root node from annotation table and shrinking the size of annotation table to the examples of test set
	ann <- ann[ind.test,-which(colnames(ann)==root)];

	## Computing FLAT Performances
	## FLAT AUC computed by precrec package
	AUC.flat <- AUROC.single.over.classes(ann, S); 	
			
	## FLAT PxRs computed by PerfMeas pacakge
	PXR.flat <- precision.at.multiple.recall.level.over.classes(ann, S);

	## F.measure: Computing Flat Examples-Measures 
	FMM.flat <- find.best.f(ann, S, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=TRUE); 

	## FLAT PRC computed by precrec package 
	PRC.flat <- AUPRC.single.over.classes(ann, S);  

	## Hierarchical Top Down Correction ####################
	## in this way we fill memory because we store two double-float matrix. Solution overwrite!! we have already calculated the flat performances..
	# S.htd <- htd(S,g,root);
	S <- htd(S,g,root);
	
	## Computing Hier Performances
	## Hierarchical AUC (average and per.class) computed by precrec package
	AUC.hier <- AUROC.single.over.classes(ann, S); 
		
	## Hierarchical PxR at fixed recall levels 
	PXR.hier <- precision.at.multiple.recall.level.over.classes(ann, S); 

	## Computing Hierarchical Examples-Measures 
	FMM.hier <- find.best.f(ann, S, n.round=n.round, f.criterion =f.criterion, verbose=FALSE, b.per.example=TRUE);	

	## Hierarchical PRC (average and per.class) computed by precrec package
	PRC.hier <- AUPRC.single.over.classes(ann, S);  

	## storing the hierarchical matrix
	S.hier <- S;
	rm(S); gc();
	
	## Storing Results #########
	if(norm){
		save(S.hier, file=paste0(hierScore.dir, flat.file, ".hierScores.htd.holdout.rda"), compress=TRUE);
		save(AUC.flat, AUC.hier, file=paste0(perf.dir, "AUC.", flat.file, ".hierScores.htd.holdout.rda"), compress=TRUE);
		save(PXR.flat, PXR.hier, file=paste0(perf.dir, "PXR.", flat.file, ".hierScores.htd.holdout.rda"), compress=TRUE);
		save(FMM.flat, FMM.hier, file=paste0(perf.dir, "FMM.", flat.file, ".hierScores.htd.holdout.rda"), compress=TRUE);
		save(PRC.flat, PRC.hier, file=paste0(perf.dir, "PRC.", flat.file, ".hierScores.htd.holdout.rda"), compress=TRUE);
	}else{
		save(S.hier, file=paste0(hierScore.dir, norm.type,".", flat.file, ".hierScores.htd.holdout.rda"), compress=TRUE);	
		save(AUC.flat, AUC.hier, file=paste0(perf.dir, "AUC.", norm.type,".", flat.file, ".hierScores.htd.holdout.rda"), compress=TRUE);	
		save(PXR.flat, PXR.hier, file=paste0(perf.dir, "PXR.", norm.type,".", flat.file, ".hierScores.htd.holdout.rda"), compress=TRUE);	
		save(FMM.flat, FMM.hier, file=paste0(perf.dir, "FMM.", norm.type,".", flat.file, ".hierScores.htd.holdout.rda"), compress=TRUE);
		save(PRC.flat, PRC.hier, file=paste0(perf.dir, "PRC.", norm.type,".", flat.file, ".hierScores.htd.holdout.rda"), compress=TRUE);
	}
}
