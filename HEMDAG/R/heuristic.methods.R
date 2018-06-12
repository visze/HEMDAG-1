###################################
##  Obozinski Heuristic Methods  ##
###################################

#' @name Heuristic-Methods
#' @aliases heuristicMAX
#' @aliases heuristicAND
#' @aliases heuristicOR
#' @title Obozinski Heuristic Methods 
#' @description Implementation of the Heuristic Methods MAX, AND, OR (\cite{Obozinski et al., Genome Biology, 2008, 
#' \href{https://genomebiology.biomedcentral.com/articles/10.1186/gb-2008-9-s1-s6}{doi:10.1186/gb-2008-9-s1-s6}})
#' @details Heuristic Methods:
#' \enumerate{
#'	\item \bold{MAX}: reports the largest logist regression (LR) value of self and all descendants: \eqn{p_i = max_{j \in descendants(i)} \hat{p_j}};
#' 	\item \bold{AND}: reports the product of LR values of all ancestors and self. This is equivalent to computing the probability that all 
#' ancestral terms are "on" assuming that, conditional on the data, all predictions are independent: \eqn{p_i = \prod_{j \in ancestors(i)} \hat{p_j}};
#'	\item \bold{OR}: computes the probability that at least one of the descendant terms is "on" assuming again that, conditional on the data, 
#' all predictions are independent: \eqn{1 - p_i = \prod_{j \in descendants(i)} (1 - \hat{p_j})};
#' }
#' @param S a named flat scores matrix with examples on rows and classes on columns
#' @param g a graph of class \code{graphNEL}. It represents the hierarchy of the classes
#' @param root name of the class that it is the top-level (root) of the hierarchy (\code{def:00})
#' @return a matrix with the scores of the classes corrected according to the chosen heuristic algorithm
#' @export
#' @examples
#' data(graph);
#' data(scores);
#' data(labels);
#' root <- root.node(g);
#' S.heuristicMAX <- heuristicMAX(S,g,root);
#' S.heuristicAND <- heuristicAND(S,g,root);
#' S.heuristicOR <- heuristicOR(S,g,root);
heuristicMAX <- function(S, g, root="00"){
	if(!(root %in% colnames(S))) {
	  max.score <- max(S);
	  z <- rep(max.score,nrow(S));
	  S <- cbind(z,S);
	  colnames(S)[1] <- root;
	}
	desc <- build.descendants(g);
	for(i in 1:length(desc)){
	  curr.nd <- S[,names(desc[i])];
	  progeny <- as.matrix(S[,desc[[i]]]);
	  curr.nd <- apply(progeny, 1, max);
	  S[,names(desc[i])] <- curr.nd;
	}
	S <- S[,-which(colnames(S)==root)]; 
	return(S);
}

#' @rdname Heuristic-Methods
#' @export 
heuristicAND <- function(S, g, root="00"){
	if(!(root %in% colnames(S))) {
	  max.score <- max(S);
	  z <- rep(max.score,nrow(S));
	  S <- cbind(z,S);
	  colnames(S)[1] <- root;
	}
	S.corr <- S;
	
	anc <- build.ancestors(g);
	for(i in 1:length(anc)){
	  if (length(anc[[i]]) > 1) {  # ancestors of i include also i
	    curr.nd <- S[,names(anc[i])];
	    forefathers <- as.matrix(S[,anc[[i]]]);
	    curr.nd <- apply(forefathers, 1, prod);
	    S.corr[,names(anc[i])] <- curr.nd;  
	  }	
	}
	S.corr <- S.corr[,-which(colnames(S.corr)==root)];
	return(S.corr);
}

#' @rdname Heuristic-Methods
#' @export 
heuristicOR <- function(S, g, root="00"){
	if(!(root %in% colnames(S))) {
	  max.score <- max(S);
	  z <- rep(max.score,nrow(S));
	  S <- cbind(z,S);
	  colnames(S)[1] <- root;
	}
	S.corr <- S;
	desc <- build.descendants(g);
	for(i in 1:length(desc)){
	  if(length(desc[[i]]) > 1){  # descendants of i include also i
		comp.progeny <- 1 - as.matrix(S[,desc[[i]]]);
		curr.nd <- apply(comp.progeny, 1, prod);		
		S.corr[,names(desc[i])] <- 1 - curr.nd;
	  }
	}
	S.corr <- S.corr[,-which(colnames(S.corr)==root)]; 
	return(S.corr);
}

##**********************##
## DO HEURISTIC METHODS ##
##**********************##

#' @title Do Heuristic Methods
#' @seealso \code{\link{Heuristic-Methods}}
#' @description High level function to compute the hierarchical heuristic methods MAX, AND, OR (Heuristic Methods MAX, AND, OR (\cite{Obozinski et al., 
#' Genome Biology, 2008})
#' @param heuristic.fun can be one of the following three values:
#' \enumerate{
#' 	\item heuristicMAX: run the MAX heuristic method;
#' 	\item heuristicAND: run theAND heuristic method;
#' 	\item heuristicOR: run theOR heuristic method;
#' }
#' @param norm boolean value: 
#' \itemize{
#' \item \code{TRUE} (def.): the flat scores matrix has been already normalized in according to a normalization method;	
#' \item \code{FALSE}: the flat scores matrix has not been normalized yet. See the parameter \code{norm.type} for which normalization can be applied.
#' }
#' @param norm.type can be one of the following three values: 
#' \enumerate{
#' 	\item heuristicMAX: run the heuristic method MAX;
#' 	\item heuristicAND: run the heuristic method AND;
#' 	\item heuristicOR: run the heuristic method OR;
#' }
#' @param flat.file name of the file containing the flat scores matrix to be normalized or already normalized (without rda extension)
#' @param ann.file name of the file containing the the label matrix of the examples (without rda extension)
#' @param dag.file name of the file containing the graph that represents the hierarchy of the classes (without rda extension)
#' @param flat.dir relative path where flat scores matrix is stored
#' @param ann.dir relative path where annotation matrix is stored
#' @param dag.dir relative path where graph is stored
#' @param flat.norm.dir relative path where flat normalized scores matrix must be strored. Use this parameter if and only if \code{norm} is
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
#' @return Five \code{rda} files stored in the rispective output directories:
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
#' Do.heuristic.methods(heuristic.fun=heuristicAND, norm=FALSE, 
#' norm.type="MaxNorm", flat.file=flat.file, ann.file=ann.file, dag.file=dag.file, 
#' flat.dir=flat.dir, ann.dir=ann.dir, dag.dir=dag.dir, flat.norm.dir=flat.norm.dir, 
#' n.round=3, f.criterion ="F", hierScore.dir=hierScore.dir, perf.dir=perf.dir);
Do.heuristic.methods <- function(heuristic.fun=heuristic.fun, norm=TRUE, norm.type= "NONE", flat.file=flat.file, 
	ann.file=ann.file, dag.file=dag.file, flat.dir=flat.dir, ann.dir=ann.dir, dag.dir=dag.dir, flat.norm.dir=NULL, n.round=3, 
	f.criterion ="F", hierScore.dir=hierScore.dir, perf.dir=perf.dir){

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

	## Obozinski's Hierarchical Heuristic Methods ####################
	S <- heuristic.fun(S, g, root);
	
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
	heuristic.name <- as.character(substitute(heuristic.fun));
	if(norm){
		save(S.hier, file=paste0(hierScore.dir, flat.file, ".hierScores",heuristic.name,".rda"), compress=TRUE);
		save(AUC.flat, AUC.hier, file=paste0(perf.dir, "AUC.", flat.file, ".hierScores.",heuristic.name,".rda"), compress=TRUE);
		save(PXR.flat, PXR.hier, file=paste0(perf.dir, "PXR.", flat.file, ".hierScores.",heuristic.name,".rda"), compress=TRUE);
		save(FMM.flat, FMM.hier, file=paste0(perf.dir, "FMM.", flat.file, ".hierScores.",heuristic.name,".rda"), compress=TRUE);
		save(PRC.flat, PRC.hier, file=paste0(perf.dir, "PRC.", flat.file, ".hierScores.",heuristic.name,".rda"), compress=TRUE);
	}else{
		save(S.hier, file=paste0(hierScore.dir, norm.type,".", flat.file, ".hierScores.",heuristic.name,".rda"), compress=TRUE);	
		save(AUC.flat, AUC.hier, file=paste0(perf.dir, "AUC.", norm.type,".", flat.file, ".hierScores.",heuristic.name,".rda"), compress=TRUE);	
		save(PXR.flat, PXR.hier, file=paste0(perf.dir, "PXR.", norm.type,".", flat.file, ".hierScores.",heuristic.name,".rda"), compress=TRUE);	
		save(FMM.flat, FMM.hier, file=paste0(perf.dir, "FMM.", norm.type,".", flat.file, ".hierScores.",heuristic.name,".rda"), compress=TRUE);
		save(PRC.flat, PRC.hier, file=paste0(perf.dir, "PRC.", norm.type,".", flat.file, ".hierScores.",heuristic.name,".rda"), compress=TRUE);
	}
}

#' @title Do Heuristic Methods holdout
#' @description High level function to compute the hierarchical heuristic methods MAX, AND, OR (Heuristic Methods MAX, AND, OR (\cite{Obozinski et al., 
#' Genome Biology, 2008}) applying a classical holdout procedure
#' @param heuristic.fun can be one of the following three values:
#' \enumerate{
#' 	\item heuristicMAX: run the heuristic method MAX;
#' 	\item heuristicAND: run the heuristic method AND;
#' 	\item heuristicOR: run the heuristic method OR;
#' }
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
#' @param flat.norm.dir relative path where flat normalized scores matrix must be strored. Use this parameter if and only if \code{norm} is
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
#' @return Five \code{rda} files stored in the rispective output directories:
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
#' Do.heuristic.methods.holdout(heuristic.fun=heuristicAND, norm=FALSE, 
#' norm.type="MaxNorm", flat.file=flat.file, ann.file=ann.file, dag.file=dag.file, 
#' ind.test.set=ind.test.set, ind.dir=ind.dir, flat.dir=flat.dir, ann.dir=ann.dir, 
#' dag.dir=dag.dir, flat.norm.dir=flat.norm.dir, n.round=3, f.criterion="F", 
#' hierScore.dir=hierScore.dir, perf.dir=perf.dir);
Do.heuristic.methods.holdout <- function(heuristic.fun=heuristic.fun, norm=TRUE, norm.type= "NONE", flat.file=flat.file, 
	ann.file=ann.file, dag.file=dag.file, ind.test.set=ind.test.set, ind.dir=ind.dir, flat.dir=flat.dir, ann.dir=ann.dir, 
	dag.dir=dag.dir, flat.norm.dir=NULL, n.round=3, f.criterion ="F", hierScore.dir=hierScore.dir, perf.dir=perf.dir){

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

	## Obozinski's Hierarchical Heuristic Methods ####################
	S <- heuristic.fun(S, g, root);
	
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
	heuristic.name <- as.character(substitute(heuristic.fun));
	if(norm){
		save(S.hier, file=paste0(hierScore.dir, flat.file, ".hierScores",heuristic.name,".holdout.rda"), compress=TRUE);
		save(AUC.flat, AUC.hier, file=paste0(perf.dir, "AUC.", flat.file, ".hierScores.",heuristic.name,".holdout.rda"), compress=TRUE);
		save(PXR.flat, PXR.hier, file=paste0(perf.dir, "PXR.", flat.file, ".hierScores.",heuristic.name,".holdout.rda"), compress=TRUE);
		save(FMM.flat, FMM.hier, file=paste0(perf.dir, "FMM.", flat.file, ".hierScores.",heuristic.name,".holdout.rda"), compress=TRUE);
		save(PRC.flat, PRC.hier, file=paste0(perf.dir, "PRC.", flat.file, ".hierScores.",heuristic.name,".holdout.rda"), compress=TRUE);
	}else{
		save(S.hier, file=paste0(hierScore.dir, norm.type,".", flat.file, ".hierScores.",heuristic.name,".holdout.rda"), compress=TRUE);	
		save(AUC.flat, AUC.hier, file=paste0(perf.dir, "AUC.", norm.type,".", flat.file, ".hierScores.",heuristic.name,".holdout.rda"), compress=TRUE);	
		save(PXR.flat, PXR.hier, file=paste0(perf.dir, "PXR.", norm.type,".", flat.file, ".hierScores.",heuristic.name,".holdout.rda"), compress=TRUE);	
		save(FMM.flat, FMM.hier, file=paste0(perf.dir, "FMM.", norm.type,".", flat.file, ".hierScores.",heuristic.name,".holdout.rda"), compress=TRUE);
		save(PRC.flat, PRC.hier, file=paste0(perf.dir, "PRC.", norm.type,".", flat.file, ".hierScores.",heuristic.name,".holdout.rda"), compress=TRUE);
	}
}
