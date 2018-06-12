##*********##
## TPR-DAG ##
##*********##

#' @name TPR-DAG 
#' @aliases tpr.threshold
#' @aliases tpr.threshold.free
#' @aliases tpr.weighted.threshold.free
#' @aliases tpr.weighted.threshold
#' @title TPR-DAG variants
#' @description Different variants of the TPR-DAG algorithm are implemented. In their more general form the TPR-DAG algorithms adopt
#' a two step learnig strategy:
#' \enumerate{
#'	\item in the first step they compute a \emph{per-level bottom-up} visit from the leaves to the root to propagate positive predictions across the hierarchy;
#'	\item in the second step they compute a \emph{per-level top-down} visit from the root to the leaves in order to assure the hierarchical 
#' 	consistency of the predictions
#' }
#' @details The \emph{vanilla} TPR-DAG adopts a per-level bottom-up traversal of the DAG to correct the flat predictions \eqn{\hat{y}_i}:
#' \deqn{
#' 	\bar{y}_i := \frac{1}{1 + |\phi_i|} (\hat{y}_i + \sum_{j \in \phi_i} \bar{y}_j)
#' }
#' where \eqn{\phi_i} are the positive children of \eqn{i}.
#' Different strategies to select the positive children \eqn{\phi_i} can be applied:
#' \enumerate{
#' 	\item \strong{Threshold-Free} strategy: the positive nodes are those children that can increment the score of the node \eqn{i}, that is those nodes 
#' 	that achieve a score higher than that of their parents:
#' 	\deqn{
#' 		\phi_i := \{ j \in child(i) | \bar{y}_j > \hat{y}_i \}
#' 	}
#' 	\item \strong{Threshold} strategy: the positive children are selected on the basis of a threshold that can ben selected in two different ways:
#' 	\enumerate{
#' 		\item for each node a constant threshold \eqn{\bar{t}} is a priori selected:
#'		\deqn{
#'			\phi_i := \{ j \in child(i) | \bar{y}_j > \bar{t} \}
#'		}
#' 		For instance if the predictions represent probabilities it could be meaningful to a priori select \eqn{\bar{t}=0.5}.
#' 		\item the threshold is selected to maximize some performance metric \eqn{\mathcal{M}} estimated on the training data, as for instance
#' 		the F-score or the AUPRC. In other words the threshold is selected to maximize some measure of accuracy of the predictions 
#' 		\eqn{\mathcal{M}(j,t)} on the training data for the class \eqn{j} with respect to the threshold \eqn{t}. 
#' 		The corresponding set of positives \eqn{\forall i \in V} is:
#' 		\deqn{
#' 			\phi_i := \{ j \in child(i) | \bar{y}_j > t_j^*,  t_j^* = \arg \max_{t} \mathcal{M}(j,t) \}
#' 		}
#' 		For instance \eqn{t_j^*} can be selected from a set of \eqn{t \in (0,1)} through internal cross-validation techniques.
#'	}
#' }
#' The weighted TPR-DAG version can be designed by adding a weight \eqn{w \in [0,1]} to balance between the 
#' contribution of the node \eqn{i} and that of its positive children \eqn{\phi}, through their convex combination:
#' \deqn{
#' 	\bar{y}_i := w \hat{y}_i + \frac{(1 - w)}{|\phi_i|} \sum_{j \in \phi_i} \bar{y}_j
#' }
#' If \eqn{w=1} no weight is attributed to the children and the TPR-DAG reduces to the HTD-DAG algorithm, since in this
#' way only the prediction for node \eqn{i} is used in the bottom-up step of the algorithm. If \eqn{w=0} only the predictors 
#' associated to the children nodes vote to predict node \eqn{i}. In the intermediate cases we attribute more importance to the predictor for the
#' node \eqn{i} or to its children depending on the values of \eqn{w}.
#' @param S a named flat scores matrix with examples on rows and classes on columns
#' @param g a graph of class \code{graphNEL}. It represents the hierarchy of the classes
#' @param root name of the class that it is on the top-level of the hierarchy (def. \code{root="00"})
#' @param t threshold for the choice of positive children (def. \code{t=0.5})
#' @param w weight to balance between the contribution of the node \eqn{i} and that of its positive children
#' @return a named matrix with the scores of the classes corrected according to the TPR-DAG algorithm.
#' @examples
#' data(graph);
#' data(scores);
#' data(labels);
#' root <- root.node(g);
#' S.tprTF <- tpr.threshold.free(S,g,root);
#' S.tprT <- tpr.threshold(S,g,root,t=0.5);
#' S.tprW <- tpr.weighted.threshold.free(S,g,root,w=0.5);
#' S.tprWT <- tpr.weighted.threshold(S,g,root,w=0.5, t=0.5);

#' @rdname TPR-DAG
#' @export 
tpr.threshold <- function(S, g, root="00", t=0.5){
	levels <- graph.levels(g,root);
	if(!(root %in% colnames(S))){
		max.score <- max(S);
		z <- rep(max.score,nrow(S));
		S <- cbind(z,S);
		colnames(S)[1] <- root;
	}
	# bottom-up visit
	chd.bup <- get.children.bottom.up(g,levels);
	for(i in 1:length(chd.bup)){
		if(length(chd.bup[[i]])!=0){
			parent <- S[,names(chd.bup[i])];
			children <- as.matrix(S[,chd.bup[[i]]]);
			# colnames(children) <- chd.bup[[i]]
			for(j in 1:length(parent)){
				child.set <- children[j,] > t;    # positive children selection
				child.pos <- children[j,][child.set];
				parent[j] <- (parent[j] + sum(child.pos))/(1+length(child.pos));  # flat scores correction
			}
			S[,names(chd.bup[i])] <- parent;
		}
	}
	# top-down visit
	par.tod <- get.parents.top.down(g,levels,root)
	for(i in 1:length(par.tod)){
		child <- S[,names(par.tod[i])];
		parents <- as.matrix(S[,par.tod[[i]]]);
		# colnames(parents) <- par.tod[[i]]
		# Note: the version with an apply and an ifelse statement is slower ...
		for(j in 1:length(child)){
			x <- min(parents[j,]);
			if(x < child[j]){
				child[j] <- x;    # hierarchical correction
			}
		}
		S[,names(par.tod[i])] <- child;
	}
	S <- S[,-which(colnames(S)==root)];
	return(S);
}

#' @rdname TPR-DAG
#' @export 
tpr.threshold.free <- function(S, g, root="00"){
	levels <- graph.levels(g,root)
	if(!(root %in% colnames(S))){
		max.score <- max(S);
		z <- rep(max.score,nrow(S));
		S <- cbind(z,S);
		colnames(S)[1] <- root;
	}
	# bottom-up visit
	chd.bup <- get.children.bottom.up(g,levels);
	for(i in 1:length(chd.bup)){
		if(length(chd.bup[[i]])!=0){
			parent <- S[,names(chd.bup[i])];
			children <- as.matrix(S[,chd.bup[[i]]]);
			# colnames(children) <- chd.bup[[i]]
			for(j in 1:length(parent)){
				child.set <- children[j,] > parent[j]; # positive children selection
				child.pos <- children[j,][child.set];
				parent[j] <- (parent[j] + sum(child.pos))/(1+length(child.pos));  # flat score correction
			}
			S[,names(chd.bup[i])] <- parent;
		} 
	}
	# top-down visit
	par.tod <- get.parents.top.down(g,levels,root);
	for(i in 1:length(par.tod)){
		child <- S[,names(par.tod[i])];
		parents <- as.matrix(S[,par.tod[[i]]]);
		# colnames(parents) <- par.tod[[i]]
		# Note: the version with an apply and an ifelse statement is slower ...
		for(j in 1:length(child)){
			x <- min(parents[j,]);
			if(x < child[j]){
				child[j] <- x;   # hierarchical correction
			}
		}
		S[,names(par.tod[i])] <- child;
	}
	S <- S[,-which(colnames(S)==root)];
	return(S);
}

#' @rdname TPR-DAG 
#' @export 
tpr.weighted.threshold.free <- function(S, g, root="00", w=0.5){
	levels <- graph.levels(g,root);
	if(!(root %in% colnames(S))){
		max.score <- max(S);
		z <- rep(max.score,nrow(S));
		S <- cbind(z,S);
		colnames(S)[1] <- root;
	}
	# bottom-up visit
	chd.bup <- get.children.bottom.up(g,levels);
	for(i in 1:length(chd.bup)){
		if(length(chd.bup[[i]])!=0){
			parent <- S[,names(chd.bup[i])];
			children <- as.matrix(S[,chd.bup[[i]]]);
			# colnames(children) <- chd.bup[[i]]
			for(j in 1:length(parent)){
				child.set <- children[j,] > parent[j];    # positive children selection
				child.pos <- children[j,][child.set];
				if(length(child.pos)!=0){
					parent[j] <- w*parent[j] + (1-w)*sum(child.pos)/length(child.pos);  # flat score correction
				}
			}
			S[,names(chd.bup[i])] <- parent;
		}
	}
	# top-down visit
	par.tod <- get.parents.top.down(g,levels,root)
	for(i in 1:length(par.tod)){
		child <- S[,names(par.tod[i])];
		parents <- as.matrix(S[,par.tod[[i]]]);
		# colnames(parents) <- par.tod[[i]]
		# Note: the version with an apply and an ifelse statement is slower ...
		for(j in 1:length(child)){
			x <- min(parents[j,]);
			if(x < child[j]){
				child[j] <- x;    # hierarchical correction
			}
		}
		S[,names(par.tod[i])] <- child;
	}
	S <- S[,-which(colnames(S)==root)];
	return(S);
}

#' @rdname TPR-DAG 
#' @export 
tpr.weighted.threshold <- function(S, g, root="00", t=0.5, w=0.5){
	levels <- graph.levels(g,root)
	if(!(root %in% colnames(S))){
		max.score <- max(S);
		z <- rep(max.score,nrow(S));
		S <- cbind(z,S);
		colnames(S)[1] <- root;
	}
	# bottom-up visit
	chd.bup <- get.children.bottom.up(g,levels);
	for(i in 1:length(chd.bup)){
		if(length(chd.bup[[i]])!=0){
			parent <- S[,names(chd.bup[i])];
			children <- as.matrix(S[,chd.bup[[i]]]);
			# colnames(children) <- chd.bup[[i]]
			for(j in 1:length(parent)){
				child.set <- children[j,] > t;    # positive children selection
				child.pos <- children[j,][child.set];
				if(length(child.pos)!=0){
					parent[j] <- w*parent[j] + (1-w)*sum(child.pos)/length(child.pos);  # flat score prediction
				}
			}
			S[,names(chd.bup[i])] <- parent;
		}
	}
	# top-down visit grafo
	par.tod <- get.parents.top.down(g,levels,root)
	for(i in 1:length(par.tod)){
		child <- S[,names(par.tod[i])];
		parents <- as.matrix(S[,par.tod[[i]]]);
		# colnames(parents) <- par.tod[[i]]
		# Note: the version with an apply and an ifelse statement is slower ...
		for(j in 1:length(child)){
			x <- min(parents[j,]);
			if(x < child[j]){
				child[j] <- x;    # hierarchical correction
			}
		}
		S[,names(par.tod[i])] <- child;
	}
	S <- S[,-which(colnames(S)==root)];
	return(S);
}

##************##
## DO TPR-DAG ##
##************##

#' @name TPR.DAG.CV
#' @seealso \code{\link{TPR-DAG}}
#' @aliases Do.tpr.threshold.cv
#' @aliases Do.tpr.weighted.threshold.free.cv
#' @aliases Do.tpr.weighted.threshold.cv
#' @title TPR-DAG cross-validation 
#' @description High level function to correct the computed scores in a hierarchy according to one of the TPR-DAG variants 
#' performing cross-validation experiments
#' @details These high level functions perform a classical cross-validation procedure to find the best threshold maximizing on F-score measure.
#' @param threshold range of threshold values to be tested in order to find the best threshold (def: \code{from:0.1}, \code{to:0.9}, \code{by:0.1}).
#' The denser the range is, the higher the probability to find the best theshold is, but obviously the execution time will be increasing.
#' @param weight range of weight values to be tested in order to find the best weight (def: \code{from:0.1}, \code{to:0.9}, \code{by:0.1}).
#' The denser the range is, the higher the probability to find the best weight is, but obviously the execution time will be increasing.
#' @param kk number of folds of the cross validation (\code{def: kk=5});
#' @param seed intialization seed for the random generator to create folds. If \code{NULL} (def.) no initialization is performed
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
#' threshold <- seq(from=0.1, to=0.5, by=0.1);
#' Do.tpr.threshold.cv(threshold=threshold, kk=5, seed=23, norm=FALSE, norm.type= "MaxNorm", 
#' flat.file=flat.file, ann.file=ann.file, dag.file=dag.file, flat.dir=flat.dir, ann.dir=ann.dir, 
#' dag.dir=dag.dir, flat.norm.dir=flat.norm.dir, n.round=3, f.criterion ="F", 
#' hierScore.dir=hierScore.dir, perf.dir=perf.dir);

#' @rdname TPR.DAG.CV
#' @export
Do.tpr.threshold.cv <- function(threshold=seq(from=0.1, to=0.9, by=0.1), kk=5, seed=NULL, norm=TRUE, norm.type= "NONE", flat.file=flat.file, 
	ann.file=ann.file, dag.file=dag.file, flat.dir=flat.dir, ann.dir=ann.dir, dag.dir=dag.dir, flat.norm.dir=NULL, n.round=3, f.criterion ="F", 
	hierScore.dir=hierScore.dir, perf.dir=perf.dir){
	
	if(kk==1) 
		stop("Smallest number of folds to define test and training set is 2. Set kk larger or equal to 2");

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

	## loading annotation matrix
	ann.path <- paste0(ann.dir, ann.file,".rda");
	ann <- get(load(ann.path));
	gc();

	## removing root node from annotation table 
	ann <- ann[,-which(colnames(ann)==root)];

	## FLAT AUC computed by precrec package
	AUC.flat <- AUROC.single.over.classes(ann, S); 
	
	# FLAT PxRs computed by PerfMeas pacakge
	PXR.flat <- precision.at.multiple.recall.level.over.classes(ann, S);

	## F.measure: Computing Flat Examples-Measures 
	FMM.flat <- find.best.f(ann, S, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=TRUE);

	## FLAT PRC computed by precrec package (more precise and accurate than PerfMeas)
	PRC.flat <- AUPRC.single.over.classes(ann, S); #saving PRC result in the same format of package PerfMeas

	## Hierarchical Correction #########################

	## splitting data in k unstratified fold
	folds <- do.unstratified.cv.data(S, kk=kk, seed=seed); 
	
	## storing the best F.max and the best threshold for each of k training set
	training.top.Fmax <- vector(mode="list", length=kk);
	names(training.top.Fmax) <- paste0(rep("fold",kk), 1:kk);

	## storing all the best average protein-centric measures (average and per.example) for each of k test set 
	best.avg.meas.test <- vector(mode="list", length=kk);
	names(best.avg.meas.test) <- paste0(rep("fold",kk), 1:kk);
	FMM.per.example <- c();

	## storing average macro AUC for each of k test set 
	AUC.average.test <- vector(mode="list", length=kk);
	names(AUC.average.test) <- paste0(rep("fold",kk), 1:kk);

	## storing per class macro AUC for each of k test set 
	AUC.class.test <- vector(mode="list", length=kk);
	names(AUC.class.test) <- paste0(rep("fold",kk), 1:kk);

	## storing average PxR for each of k test set 
	PXR.average.test <- vector(mode="list", length=kk);
	names(PXR.average.test) <- paste0(rep("fold",kk), 1:kk);

	## storing per class PxR for each of k test set 
	PXR.class.test <- vector(mode="list", length=kk);
	names(PXR.class.test) <- paste0(rep("fold",kk), 1:kk);

	## storing average macro AUC for each of k test set 
	PRC.average.test <- vector(mode="list", length=kk);
	names(PRC.average.test) <- paste0(rep("fold",kk), 1:kk);

	## storing per class macro AUC for each of k test set 
	PRC.class.test <- vector(mode="list", length=kk);
	names(PRC.class.test) <- paste0(rep("fold",kk), 1:kk);

	## variable for hosting the k sub-matrix assembled 
	S.hier <- c();

	## Let's start k-fold crossing validation for choosing best threshold maximizing on F.max...
	for(k in 1:kk){
		test <- S[folds[[k]],];		# test set: 1 out of k folds (e.g.: if k=5, testing set is 1/5 of data)
		training <- S[!rownames(S) %in% rownames(test),];	# training set: (k-1) out of k folds (e.g.: if k=5, training set is 4/5 of data)
		
		top.Fmax <- 0;
		best.Fmaxt <- 0;
		for(t in threshold){
			pred.training <- tpr.threshold(training, g, root=root, t=t);	## training set hierarchical correction...
			target.training <- ann[rownames(pred.training),colnames(pred.training)];
			training.avg.meas <- find.best.f(target.training, pred.training, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=FALSE); 
			training.Fmax <- training.avg.meas[4];	## F.max maximization...
			if(training.Fmax > top.Fmax){
				top.Fmax <- training.Fmax;
				best.Fmaxt <- t;
				training.top.Fmax[[k]] <- c(best.Fmax=top.Fmax, best.thres=best.Fmaxt);
				cat("training fold:",k, "better Fmax.avg found:",top.Fmax, "best threshold:",best.Fmaxt, sep="\t", "\n");
			}
		}
		pred.test <- tpr.threshold(test, g, root=root, t=training.top.Fmax[[k]][2]);
		target.test <- ann[rownames(pred.test),colnames(pred.test)];
		test.avg.meas <- find.best.f(target.test, pred.test, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=TRUE);	
		best.avg.meas.test[[k]] <- test.avg.meas$average;
		FMM.per.example <- rbind(FMM.per.example,test.avg.meas$per.example);

		## AUC (average and per.class) computed with PerfMeas package
		AUC.average.test[[k]] <- AUROC.single.over.classes(target.test,pred.test)$average;	## storing average AUC of each k testing set
		AUC.class.test[[k]] <- AUROC.single.over.classes(target.test,pred.test)$per.class;	## storing AUC per.class of each k testing set

		## PxR at fixed recall levels (average and per.class) computed with PerfMeas package
		PXR.average.test[[k]] <- precision.at.multiple.recall.level.over.classes(target.test,pred.test)$avgPXR;	## storing average PxR of each k testing set
		PXR.class.test[[k]] <- precision.at.multiple.recall.level.over.classes(target.test,pred.test)$PXR;	## storing PxR per.class of each k testing set

		## PRC (average and per.class) computed with precrec package
		PRC.average.test[[k]] <- AUPRC.single.over.classes(target.test,pred.test)$average;	## storing average PRC of each k testing set
		PRC.class.test[[k]] <- AUPRC.single.over.classes(target.test,pred.test)$per.class;	## storing PRC per.class of each k testing set

		## assembling of all the hierarchical scores of each k sub-matrix to build the full matrix 
		S.hier <- rbind(S.hier, pred.test);
	}
	## put the rows (i.e. genes) of assembled k sub-matrix in the same order of the full beginning matrix
	S.hier <- S.hier[rownames(S),];

	## remove no longer useful variables..
	rm(S, folds, pred.test); gc();

	## Averaging Performances Measures across k testing set ###################
	## averaging protein-centric measures across k testing sets
	F.cv <- Reduce("+", best.avg.meas.test)/kk;
	F.meas <-  apply(FMM.per.example,2,mean);
	names(F.meas)[4] <- "avF";
	F.max <- 2*(F.meas[["P"]] * F.meas[["R"]])/((F.meas[["P"]] + F.meas[["R"]]));
	names(F.max) <- "F";
	FMM.avg <- append(F.meas, F.max, after=3);
	FMM.avg <- append(FMM.avg,F.cv["T"]);
	FMM.hier <- list(average=FMM.avg, per.example=FMM.per.example);

	## averaging AUC (average and per.class) across k testing sets 
	AUC.average.tpr.over.test <- Reduce("+", AUC.average.test)/kk;
	AUC.class.tpr.over.test <- Reduce("+", AUC.class.test)/kk;
	AUC.hier <- list(average=AUC.average.tpr.over.test, per.class=AUC.class.tpr.over.test);

	## averaging PxR (average and per.class) across k testing sets 
	PXR.average.tpr.over.test <- Reduce("+", PXR.average.test)/kk;
	PXR.class.tpr.over.test <- Reduce("+", PXR.class.test)/kk;
	PXR.hier <- list(avgPXR=PXR.average.tpr.over.test, PXR=PXR.class.tpr.over.test);

	## averaging PRC (average and per.class) across k testing sets 
	PRC.average.tpr.over.test <- Reduce("+", PRC.average.test)/kk;
	PRC.class.tpr.over.test <- Reduce("+", PRC.class.test)/kk;
	PRC.hier <- list(average=PRC.average.tpr.over.test, per.class=PRC.class.tpr.over.test);

	## Storing Results #########
	if(norm){
		save(S.hier, file=paste0(hierScore.dir, flat.file, ".hierScores.tprT.rda"), compress=TRUE);
		save(AUC.flat, AUC.hier, file=paste0(perf.dir, "AUC.", flat.file, ".hierScores.tprT.rda"), compress=TRUE);
		save(PXR.flat, PXR.hier, file=paste0(perf.dir, "PXR.", flat.file, ".hierScores.tprT.rda"), compress=TRUE);
		save(FMM.flat, FMM.hier, file=paste0(perf.dir, "FMM.", flat.file, ".hierScores.tprT.rda"), compress=TRUE);
		save(PRC.flat, PRC.hier, file=paste0(perf.dir, "PRC.", flat.file, ".hierScores.tprT.rda"), compress=TRUE);
	}else{
		save(S.hier, file=paste0(hierScore.dir, norm.type,".", flat.file, ".hierScores.tprT.rda"), compress=TRUE);
		save(AUC.flat, AUC.hier, file=paste0(perf.dir, "AUC.", norm.type, ".", flat.file, ".hierScores.tprT.rda"), compress=TRUE);
		save(PXR.flat, PXR.hier, file=paste0(perf.dir, "PXR.", norm.type, ".", flat.file, ".hierScores.tprT.rda"), compress=TRUE);
		save(FMM.flat, FMM.hier, file=paste0(perf.dir, "FMM.", norm.type, ".", flat.file, ".hierScores.tprT.rda"), compress=TRUE);
		save(PRC.flat, PRC.hier, file=paste0(perf.dir, "PRC.", norm.type, ".", flat.file, ".hierScores.tprT.rda"), compress=TRUE);
	}
}

#' @name TPR.DAG.HOLDOUT
#' @seealso \code{\link{TPR-DAG}}
#' @aliases Do.tpr.threshold.holdout
#' @aliases Do.tpr.weighted.threshold.free.holdout
#' @aliases Do.tpr.weighted.threshold.holdout
#' @title TPR-DAG holdout 
#' @description High level function to correct the computed scores in a hierarchy according to one of the TPR-DAG variants 
#' performing a classical holdout procedure.
#' @details These high level functions perform a classical holdout procedure to find the best threshold maximizing on F.score measure.
#' @param threshold range of threshold values to be tested in order to find the best threshold (def: \code{from:0.1}, \code{to:0.9}, \code{by:0.1}).
#' The denser the range is, the higher the probability to find the best theshold is, but obviously the execution time will be increasing.
#' @param weight range of weight values to be tested in order to find the best weight (def: \code{from:0.1}, \code{to:0.9}, \code{by:0.1}).
#' The denser the range is, the higher the probability to find the best weight is, but obviously the execution time will be increasing.
#' @param kk number of folds of the cross validation (\code{def: kk=5});
#' @param seed intialization seed for the random generator to create folds. If \code{NULL} (def.) no initialization is performed
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
#' threshold <- seq(from=0.1, to=0.3, by=0.1);
#' Do.tpr.threshold.holdout(threshold=threshold, kk=5, seed=23, norm=FALSE, norm.type="MaxNorm", 
#' flat.file=flat.file, ann.file=ann.file, dag.file=dag.file, ind.test.set=ind.test.set, 
#' ind.dir=ind.dir, flat.dir=flat.dir, ann.dir=ann.dir, dag.dir=dag.dir, flat.norm.dir=flat.norm.dir, 
#' n.round=3, f.criterion ="F", hierScore.dir=hierScore.dir, perf.dir=perf.dir);

#' @rdname TPR.DAG.HOLDOUT
#' @export
Do.tpr.threshold.holdout <- function(threshold=seq(from=0.1, to=0.9, by=0.1), kk=5, seed=NULL, norm=TRUE, norm.type="NONE", 
	flat.file=flat.file, ann.file=ann.file, dag.file=dag.file, ind.test.set=ind.test.set, ind.dir=ind.dir,flat.dir=flat.dir, 
	ann.dir=ann.dir,dag.dir=dag.dir, flat.norm.dir=NULL, n.round=3,	f.criterion ="F", hierScore.dir=hierScore.dir, perf.dir=perf.dir){
	
	## Loading Data ############
	## loading examples indices of the test set
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

	## scores flat matrix shrinked to test and training test respectively
	S.test <- S[ind.test,];
	S.training <- S[-ind.test,];
	rm(S);

	## loading annotation matrix
	ann.path <- paste0(ann.dir, ann.file,".rda");
	ann <- get(load(ann.path));
	gc();

	## removing root node from annotation table and shrinking annotation table to tet set
	ann.test <- ann[ind.test,-which(colnames(ann)==root)];

	## FLAT AUC computed by precrec package
	AUC.flat <- AUROC.single.over.classes(ann.test, S.test); 

	# FLAT PxRs computed by PerfMeas pacakge
	PXR.flat <- precision.at.multiple.recall.level.over.classes(ann.test, S.test);

	## F.measure: Computing Flat Examples-Measures 
	FMM.flat <- find.best.f(ann.test, S.test, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=TRUE);

	## FLAT PRC computed by precrec package (more precise and accurate than PerfMeas)
	PRC.flat <- AUPRC.single.over.classes(ann.test, S.test); 

	## Hierarchical Correction #########################
	## TPR-Threshold correction using k-fold cross-validation to find best threshold on training set
	folds <- do.unstratified.cv.data(S.training, kk=kk, seed=seed);
			
	## Let's start k-fold crossing validation for choosing best threshold maximizing on the basis of F.max value...
	for(k in 1:kk){
		training <- S.training[folds[[k]],];
		top.Fmax <- 0;
		best.Fmaxt <- 0;
		for(t in threshold){
			pred.training <- tpr.threshold(training, g, root=root, t=t);	## training set hierarchical correction...
			target.training <- ann[rownames(pred.training),colnames(pred.training)];
			training.avg.meas <- find.best.f(target.training, pred.training, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=FALSE);
			training.Fmax <- training.avg.meas[4];	## F.max maximization...
			if(training.Fmax > top.Fmax){
				top.Fmax <- training.Fmax;
				best.Fmaxt <- t;
				cat("training fold:",k, "better Fmax.avg found:",top.Fmax, "best threshold:",best.Fmaxt, sep="\t", "\n");
			}
		}
	}
	S.test <- tpr.threshold(S.test, g, root=root, t=best.Fmaxt);	## testing set hierarchical correction..
		
	## AUC (average and per.class) computed with PerfMeas package on the test set
	AUC.hier <- AUROC.single.over.classes(ann.test, S.test);

	## PxR at fixed recall levels (average and per.class) computed with PerfMeas package on the test set
	PXR.hier <- precision.at.multiple.recall.level.over.classes(ann.test, S.test);

	## F.measure: Computing Hierarchical Examples-Measures 
	FMM.hier <- find.best.f(ann.test, S.test, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=TRUE);	

	## Hierarchical PRC (average and per.class) computed by precrec package
	PRC.hier <- AUPRC.single.over.classes(ann.test, S.test); 

	## storing the hierarchical matrix
	S.hier <- S.test;
	rm(S.test, S.training); gc();

	## Storing Results #########
	if(norm){
		save(S.hier, file=paste0(hierScore.dir, flat.file, ".hierScores.tprT.holdout.rda"), compress=TRUE);
		save(AUC.flat, AUC.hier, file=paste0(perf.dir, "AUC.", flat.file, ".hierScores.tprT.holdout.rda"), compress=TRUE);
		save(PXR.flat, PXR.hier, file=paste0(perf.dir, "PXR.", flat.file, ".hierScores.tprT.holdout.rda"), compress=TRUE);
		save(FMM.flat, FMM.hier, file=paste0(perf.dir, "FMM.", flat.file, ".hierScores.tprT.holdout.rda"), compress=TRUE);
		save(PRC.flat, PRC.hier, file=paste0(perf.dir, "PRC.", flat.file, ".hierScores.tprT.holdout.rda"), compress=TRUE);
	}else{
		save(S.hier, file=paste0(hierScore.dir, norm.type,".", flat.file, ".hierScores.tprT.holdout.rda"), compress=TRUE);
		save(AUC.flat, AUC.hier, file=paste0(perf.dir, "AUC.", norm.type, ".", flat.file, ".hierScores.tprT.holdout.rda"), compress=TRUE);
		save(PXR.flat, PXR.hier, file=paste0(perf.dir, "PXR.", norm.type, ".", flat.file, ".hierScores.tprT.holdout.rda"), compress=TRUE);
		save(FMM.flat, FMM.hier, file=paste0(perf.dir, "FMM.", norm.type, ".", flat.file, ".hierScores.tprT.holdout.rda"), compress=TRUE);
		save(PRC.flat, PRC.hier, file=paste0(perf.dir, "PRC.", norm.type, ".", flat.file, ".hierScores.tprT.holdout.rda"), compress=TRUE);
	}
}

#' @rdname TPR.DAG.CV
#' @export
Do.tpr.weighted.threshold.free.cv <- function(weight=seq(from=0.1, to=1, by=0.1), kk=5, seed=NULL, norm=TRUE, norm.type= "NONE", 
	flat.file=flat.file, ann.file=ann.file, dag.file=dag.file, flat.dir=flat.dir, ann.dir=ann.dir,dag.dir=dag.dir, 
	flat.norm.dir=NULL, n.round=3, f.criterion ="F", hierScore.dir=hierScore.dir, perf.dir=perf.dir){
	
	if(kk==1) 
		stop("Smallest number of folds to define test and training set is 2. Set kk larger or equal to 2");

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

	## loading annotation matrix
	ann.path <- paste0(ann.dir, ann.file,".rda");
	ann <- get(load(ann.path));
	gc();

	## removing root node from annotation table 
	ann <- ann[,-which(colnames(ann)==root)];

	## FLAT AUC computed by precrec package
	AUC.flat <- AUROC.single.over.classes(ann, S);

	# FLAT PxRs computed by PerfMeas pacakge
	PXR.flat <- precision.at.multiple.recall.level.over.classes(ann, S);

	## F.measure: Computing Flat Examples-Measures 
	FMM.flat <- find.best.f(ann, S, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=TRUE); 

	## FLAT PRC computed by precrec package (more precise and accurate than PerfMeas)
	PRC.flat <- AUPRC.single.over.classes(ann, S); gc();

	## Hierarchical Correction #########################
	## TPR-Weighted-Threshold-Free correction using k-fold cross-validation to find best weight 

	## splitting data in k stratified fold
	folds <- do.unstratified.cv.data(S, kk=kk, seed=seed); 

	## storing the best F.max and the best threshold for each of k training set
	training.top.Fmax <- vector(mode="list", length=kk);
	names(training.top.Fmax) <- paste0(rep("fold",kk), 1:kk);

	## storing all the best protein-centric measures for each of k test set 
	best.avg.meas.test <- vector(mode="list", length=kk);
	names(best.avg.meas.test) <- paste0(rep("fold",kk), 1:kk);
	FMM.per.example <- c();

	## storing average macro AUC for each of k test set 
	AUC.average.test <- vector(mode="list", length=kk);
	names(AUC.average.test) <- paste0(rep("fold",kk), 1:kk);

	## storing per class macro AUC for each of k test set 
	AUC.class.test <- vector(mode="list", length=kk);
	names(AUC.class.test) <- paste0(rep("fold",kk), 1:kk);

	## storing average PxR for each of k test set 
	PXR.average.test <- vector(mode="list", length=kk);
	names(PXR.average.test) <- paste0(rep("fold",kk), 1:kk);

	## storing per class PxR for each of k test set 
	PXR.class.test <- vector(mode="list", length=kk);
	names(PXR.class.test) <- paste0(rep("fold",kk), 1:kk);

	## storing average macro AUC for each of k test set 
	PRC.average.test <- vector(mode="list", length=kk);
	names(PRC.average.test) <- paste0(rep("fold",kk), 1:kk);

	## storing per class macro AUC for each of k test set 
	PRC.class.test <- vector(mode="list", length=kk);
	names(PRC.class.test) <- paste0(rep("fold",kk), 1:kk);

	## variable for hosting the k sub-matrix assembled 
	S.hier <- c();

	## Let's start k-fold crossing validation for choosing best threshold maximizing on F.max...
	for(k in 1:kk){
		test <- S[folds[[k]],];		# test set: 1 out of k folds (e.g.: if k=5, testing set is 1/5 of data)
		training <- S[!rownames(S) %in% rownames(test),];		# training set: (k-1) out of k folds (e.g.: if k=5, training set is 4/5 of data)
		gc();

		top.Fmax <- 0;
		best.Fmaxw <- 0;
		for(w in weight){
			pred.training <- tpr.weighted.threshold.free(training, g, root=root, w=w);	## training set hierarchical correction...
			target.training <- ann[rownames(pred.training),colnames(pred.training)];
			training.avg.meas <- find.best.f(target.training, pred.training, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=FALSE);	 
			training.Fmax <- training.avg.meas[4];	## F.max maximization...
			if(training.Fmax > top.Fmax){
				top.Fmax <- training.Fmax;
				best.Fmaxw <- w;
				training.top.Fmax[[k]] <- c(best.Fmax=top.Fmax, best.weigh=best.Fmaxw);
				cat("training fold:",k, "better Fmax.avg found:",top.Fmax, "best weight:",best.Fmaxw, sep="\t", "\n");
			}
		}
		pred.test <- tpr.weighted.threshold.free(test, g, root=root, w=training.top.Fmax[[k]][2]);	## testing set hierarchical correction...
		target.test <- ann[rownames(pred.test),colnames(pred.test)];
		test.avg.meas <- find.best.f(target.test, pred.test, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=TRUE);	
		best.avg.meas.test[[k]] <- test.avg.meas$average;
		FMM.per.example <- rbind(FMM.per.example,test.avg.meas$per.example);

		## Hierarchical AUC (average and per.class) computed with precrec package
		AUC.average.test[[k]] <- AUROC.single.over.classes(target.test,pred.test)$average;	## storing average AUC of each k testing set
		AUC.class.test[[k]] <- AUROC.single.over.classes(target.test,pred.test)$per.class;	## storing AUC per.class of each k testing set

		## Hierarchical PxR at fixed recall levels (average and per.class) computed with PerfMeas package
		PXR.average.test[[k]] <- precision.at.multiple.recall.level.over.classes(target.test,pred.test)$avgPXR;	## storing average PxR of each k testing set
		PXR.class.test[[k]] <- precision.at.multiple.recall.level.over.classes(target.test,pred.test)$PXR;	## storing PxR per.class of each k testing set

		## PRC (average and per.class) computed with precrec package
		PRC.average.test[[k]] <- AUPRC.single.over.classes(target.test,pred.test)$average;	## storing average PRC of each k testing set
		PRC.class.test[[k]] <- AUPRC.single.over.classes(target.test,pred.test)$per.class;	## storing PRC per.class of each k testing set

		## assembling of all the hierarchical scores of each k sub-matrix to build the full matrix 
		S.hier <- rbind(S.hier, pred.test);
	}
	## put the rows (i.e. genes) of assembled k sub-matrix in the same order of the full beginning matrix
	S.hier <- S.hier[rownames(S),];

	## remove no longer useful variables..
	rm(S, folds, pred.test); gc();

	## Averaging Performances Measures across k testing set ###################
	## averaging protein-centric measures across k testing sets
	F.cv <- Reduce("+", best.avg.meas.test)/kk;
	F.meas <-  apply(FMM.per.example,2,mean);
	names(F.meas)[4] <- "avF";
	F.max <- 2*(F.meas[["P"]] * F.meas[["R"]])/((F.meas[["P"]] + F.meas[["R"]]));
	names(F.max) <- "F";
	FMM.avg <- append(F.meas, F.max, after=3);
	FMM.avg <- append(FMM.avg,F.cv["T"]);
	FMM.hier <- list(average=FMM.avg, per.example=FMM.per.example);

	## averaging macro-AUC (average and per.class) across k testing sets 
	AUC.average.tpr.over.test <- Reduce("+", AUC.average.test)/kk;
	AUC.class.tpr.over.test <- Reduce("+", AUC.class.test)/kk;
	AUC.hier <- list(average=AUC.average.tpr.over.test, per.class=AUC.class.tpr.over.test);

	## averaging PxR (average and per.class) across k testing sets 
	PXR.average.tpr.over.test <- Reduce("+", PXR.average.test)/kk;
	PXR.class.tpr.over.test <- Reduce("+", PXR.class.test)/kk;
	PXR.hier <- list(avgPXR=PXR.average.tpr.over.test, PXR=PXR.class.tpr.over.test);

	## averaging PRC (average and per.class) across k testing sets 
	PRC.average.tpr.over.test <- Reduce("+", PRC.average.test)/kk;
	PRC.class.tpr.over.test <- Reduce("+", PRC.class.test)/kk;
	PRC.hier <- list(average=PRC.average.tpr.over.test, per.class=PRC.class.tpr.over.test);

	## Storing Results #########
	if(norm){
		save(S.hier, file=paste0(hierScore.dir, flat.file, ".hierScores.tprW.rda"), compress=TRUE);
		save(AUC.flat, AUC.hier, file=paste0(perf.dir, "AUC.", flat.file, ".hierScores.tprW.rda"), compress=TRUE);
		save(PXR.flat, PXR.hier, file=paste0(perf.dir, "PXR.", flat.file, ".hierScores.tprW.rda"), compress=TRUE);
		save(FMM.flat, FMM.hier, file=paste0(perf.dir, "FMM.", flat.file, ".hierScores.tprW.rda"), compress=TRUE);
		save(PRC.flat, PRC.hier, file=paste0(perf.dir, "PRC.", flat.file, ".hierScores.tprW.rda"), compress=TRUE);
	}else{
		save(S.hier, file=paste0(hierScore.dir, norm.type,".", flat.file, ".hierScores.tprW.rda"), compress=TRUE);
		save(AUC.flat, AUC.hier, file=paste0(perf.dir, "AUC.", norm.type, ".", flat.file, ".hierScores.tprW.rda"), compress=TRUE);
		save(PXR.flat, PXR.hier, file=paste0(perf.dir, "PXR.", norm.type, ".", flat.file, ".hierScores.tprW.rda"), compress=TRUE);
		save(FMM.flat, FMM.hier, file=paste0(perf.dir, "FMM.", norm.type, ".", flat.file, ".hierScores.tprW.rda"), compress=TRUE);
		save(PRC.flat, PRC.hier, file=paste0(perf.dir, "PRC.", norm.type, ".", flat.file, ".hierScores.tprW.rda"), compress=TRUE);
	}
}

#' @rdname TPR.DAG.HOLDOUT
#' @export
Do.tpr.weighted.threshold.free.holdout <- function(weight=seq(from=0.1, to=1, by=0.1), kk=5, seed=NULL, norm=TRUE, norm.type= "NONE", 
	flat.file=flat.file, ann.file=ann.file, dag.file=dag.file, ind.test.set=ind.test.set, ind.dir=ind.dir, flat.dir=flat.dir, ann.dir=ann.dir, 
	dag.dir=dag.dir, flat.norm.dir=NULL, n.round=3, f.criterion ="F", hierScore.dir=hierScore.dir, perf.dir=perf.dir){
	
	## Loading Data ############
	## loading examples indices of the test set
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

	## scores flat matrix shrinked to test and training test respectively
	S.test <- S[ind.test,];
	S.training <- S[-ind.test,];
	rm(S);

	## loading annotation matrix
	ann.path <- paste0(ann.dir, ann.file,".rda");
	ann <- get(load(ann.path));
	gc();

	## removing root node from annotation table and shrinking annotation table to tet set
	ann.test <- ann[ind.test,-which(colnames(ann)==root)];

	## FLAT AUC computed by precrec package
	AUC.flat <- AUROC.single.over.classes(ann.test, S.test); 
	
	# FLAT PxRs computed by PerfMeas pacakge
	PXR.flat <- precision.at.multiple.recall.level.over.classes(ann.test, S.test);

	## F.measure: Computing Flat Examples-Measures 
	FMM.flat <- find.best.f(ann.test, S.test, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=TRUE);

	## FLAT PRC computed by precrec package (more precise and accurate than PerfMeas)
	PRC.flat <- AUPRC.single.over.classes(ann.test, S.test); 
	
	## Hierarchical Correction #########################
	folds <- do.unstratified.cv.data(S.training, kk=kk, seed=seed);
			
	## Let's start k-fold crossing validation for choosing best threshold maximizing on the basis of F.max value...
	for(k in 1:kk){
		training <- S.training[folds[[k]],];	
		top.Fmax <- 0;
		best.Fmaxw <- 0;
		for(w in weight){
			pred.training <- tpr.weighted.threshold.free(training, g, root=root, w=w);	## training set hierarchical correction...
			target.training <- ann[rownames(pred.training),colnames(pred.training)];
			training.avg.meas <- find.best.f(target.training, pred.training, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=FALSE); 
			training.Fmax <- training.avg.meas[4];	## F.max maximization...
			if(training.Fmax > top.Fmax){
				top.Fmax <- training.Fmax;
				best.Fmaxw <- w;
				cat("training fold:",k, "better Fmax.avg found:",top.Fmax, "best threshold:",best.Fmaxw, sep="\t", "\n");
			}
		}
	}
	S.test <- tpr.weighted.threshold.free(S.test, g, root=root, w=best.Fmaxw);	## testing set hierarchical correction..
		
	## AUC (average and per.class) computed with PerfMeas package on the test set
	AUC.hier <- AUROC.single.over.classes(ann.test, S.test);	

	## PxR at fixed recall levels (average and per.class) computed with PerfMeas package on the test set
	PXR.hier <- precision.at.multiple.recall.level.over.classes(ann.test, S.test);

	## F.measure: Computing Hierarchical Examples-Measures 
	FMM.hier <- find.best.f(ann.test, S.test, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=TRUE);

	## Hierarchical PRC (average and per.class) computed by precrec package
	PRC.hier <- AUPRC.single.over.classes(ann.test, S.test); 

	## storing the hierarchical matrix
	S.hier <- S.test;
	rm(S.test, S.training);

	## Storing Results #########
	if(norm){
		save(S.hier, file=paste0(hierScore.dir, flat.file, ".hierScores.tprW.holdout.rda"), compress=TRUE);
		save(AUC.flat, AUC.hier, file=paste0(perf.dir, "AUC.", flat.file, ".hierScores.tprW.holdout.rda"), compress=TRUE);
		save(PXR.flat, PXR.hier, file=paste0(perf.dir, "PXR.", flat.file, ".hierScores.tprW.holdout.rda"), compress=TRUE);
		save(FMM.flat, FMM.hier, file=paste0(perf.dir, "FMM.", flat.file, ".hierScores.tprW.holdout.rda"), compress=TRUE);
		save(PRC.flat, PRC.hier, file=paste0(perf.dir, "PRC.", flat.file, ".hierScores.tprW.holdout.rda"), compress=TRUE);
	}else{
		save(S.hier, file=paste0(hierScore.dir, norm.type,".", flat.file, ".hierScores.tprW.holdout.rda"), compress=TRUE);
		save(AUC.flat, AUC.hier, file=paste0(perf.dir, "AUC.", norm.type, ".", flat.file, ".hierScores.tprW.holdout.rda"), compress=TRUE);
		save(PXR.flat, PXR.hier, file=paste0(perf.dir, "PXR.", norm.type, ".", flat.file, ".hierScores.tprW.holdout.rda"), compress=TRUE);
		save(FMM.flat, FMM.hier, file=paste0(perf.dir, "FMM.", norm.type, ".", flat.file, ".hierScores.tprW.holdout.rda"), compress=TRUE);
		save(PRC.flat, PRC.hier, file=paste0(perf.dir, "PRC.", norm.type, ".", flat.file, ".hierScores.tprW.holdout.rda"), compress=TRUE);
	}
}

#' @rdname TPR.DAG.CV
#' @export
Do.tpr.weighted.threshold.cv <- function(threshold=seq(from=0.1, to=0.9, by=0.1), weight=seq(from=0.1, to=1, by=0.1), kk=5, seed=NULL, 
	norm=TRUE, norm.type= "NONE", flat.file=flat.file, ann.file=ann.file, dag.file=dag.file, flat.dir=flat.dir, ann.dir=ann.dir, 
	dag.dir=dag.dir, flat.norm.dir=NULL, n.round=3,  f.criterion ="F", hierScore.dir=hierScore.dir, perf.dir=perf.dir){
	
	if(kk==1) 
		stop("Smallest number of folds to define test and training set is 2. Set kk larger or equal to 2");

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

	## loading annotation matrix
	ann.path <- paste0(ann.dir, ann.file,".rda");
	ann <- get(load(ann.path));
	gc();

	## removing root node from annotation table 
	ann <- ann[,-which(colnames(ann)==root)];

	## FLAT AUC computed by precrec package
	AUC.flat <- AUROC.single.over.classes(ann, S); 

	# FLAT PxRs computed by PerfMeas pacakge
	PXR.flat <- precision.at.multiple.recall.level.over.classes(ann, S);

	## F.measure: Computing Flat Examples-Measures 
	FMM.flat <- find.best.f(ann, S, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=TRUE);

	## FLAT PRC computed by precrec package (more precise and accurate than PerfMeas)
	PRC.flat <- AUPRC.single.over.classes(ann, S); 

	## Hierarchical Correction #########################
	## TPR-Weighted-Threshold correction using k-fold cross-validation to find best threshold and weight 

	## splitting data in k unstratified folds 
	folds <- do.unstratified.cv.data(S, kk=kk, seed=seed); 
	
	## storing the best F.max and the best threshold and weight for each of k training set
	training.top.Fmax <- vector(mode="list", length=kk);
	names(training.top.Fmax) <- paste0(rep("fold",kk), 1:kk);

	## storing all the best protein-centric measures for each of k test set 
	best.avg.meas.test <- vector(mode="list", length=kk);
	names(best.avg.meas.test) <- paste0(rep("fold",kk), 1:kk);
	FMM.per.example <- c();

	## storing average macro AUC for each of k test set 
	AUC.average.test <- vector(mode="list", length=kk);
	names(AUC.average.test) <- paste0(rep("fold",kk), 1:kk);

	## storing per class macro AUC for each of k test set 
	AUC.class.test <- vector(mode="list", length=kk);
	names(AUC.class.test) <- paste0(rep("fold",kk), 1:kk);

	## storing average PxR for each of k test set 
	PXR.average.test <- vector(mode="list", length=kk);
	names(PXR.average.test) <- paste0(rep("fold",kk), 1:kk);

	## storing per class PxR for each of k test set 
	PXR.class.test <- vector(mode="list", length=kk);
	names(PXR.class.test) <- paste0(rep("fold",kk), 1:kk);

	## storing average macro PRCC for each of k test set 
	PRC.average.test <- vector(mode="list", length=kk);
	names(PRC.average.test) <- paste0(rep("fold",kk), 1:kk);

	## storing per class macro PRC for each of k test set 
	PRC.class.test <- vector(mode="list", length=kk);
	names(PRC.class.test) <- paste0(rep("fold",kk), 1:kk);

	## variable for hosting the k sub-matrix assembled 
	S.hier <- c();

	## Let's start k-fold crossing validation for choosing best threshold and weight maximizing on F.max...
	for(k in 1:kk){
		test <- S[folds[[k]],];	# test set: 1 out of k folds (e.g.: if k=5, testing set is 1/5 of data)
		training <- S[!rownames(S) %in% rownames(test),];		# training set: (k-1) out of k folds (e.g.: if k=5, training set is 4/5 of data)
		gc();

		top.Fmax <- 0;
		best.Fmaxt <- 0;
		best.Fmaxw <- 0;
		for(t in threshold){
			for(w in weight){
				pred.training <- tpr.weighted.threshold(training, g, root=root, w=w, t=t);	## training set hierarchical correction...
				target.training <- ann[rownames(pred.training),colnames(pred.training)];
				training.avg.meas <- find.best.f(target.training, pred.training, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=FALSE);	 
				training.Fmax <- training.avg.meas[4];	## F.max maximization...
				if(training.Fmax > top.Fmax){
					top.Fmax <- training.Fmax;
					best.Fmaxt <- t;
					best.Fmaxw <- w;
					training.top.Fmax[[k]] <- c(best.Fmax=top.Fmax, best.thres=best.Fmaxt, best.weight=best.Fmaxw);
					cat("training fold:",k, "better Fmax.avg found:",top.Fmax, "best threshold:",best.Fmaxt, "best weight:", best.Fmaxw,sep="\t", "\n");
				}
			}
		}
		pred.test <- tpr.weighted.threshold(test, g, root=root, t=training.top.Fmax[[k]][2], w=training.top.Fmax[[k]][3]);
		target.test <- ann[rownames(pred.test),colnames(pred.test)];
		test.avg.meas <- find.best.f(target.test, pred.test, n.round=n.round, f.criterion =f.criterion, verbose=FALSE, b.per.example=TRUE);	
		best.avg.meas.test[[k]] <- test.avg.meas$average;
		FMM.per.example <- rbind(FMM.per.example,test.avg.meas$per.example);

		## Hierarchical AUC (average and per.class) computed with PerfMeas package
		AUC.average.test[[k]] <- AUROC.single.over.classes(target.test,pred.test)$average;	## storing average AUC of each k testing set
		AUC.class.test[[k]] <- AUROC.single.over.classes(target.test,pred.test)$per.class;	## storing AUC per.class of each k testing set

		## Hierarchical PxR at fixed recall levels (average and per.class) computed with PerfMeas package
		PXR.average.test[[k]] <- precision.at.multiple.recall.level.over.classes(target.test,pred.test)$avgPXR;	## storing average PxR of each k testing set
		PXR.class.test[[k]] <- precision.at.multiple.recall.level.over.classes(target.test,pred.test)$PXR;	## storing PxR per.class of each k testing set

		## PRC (average and per.class) computed with precrec package
		PRC.average.test[[k]] <- AUPRC.single.over.classes(target.test,pred.test)$average;	## storing average AUC of each k testing set
		PRC.class.test[[k]] <- AUPRC.single.over.classes(target.test,pred.test)$per.class;	## storing AUC per.class of each k testing set

		## assembling of all the hierarchical scores of each k sub-matrix to build the full matrix 
		S.hier <- rbind(S.hier, pred.test);
	}
	## put the rows (i.e. genes) of assembled k sub-matrix in the same order of the full beginning matrix
	S.hier <- S.hier[rownames(S),];

	## remove no longer useful variables..
	rm(S, folds, pred.test);

	## Averaging Performances Measures across k testing set ###################
	## averaging protein-centric measures across k testing sets
	F.cv <- Reduce("+", best.avg.meas.test)/kk;
	F.meas <-  apply(FMM.per.example,2,mean);
	names(F.meas)[4] <- "avF";
	F.max <- 2*(F.meas[["P"]] * F.meas[["R"]])/((F.meas[["P"]] + F.meas[["R"]]));
	names(F.max) <- "F";
	FMM.avg <- append(F.meas, F.max, after=3);
	FMM.avg <- append(FMM.avg,F.cv["T"]);
	FMM.hier <- list(average=FMM.avg, per.example=FMM.per.example);

	## averaging macro-AUC (average and per.class) across k testing sets 
	AUC.average.tpr.over.test <- Reduce("+", AUC.average.test)/kk;
	AUC.class.tpr.over.test <- Reduce("+", AUC.class.test)/kk;
	AUC.hier <- list(average=AUC.average.tpr.over.test, per.class=AUC.class.tpr.over.test);

	## averaging PxR (average and per.class) across k testing sets 
	PXR.average.tpr.over.test <- Reduce("+", PXR.average.test)/kk;
	PXR.class.tpr.over.test <- Reduce("+", PXR.class.test)/kk;
	PXR.hier <- list(avgPXR=PXR.average.tpr.over.test, PXR=PXR.class.tpr.over.test);

	## averaging PRC (average and per.class) across k testing sets 
	PRC.average.tpr.over.test <- Reduce("+", PRC.average.test)/kk;
	PRC.class.tpr.over.test <- Reduce("+", PRC.class.test)/kk;
	PRC.hier <- list(average=PRC.average.tpr.over.test, per.class=PRC.class.tpr.over.test);

	## Storing Results #########
	if(norm){
		save(S.hier, file=paste0(hierScore.dir, flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
		save(AUC.flat, AUC.hier, file=paste0(perf.dir, "AUC.", flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
		save(PXR.flat, PXR.hier, file=paste0(perf.dir, "PXR.", flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
		save(FMM.flat, FMM.hier, file=paste0(perf.dir, "FMM.", flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
		save(PRC.flat, PRC.hier, file=paste0(perf.dir, "PRC.", flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
	}else{
		save(S.hier, file=paste0(hierScore.dir, norm.type,".", flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
		save(AUC.flat, AUC.hier, file=paste0(perf.dir, "AUC.", norm.type, ".", flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
		save(PXR.flat, PXR.hier, file=paste0(perf.dir, "PXR.", norm.type, ".", flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
		save(FMM.flat, FMM.hier, file=paste0(perf.dir, "FMM.", norm.type, ".", flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
		save(PRC.flat, PRC.hier, file=paste0(perf.dir, "PRC.", norm.type, ".", flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
	}
}

#' @rdname TPR.DAG.HOLDOUT
#' @export
Do.tpr.weighted.threshold.holdout <- function(threshold=seq(from=0.1, to=0.9, by=0.1), weight=seq(from=0.1, to=1, by=0.1), kk=5, seed=NULL, 
	norm=TRUE, norm.type= "NONE", flat.file=flat.file, ann.file=ann.file, dag.file=dag.file, ind.test.set=ind.test.set, ind.dir=ind.dir,
	flat.dir=flat.dir, ann.dir=ann.dir,dag.dir=dag.dir, flat.norm.dir=NULL, n.round=3, f.criterion ="F", hierScore.dir=hierScore.dir, perf.dir=perf.dir){

	## Loading Data ############
	## loading examples indices of the test set
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

	## scores flat matrix shrinked to test and training test respectively
	S.test <- S[ind.test,];
	S.training <- S[-ind.test,];
	rm(S);

	## loading annotation matrix
	ann.path <- paste0(ann.dir, ann.file,".rda");
	ann <- get(load(ann.path));
	gc();

	## removing root node from annotation table and shrinking annotation table to tet set
	ann.test <- ann[ind.test,-which(colnames(ann)==root)];

	## FLAT AUC computed by precrec package
	AUC.flat <- AUROC.single.over.classes(ann.test, S.test); 

	# FLAT PxRs computed by PerfMeas pacakge
	PXR.flat <- precision.at.multiple.recall.level.over.classes(ann.test, S.test);

	## F.measure: Computing Flat Examples-Measures 
	FMM.flat <- find.best.f(ann.test, S.test, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=TRUE);

	## FLAT PRC computed by precrec package (more precise and accurate than PerfMeas)
	PRC.flat <- AUPRC.single.over.classes(ann.test, S.test); 

	## Hierarchical Correction #########################
	folds <- do.unstratified.cv.data(S.training, kk=kk, seed=seed);
			
	## Let's start k-fold crossing validation for choosing best threshold maximizing on the basis of F.max value...
	for(k in 1:kk){
		training <- S.training[folds[[k]],];	
		
		top.Fmax <- 0;
		best.Fmaxt <- 0;
		best.Fmaxw <- 0;
		for(t in threshold){
			for(w in weight){
				pred.training <- tpr.weighted.threshold(training, g, root=root, w=w, t=t); ## training set hierarchical correction...
				target.training <- ann[rownames(pred.training),colnames(pred.training)];
				training.avg.meas <- find.best.f(target.training, pred.training, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=FALSE); 
				training.Fmax <- training.avg.meas[4];	## F.max maximization...
				if(training.Fmax > top.Fmax){
					top.Fmax <- training.Fmax;
					best.Fmaxt <- t;
					best.Fmaxw <- w;
					cat("training fold:",k, "better Fmax.avg found:",top.Fmax, "best threshold:",best.Fmaxt, "best weight:",best.Fmaxw, sep="\t","\n");
				}
			}
		}
	}
	S.test <- tpr.weighted.threshold(S.test, g, root=root, t=best.Fmaxt, w=best.Fmaxw);	## testing set hierarchical correction..
		
	## AUC (average and per.class) computed with PerfMeas package on the test set
	AUC.hier <- AUROC.single.over.classes(ann.test, S.test);

	## PxR at fixed recall levels (average and per.class) computed with PerfMeas package on the test set
	PXR.hier <- precision.at.multiple.recall.level.over.classes(ann.test, S.test);

	## F.measure: Computing Hierarchical Examples-Measures 
	FMM.hier <- find.best.f(ann.test, S.test, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=TRUE);	

	## Hierarchical PRC (average and per.class) computed by precrec package
	PRC.hier <- AUPRC.single.over.classes(ann.test, S.test); 

	## storing the hierarchical matrix
	S.hier <- S.test;
	rm(S.test, S.training);

	## Storing Results #########
	if(norm){
		save(S.hier, file=paste0(hierScore.dir, flat.file, ".hierScores.tprWT.holdout.rda"), compress=TRUE);
		save(AUC.flat, AUC.hier, file=paste0(perf.dir, "AUC.", flat.file, ".hierScores.tprWT.holdout.rda"), compress=TRUE);
		save(PXR.flat, PXR.hier, file=paste0(perf.dir, "PXR.", flat.file, ".hierScores.tprWT.holdout.rda"), compress=TRUE);
		save(FMM.flat, FMM.hier, file=paste0(perf.dir, "FMM.", flat.file, ".hierScores.tprWT.holdout.rda"), compress=TRUE);
		save(PRC.flat, PRC.hier, file=paste0(perf.dir, "PRC.", flat.file, ".hierScores.tprWT.holdout.rda"), compress=TRUE);
	}else{
		save(S.hier, file=paste0(hierScore.dir, norm.type,".", flat.file, ".hierScores.tprWT.holdout.rda"), compress=TRUE);
		save(AUC.flat, AUC.hier, file=paste0(perf.dir, "AUC.", norm.type, ".", flat.file, ".hierScores.tprWT.holdout.rda"), compress=TRUE);
		save(PXR.flat, PXR.hier, file=paste0(perf.dir, "PXR.", norm.type, ".", flat.file, ".hierScores.tprWT.holdout.rda"), compress=TRUE);
		save(FMM.flat, FMM.hier, file=paste0(perf.dir, "FMM.", norm.type, ".", flat.file, ".hierScores.tprWT.holdout.rda"), compress=TRUE);
		save(PRC.flat, PRC.hier, file=paste0(perf.dir, "PRC.", norm.type, ".", flat.file, ".hierScores.tprWT.holdout.rda"), compress=TRUE);
	}
}

#' @title TPR-DAG vanilla
#' @seealso \code{\link{TPR-DAG}}
#' @description High level function to correct the computed scores in a hierarchy according to the TPR-DAG vanilla algorithm
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
#' Do.tpr.threshold.free(norm=FALSE, norm.type= "MaxNorm", flat.file=flat.file, 
#' ann.file=ann.file, dag.file=dag.file, flat.dir=flat.dir, ann.dir=ann.dir, 
#' dag.dir=dag.dir, flat.norm.dir=flat.norm.dir, n.round=3, f.criterion ="F", 
#' hierScore.dir=hierScore.dir, perf.dir=perf.dir);
Do.tpr.threshold.free <- function(norm=TRUE, norm.type= "NONE", flat.file=flat.file, ann.file=ann.file, dag.file=dag.file, flat.dir=flat.dir, 
	ann.dir=ann.dir, dag.dir=dag.dir, flat.norm.dir=NULL, n.round=3, f.criterion="F", hierScore.dir=hierScore.dir, perf.dir=perf.dir){
	
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

	## loading annotation matrix
	ann.path <- paste0(ann.dir, ann.file,".rda");
	ann <- get(load(ann.path));
	gc();

	## removing root node from annotation table 
	ann <- ann[,-which(colnames(ann)==root)];

	## Computing Flat Perfromances
	## FLAT AUC computed by precrec package
	AUC.flat <- AUROC.single.over.classes(ann, S); 
	
	# FLAT PxRs computed by PerfMeas pacakge
	PXR.flat <- precision.at.multiple.recall.level.over.classes(ann, S);

	## F.measure: Computing Flat Examples-Measures 
	FMM.flat <- find.best.f(ann, S, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=TRUE);

	## FLAT PRC computed by precrec package (more precise and accurate than PerfMeas)
	PRC.flat <- AUPRC.single.over.classes(ann, S); 

	## Hierarchical Correction #########################
	## TPR-Threshold-free correction 
	S <- tpr.threshold.free(S, g, root=root);

	## Computing Hierarchical Performances
	## Hierarchical AUC (average and per.class) computed by PerfMeas package
	AUC.hier <- AUROC.single.over.classes(ann, S);
	
	## Hierarchical PxR at fixed recall levels (average and per.class) computed by PerfMeas package
	PXR.hier <- precision.at.multiple.recall.level.over.classes(ann, S);

	## Computing Hierarchical Examples-Measures 
	FMM.hier <- find.best.f(ann, S, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=TRUE);

	## Hierarchical PRC (average and per.class) computed by precrec package
	PRC.hier <- AUPRC.single.over.classes(ann, S); 

	## storing the hierarchical matrix
	S.hier <- S;
	rm(S);

	## Storing Results #########
	if(norm){
		save(S.hier, file=paste0(hierScore.dir, flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
		save(AUC.flat, AUC.hier, file=paste0(perf.dir, "AUC.", flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
		save(PXR.flat, PXR.hier, file=paste0(perf.dir, "PXR.", flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
		save(FMM.flat, FMM.hier, file=paste0(perf.dir, "FMM.", flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
		save(PRC.flat, PRC.hier, file=paste0(perf.dir, "PRC.", flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
	}else{
		save(S.hier, file=paste0(hierScore.dir, norm.type,".", flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
		save(AUC.flat, AUC.hier, file=paste0(perf.dir, "AUC.", norm.type, ".", flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
		save(PXR.flat, PXR.hier, file=paste0(perf.dir, "PXR.", norm.type, ".", flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
		save(FMM.flat, FMM.hier, file=paste0(perf.dir, "FMM.", norm.type, ".", flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
		save(PRC.flat, PRC.hier, file=paste0(perf.dir, "PRC.", norm.type, ".", flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
	}
}

#' @title TPR-DAG vanilla holdout
#' @seealso \code{\link{TPR-DAG}}
#' @description High level function to correct the computed scores in a hierarchy according to the TPR-DAG vanilla algorithm 
#' applying a classical holdout procedure
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
#' Do.tpr.threshold.free.holdout(norm=FALSE, norm.type= "MaxNorm", flat.file=flat.file, 
#' ann.file=ann.file, dag.file=dag.file, ind.test.set=ind.test.set, ind.dir=ind.dir, 
#' flat.dir=flat.dir, ann.dir=ann.dir, dag.dir=dag.dir, flat.norm.dir=flat.norm.dir, 
#' n.round=3, f.criterion ="F", hierScore.dir=hierScore.dir, perf.dir=perf.dir);
Do.tpr.threshold.free.holdout <- function(norm=TRUE, norm.type= "NONE", flat.file=flat.file, ann.file=ann.file, dag.file=dag.file, flat.dir=flat.dir, 
	ann.dir=ann.dir, dag.dir=dag.dir, flat.norm.dir=NULL, ind.test.set=ind.test.set, ind.dir=ind.dir, n.round=3, f.criterion ="F",
	hierScore.dir=hierScore.dir, perf.dir=perf.dir){

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
	if(norm==TRUE){
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
	gc();

	## removing root node from annotation table 
	ann <- ann[ind.test,-which(colnames(ann)==root)];

	## Computing Flat Perfromances
	## FLAT AUC computed by PerfMeas package
	AUC.flat <- AUROC.single.over.classes(ann, S); 

	# FLAT PxRs computed by PerfMeas pacakge
	PXR.flat <- precision.at.multiple.recall.level.over.classes(ann, S);

	## F.measure: Computing Flat Examples-Measures 
	FMM.flat <- find.best.f(ann, S, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=TRUE);

	## FLAT PRC computed by precrec package (more precise and accurate than PerfMeas)
	PRC.flat <- AUPRC.single.over.classes(ann, S); 

	## Hierarchical Correction #########################
	## TPR-Threshold-free correction 
	S <- tpr.threshold.free(S, g, root=root);

	## Computing Hierarchical Performances
	## Hierarchical AUC (average and per.class) computed by PerfMeas package
	AUC.hier <- AUROC.single.over.classes(ann, S);
	
	## Hierarchical PxR at fixed recall levels (average and per.class) computed by PerfMeas package
	PXR.hier <- precision.at.multiple.recall.level.over.classes(ann, S);

	## Computing Hierarchical Examples-Measures 
	FMM.hier <- find.best.f(ann, S, n.round=n.round, f.criterion =f.criterion, verbose=FALSE, b.per.example=TRUE);

	## Hierarchical PRC (average and per.class) computed by precrec package
	PRC.hier <- AUPRC.single.over.classes(ann, S); 

	## storing the hierarchical matrix
	S.hier <- S;
	rm(S);

	## Storing Results #########
	if(norm){
		save(S.hier, file=paste0(hierScore.dir, flat.file, ".hierScores.tprTF.holdout.rda"), compress=TRUE);
		save(AUC.flat, AUC.hier, file=paste0(perf.dir, "AUC.", flat.file, ".hierScores.tprTF.holdout.rda"), compress=TRUE);
		save(PXR.flat, PXR.hier, file=paste0(perf.dir, "PXR.", flat.file, ".hierScores.tprTF.holdout.rda"), compress=TRUE);
		save(FMM.flat, FMM.hier, file=paste0(perf.dir, "FMM.", flat.file, ".hierScores.tprTF.holdout.rda"), compress=TRUE);
		save(PRC.flat, PRC.hier, file=paste0(perf.dir, "PRC.", flat.file, ".hierScores.tprTF.holdout.rda"), compress=TRUE);
	}else{
		save(S.hier, file=paste0(hierScore.dir, norm.type,".", flat.file, ".hierScores.tprTF.holdout.rda"), compress=TRUE);
		save(AUC.flat, AUC.hier, file=paste0(perf.dir, "AUC.", norm.type, ".", flat.file, ".hierScores.tprTF.holdout.rda"), compress=TRUE);
		save(PXR.flat, PXR.hier, file=paste0(perf.dir, "PXR.", norm.type, ".", flat.file, ".hierScores.tprTF.holdout.rda"), compress=TRUE);
		save(FMM.flat, FMM.hier, file=paste0(perf.dir, "FMM.", norm.type, ".", flat.file, ".hierScores.tprTF.holdout.rda"), compress=TRUE);
		save(PRC.flat, PRC.hier, file=paste0(perf.dir, "PRC.", norm.type, ".", flat.file, ".hierScores.tprTF.holdout.rda"), compress=TRUE);
	}
}
