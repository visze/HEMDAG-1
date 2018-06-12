##******************##
## TPR-DAG VARIANTS ##
##******************##

#' @name TPR-DAG-variants
#' @seealso \code{\link{TPR-DAG}}, \code{\link{DESCENS}}, \code{\link{HTD-DAG}}
#' @title TPR-DAG Variants
#' @description \code{TPR-DAG} is a user-friendly function gathering all the hierarchical ensemble algorithms
#' @param S a named flat scores matrix with examples on rows and classes on columns
#' @param g a graph of class \code{graphNEL}. It represents the hierarchy of the classes
#' @param root name of the class that it is on the top-level of the hierarchy (def. \code{root="00"})
#' @param positives choice of the \emph{positive} nodes to be considered in the bottom-up strategy. Can be one of the following values:
#' \itemize{
#' 	\item \code{children}: for each node are considered its positive children (\code{def.});
#' 	\item \code{descendants}: for each node are considered its positive descendants;
#' }
#' @param bottomup strategy to enhance the flat predictions by propagating the positive predictions from leaves to root. 
#' It can be one of the following values:
#' \itemize{
#' 	\item \code{threshold.free}: positive nodes are selected on the basis of the \code{threshold.free} strategy (\code{def.});
#' 	\item \code{threshold}: positive nodes are selected on the basis of the \code{threshold} strategy;
#' 	\item \code{weighted.threshold.free}: positive nodes are selected on the basis of the \code{weighted.threshold.free} strategy;
#' 	\item \code{weighted.threshold}: positive nodes are selected on the basis of the \code{weighted.threshold} strategy;
#' 	\item \code{tau}: positive nodes are selected on the basis of the \code{tau} strategy. 
#'	NOTE: \code{tau} is only a \code{DESCENS} variants. If you use \code{tau} strategy you must set the parameter \code{positives} 
#' to \code{descendants};
#' }
#' @param t threshold for the choice of positive nodes (def. \code{t=0.5}). Set \code{t} only for the variants that requiring 
#' a threshold for the selection of the positive nodes, otherwise set \code{t} to zero
#' @param w weight to balance between the contribution of the node \eqn{i} and that of its positive nodes. Set \code{w} only for the
#' \emph{weighted} variants, otherwise set \code{w} to zero
#' @return a named matrix with the scores of the classes corrected according to the chosen algorithm
#' @export 
#' @examples
#' data(graph);
#' data(scores);
#' data(labels);
#' root <- root.node(g);
#' S.hier <- TPR.DAG(S, g, root, positives="children", bottomup="threshold.free", t=0, w=0);
TPR.DAG <- function(S, g, root="00", positives="children", bottomup="threshold.free", t=0, w=0){
	
	## Setting Check
	if(positives!="children" && positives!="descendants" || bottomup!="threshold" && bottomup!="threshold.free" && 
		bottomup!="weighted.threshold" && bottomup!="weighted.threshold.free" && bottomup!="tau"){
		stop("TPR.DAG: positives or bottomup value misspelled. Check them");
	}
	if(positives=="children" && bottomup=="tau"){
		stop("tau is a descendants variants. Please set positives to descendants");
	}
	if(bottomup=="threshold" || bottomup=="tau"){w<-0;}
	if(bottomup=="threshold.free"){t<-0; w<-0;}
	if(bottomup=="weighted.threshold.free"){t<-0;}

	## computing graph levels
	levels <- graph.levels(g,root);
	if(!(root %in% colnames(S))){
		max.score <- max(S);
		z <- rep(max.score,nrow(S));
		S <- cbind(z,S);
		colnames(S)[1] <- root;
	}
	## bottom-up visit: positive children selection
	if(positives=="children"){
		chd.bup <- get.children.bottom.up(g,levels);
		for(i in 1:length(chd.bup)){
			if(length(chd.bup[[i]])!=0){
				parent <- S[,names(chd.bup[i])];
				children <- as.matrix(S[,chd.bup[[i]]]);
				# colnames(children) <- chd.bup[[i]]
				for(j in 1:length(parent)){
					if(bottomup=="threshold"){
						child.set <- children[j,] > t;    # positive children selection
						child.pos <- children[j,][child.set];
						parent[j] <- (parent[j] + sum(child.pos))/(1+length(child.pos));  # flat scores correction
					}else if(bottomup=="threshold.free"){
						child.set <- children[j,] > parent[j]; # positive children selection
						child.pos <- children[j,][child.set];
						parent[j] <- (parent[j] + sum(child.pos))/(1+length(child.pos));  # flat score correction
					}else if(bottomup=="weighted.threshold.free"){
						child.set <- children[j,] > parent[j];    # positive children selection
						child.pos <- children[j,][child.set];
						if(length(child.pos)!=0){
							parent[j] <- w*parent[j] + (1-w)*sum(child.pos)/length(child.pos);  # flat score correction
						}
					}else if(bottomup=="weighted.threshold"){
						child.set <- children[j,] > t;    # positive children selection
						child.pos <- children[j,][child.set];
						if(length(child.pos)!=0){
							parent[j] <- w*parent[j] + (1-w)*sum(child.pos)/length(child.pos);  # flat score prediction
						}
					}
				}	
				S[,names(chd.bup[i])] <- parent;
			}
		}
	## bottom-up visit: positive descendants selection
	}else if(positives=="descendants"){
		if(bottomup=="tau"){
			chd.bup <- get.children.bottom.up(g,levels);
		}
		desc.bup <- build.descendants.bottom.up(g,levels);
		nodes <- names(desc.bup);
		for(i in 1:length(desc.bup)){
			if(length(desc.bup[[i]])!=1){
				node.curr <- nodes[i];
				parent <- S[,names(desc.bup[i])];
				tmp <- setdiff(desc.bup[[i]],node.curr);
				if(bottomup=="tau"){
					delta <- setdiff(tmp, chd.bup[[i]]);  # descendants without children 
					children <- as.matrix(S[,chd.bup[[i]]]);	# genes considering children node 
					desc <-  as.matrix(S[,delta]);		# genes considering descendants nodes without children
				}else{
					desc <- as.matrix(S[,tmp]);
				}
				for(j in 1:length(parent)){
					if(bottomup=="threshold"){
						desc.set <- desc[j,] > t;    # positive descendants selection
						desc.pos <- desc[j,][desc.set];
						parent[j] <- (parent[j] + sum(desc.pos))/(1+length(desc.pos));   # flat scores correction
					}else if(bottomup=="threshold.free"){
						desc.set <- desc[j,] > parent[j];	# positive descendants selection
						desc.pos <- desc[j,][desc.set];
						parent[j] <- (parent[j] + sum(desc.pos))/(1+length(desc.pos));   # flat scores correction
					}else if(bottomup=="weighted.threshold.free"){
						desc.set <- desc[j,] > parent[j];
						desc.pos <- desc[j,][desc.set];
						if(length(desc.pos)!=0){
							parent[j] <- w*parent[j] + (1-w)*sum(desc.pos)/length(desc.pos);  # flat scores correction
						}
					}else if(bottomup=="weighted.threshold"){
						desc.set <- desc[j,] > t;
						desc.pos <- desc[j,][desc.set];
						if(length(desc.pos)!=0){
							parent[j] <- w*parent[j] + (1-w)*sum(desc.pos)/length(desc.pos);  # flat scores correction
						}
					}else if(bottomup=="tau"){
						desc.set <- desc[j,] > parent[j];			# positive descendants (without children) selection
						desc.pos <- desc[j,][desc.set];
						child.set <- children[j,] > parent[j];  	# positive children selection
						child.pos <- children[j,][child.set];
						parent[j] <- t * ((parent[j] + sum(child.pos))/(1+length(child.pos))) + (1-t) * ((parent[j] + sum(desc.pos))/(1+length(desc.pos)));
					}
					S[,names(desc.bup[i])] <- parent;
				}
			}
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
				child[j] <- x;    # hierarchical correction
			}
		}
		S[,names(par.tod[i])] <- child;
	}
	S <- S[,-which(colnames(S)==root)];
	return(S);
}


#' @name TPR-DAG-cross-validation
#' @title TPR-DAG cross-validation experiments
#' @seealso \code{\link{TPR-DAG-variants}}
#' @description High level function to correct the computed scores in a hierarchy according to the chosen ensemble algorithm
#' through a k-cross-validation 
#' @details The variants choosing the positives nodes on the basis of a parameter are cross-validated by maximizing 
#' on F-measure (\code{\link{Multilabel.F.measure}})
#' @param threshold range of threshold values to be tested in order to find the best threshold (def: \code{from:0.1}, \code{to:0.9}, \code{by:0.1}).
#' The denser the range is, the higher the probability to find the best theshold is, but obviously the execution time will be higher.
#' Set the parameter \code{threshold} only for the variants that requiring a threshold for the positive nodes selection, 
#' otherwise set the parameter \code{threshold} to zero
#' @param weight range of weight values to be tested in order to find the best weight (def: \code{from:0.1}, \code{to:0.9}, \code{by:0.1}).
#' The denser the range is, the higher the probability to find the best theshold is, but obviously the execution time will be higher.
#' Set the parameter \code{weight} only for the \emph{weighted} variants, otherwise set the parameter \code{weight} to zero
#' @param kk number of folds of the cross validation (\code{def: kk=5});
#' @param seed intialization seed for the random generator to create folds. If \code{NULL} (def.) no initialization is performed
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
#' @param positives choice of the \emph{positive} nodes to be considered in the bottom-up strategy. Can be one of the following values:
#' \itemize{
#' 	\item \code{children}: for each node are considered its positive children (\code{def.});
#' 	\item \code{descendants}: for each node are considered its positive descendants;
#' }
#' @param bottomup strategy to enhance the flat predictions by propagating the positive predictions from leaves to root. 
#' It can be one of the following values:
#' \itemize{
#' 	\item \code{threshold.free}: positive nodes are selected on the basis of the \code{threshold.free} strategy (\code{def.});
#' 	\item \code{threshold}: positive nodes are selected on the basis of the \code{threshold} strategy;
#' 	\item \code{weighted.threshold.free}: positive nodes are selected on the basis of the \code{weighted.threshold.free} strategy;
#' 	\item \code{weighted.threshold}: positive nodes are selected on the basis of the \code{weighted.threshold} strategy;
#' 	\item \code{tau}: positive nodes are selected on the basis of the \code{tau} strategy. 
#'	NOTE: \code{tau} is only a \code{DESCENS} variants. If you use \code{tau} strategy you must set the parameter \code{positives} 
#' to \code{descendants};
#' }
#' @param flat.file name of the file containing the flat scores matrix to be normalized or already normalized (without rda extension)
#' @param ann.file name of the file containing the the label matrix of the examples (without rda extension)
#' @param dag.file name of the file containing the graph that represents the hierarchy of the classes (without rda extension)
#' @param flat.dir relative path where flat scores matrix is stored
#' @param ann.dir relative path where annotation matrix is stored
#' @param dag.dir relative path where graph is stored
#' @param flat.norm.dir relative path where flat normalized scores matrix must be stored. Use this parameter if and only if \code{norm} is
#' set to \code{FALSE}, otherwise set \code{flat.norm.dir} to \code{NULL} (def.)
#' @param n.round number of rounding digits to be applied to the hierarchical scores matrix (\code{def. 3}). 
#' It is used for choosing the best threshold on the basis of the best F-measure
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
#' threshold <- weight <- 0;
#' positives <- "children";
#' bottomup <- "threshold.free";
#' Do.TPR.DAG(threshold=threshold, weight=weight, kk=5, seed=23, norm=FALSE, 
#' norm.type="MaxNorm", positives=positives, bottomup=bottomup,
#' flat.file=flat.file, ann.file=ann.file, dag.file=dag.file, flat.dir=flat.dir, 
#' ann.dir=ann.dir, dag.dir=dag.dir, flat.norm.dir=flat.norm.dir, n.round=3, 
#' f.criterion="F", hierScore.dir=hierScore.dir, perf.dir=perf.dir);
Do.TPR.DAG <- function(threshold=seq(from=0.1, to=0.9, by=0.1), weight=seq(from=0.1, to=0.9, by=0.1), kk=5, 
	seed=NULL, norm=TRUE, norm.type=NULL, positives="children", bottomup="threshold.free",
	flat.file=flat.file, ann.file=ann.file, dag.file=dag.file, flat.dir=flat.dir, flat.norm.dir=NULL, 
	ann.dir=ann.dir, dag.dir=dag.dir, n.round=3, f.criterion="F", hierScore.dir=hierScore.dir, 
	perf.dir=perf.dir){
	
	## Setting Check
	if(positives!="children" && positives!="descendants" || bottomup!="threshold" && bottomup!="threshold.free" && 
		bottomup!="weighted.threshold" && bottomup!="weighted.threshold.free" && bottomup!="tau"){
		stop("TPR.DAG: positives or bottomup value misspelled");
	}
	if(positives=="children" && bottomup=="tau"){
		stop("tau is a descendants variants. Please set positives to descendants");
	}
	if(bottomup=="threshold" || bottomup=="tau"){weight<-0;}
	if(bottomup=="threshold.free"){threshold<-0; weight<-0;}
	if(bottomup=="weighted.threshold.free"){threshold<-0;}

	if(norm==FALSE && length(norm.type)==0){
		stop("If norm is set to FALSE, you need also to specify a normalization method among those available");
	}
	if(norm==TRUE && length(norm.type)!=0){
		stop("If norm is set to TRUE, the input flat matrix is already normalized. Set norm.type' to NULL (without quote)");
	}

	if(kk==1){
		stop("Smallest number of folds to define test and training set is 2. Set kk larger or equal to 2");
	}
	
	## loading dag
	dag.path <- paste0(dag.dir, dag.file,".rda");
	g <- get(load(dag.path));
	root <- root.node(g);
	
	## loading flat matrix
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
	if(bottomup=="threshold.free"){
		S <- TPR.DAG(S, g, root=root, positives=positives, bottomup=bottomup, t=0, w=0);
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
	}else{
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
					pred.training <- TPR.DAG(training, g, root=root, positives=positives, bottomup=bottomup, w=w, t=t);
					target.training <- ann[rownames(pred.training),colnames(pred.training)];
					training.avg.meas <- find.best.f(target.training, pred.training, n.round=n.round, f.criterion=f.criterion, 
						verbose=FALSE, b.per.example=FALSE);	 
					training.Fmax <- training.avg.meas[4];	## F.max maximization...
					if(training.Fmax > top.Fmax){
						top.Fmax <- training.Fmax;
						best.Fmaxt <- t;
						best.Fmaxw <- w;
						training.top.Fmax[[k]] <- c(best.Fmax=top.Fmax, best.thres=best.Fmaxt, best.weight=best.Fmaxw);
						if(bottomup=="threshold" || bottomup=="tau"){
							cat("training fold:",k, "better Fmax.avg found:",top.Fmax, "best threshold:",best.Fmaxt, sep="\t", "\n");
						}else if(bottomup=="weighted.threshold.free"){
							cat("training fold:",k, "better Fmax.avg found:",top.Fmax, "best weight:", best.Fmaxw, sep="\t", "\n");
						}
						else{
							cat("training fold:",k, "better Fmax.avg found:",top.Fmax, "best threshold:",best.Fmaxt, "best weight:", best.Fmaxw, sep="\t", "\n");
						}
					}
				}
			}
			pred.test <- TPR.DAG(test, g, root=root, positives=positives, bottomup=bottomup, t=training.top.Fmax[[k]][2], w=training.top.Fmax[[k]][3]);
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
		rm(S, folds, pred.test, test, training, target.test, target.training);

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
	}

	## Storing Results #########
	if(positives=="children" && bottomup=="threshold.free"){
		meth.name <- "tprTF";
	}else if(positives=="children" && bottomup=="threshold"){
		meth.name <- "tprT";
	}else if(positives=="children" && bottomup=="weighted.threshold.free"){
		meth.name <- "tprW";
	}else if(positives=="children" && bottomup=="weighted.threshold"){
		meth.name <- "tprWT";
	}
	if(positives=="descendants" && bottomup=="threshold.free"){
		meth.name <- "descensTF";
	}else if(positives=="descendants" && bottomup=="threshold"){
		meth.name <- "descensT";
	}else if(positives=="descendants" && bottomup=="weighted.threshold.free"){
		meth.name <- "descensW";
	}else if(positives=="descendants" && bottomup=="weighted.threshold"){
		meth.name <- "descensWT";
	}else if(positives=="descendants" && bottomup=="tau"){
		meth.name <- "descensTAU";
	}
	if(norm){
		save(S.hier, file=paste0(hierScore.dir, flat.file, ".hierScores.",meth.name,".rda"), compress=TRUE);
		save(AUC.flat, AUC.hier, file=paste0(perf.dir, "AUC.", flat.file, ".hierScores.",meth.name,".rda"), compress=TRUE);
		save(PXR.flat, PXR.hier, file=paste0(perf.dir, "PXR.", flat.file, ".hierScores.",meth.name,".rda"), compress=TRUE);
		save(FMM.flat, FMM.hier, file=paste0(perf.dir, "FMM.", flat.file, ".hierScores.",meth.name,".rda"), compress=TRUE);
		save(PRC.flat, PRC.hier, file=paste0(perf.dir, "PRC.", flat.file, ".hierScores.",meth.name,".rda"), compress=TRUE);
	}else{
		save(S.hier, file=paste0(hierScore.dir, norm.type,".", flat.file, ".hierScores.",meth.name,".rda"), compress=TRUE);
		save(AUC.flat, AUC.hier, file=paste0(perf.dir, "AUC.", norm.type, ".", flat.file, ".hierScores.",meth.name,".rda"), compress=TRUE);
		save(PXR.flat, PXR.hier, file=paste0(perf.dir, "PXR.", norm.type, ".", flat.file, ".hierScores.",meth.name,".rda"), compress=TRUE);
		save(FMM.flat, FMM.hier, file=paste0(perf.dir, "FMM.", norm.type, ".", flat.file, ".hierScores.",meth.name,".rda"), compress=TRUE);
		save(PRC.flat, PRC.hier, file=paste0(perf.dir, "PRC.", norm.type, ".", flat.file, ".hierScores.",meth.name,".rda"), compress=TRUE);
	}
}

#' @name TPR-DAG-holdout
#' @title TPR-DAG holdout experiments
#' @seealso \code{\link{TPR-DAG-variants}}
#' @description High level function to correct the computed scores in a hierarchy according to the chosen ensemble algorithm 
#' through an hold-out procedure
#' @details The variants choosing the positives nodes on the basis of a parameter are cross-validated by maximizing 
#' on F-measure (\code{\link{Multilabel.F.measure}})
#' @param threshold range of threshold values to be tested in order to find the best threshold (def: \code{from:0.1}, \code{to:0.9}, \code{by:0.1}).
#' The denser the range is, the higher the probability to find the best theshold is, but obviously the execution time will be higher.
#' Set the parameter \code{threshold} only for the variants that requiring a threshold for the positive nodes selection, 
#' otherwise set the parameter \code{threshold} to zero
#' @param weight range of weight values to be tested in order to find the best weight (def: \code{from:0.1}, \code{to:0.9}, \code{by:0.1}).
#' The denser the range is, the higher the probability to find the best theshold is, but obviously the execution time will be higher.
#' Set the parameter \code{weight} only for the \emph{weighted} variants, otherwise set the parameter \code{weight} to zero
#' @param kk number of folds of the cross validation (\code{def: kk=5});
#' @param seed intialization seed for the random generator to create folds. If \code{NULL} (def.) no initialization is performed
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
#' @param positives choice of the \emph{positive} nodes to be considered in the bottom-up strategy. Can be one of the following values:
#' \itemize{
#' 	\item \code{children}: for each node are considered its positive children (\code{def.});
#' 	\item \code{descendants}: for each node are considered its positive descendants;
#' }
#' @param bottomup strategy to enhance the flat predictions by propagating the positive predictions from leaves to root. 
#' It can be one of the following values:
#' \itemize{
#' 	\item \code{threshold.free}: positive nodes are selected on the basis of the \code{threshold.free} strategy (\code{def.});
#' 	\item \code{threshold}: positive nodes are selected on the basis of the \code{threshold} strategy;
#' 	\item \code{weighted.threshold.free}: positive nodes are selected on the basis of the \code{weighted.threshold.free} strategy;
#' 	\item \code{weighted.threshold}: positive nodes are selected on the basis of the \code{weighted.threshold} strategy;
#' 	\item \code{tau}: positive nodes are selected on the basis of the \code{tau} strategy. 
#'	NOTE: \code{tau} is only a \code{DESCENS} variants. If you use \code{tau} strategy you must set the parameter \code{positives} 
#' to \code{descendants};
#' }
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
#' dag.dir <- flat.dir <- flat.norm.dir <- ann.dir <- "data/";
#' ind.test.set <- "test.index";
#' dag.file <- "graph";
#' flat.file <- "scores";
#' ann.file <- "labels";
#' threshold <- weight <- 0;
#' positives <- "children";
#' bottomup <- "threshold.free";
#' Do.TPR.DAG.holdout(threshold=threshold, weight=weight, kk=5, seed=23, norm=FALSE, 
#' norm.type="MaxNorm", positives=positives, bottomup=bottomup,flat.file=flat.file, 
#' ann.file=ann.file, dag.file=dag.file, ind.test.set=ind.test.set, 
#' ind.dir=ind.dir, flat.dir=flat.dir, ann.dir=ann.dir, dag.dir=dag.dir, 
#' flat.norm.dir=flat.norm.dir, n.round=3, f.criterion="F", 
#' hierScore.dir=hierScore.dir, perf.dir=perf.dir);
Do.TPR.DAG.holdout <- function(threshold=seq(from=0.1, to=0.9, by=0.1), weight=seq(from=0.1, to=1, by=0.1), kk=5, 
	seed=NULL, norm=TRUE, norm.type=NULL, positives="children", bottomup="threshold.free", flat.file=flat.file, 
	ann.file=ann.file, dag.file=dag.file, ind.test.set=ind.test.set, ind.dir=ind.dir, flat.dir=flat.dir, 
	ann.dir=ann.dir, dag.dir=dag.dir, flat.norm.dir=NULL, n.round=3, f.criterion="F", 
	hierScore.dir=hierScore.dir, perf.dir=perf.dir){

	## Setting Check
	if(positives!="children" && positives!="descendants" || bottomup!="threshold" && bottomup!="threshold.free" && 
		bottomup!="weighted.threshold" && bottomup!="weighted.threshold.free" && bottomup!="tau"){
		stop("TPR.DAG: positives or bottomup value misspelled");
	}
	if(positives=="children" && bottomup=="tau"){
		stop("tau is a descendants variants. Please set positives to descendants");
	}
	if(bottomup=="threshold" || bottomup=="tau"){weight<-0;}
	if(bottomup=="threshold.free"){threshold<-0; weight<-0;}
	if(bottomup=="weighted.threshold.free"){threshold<-0;}

	if(norm==FALSE && length(norm.type)==0){
		stop("If norm is set to FALSE, you need also to specify a normalization method among those available");
	}
	if(norm==TRUE && length(norm.type)!=0){
		stop("If norm is set to TRUE, the input flat matrix is already normalized. Set norm.type' to NULL (without quote)");
	}

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

	## loading annotation matrix
	ann.path <- paste0(ann.dir, ann.file,".rda");
	ann <- get(load(ann.path));

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
	if(bottomup=="threshold.free"){
		S.test <- TPR.DAG(S.test, g, root=root, positives=positives, bottomup=bottomup, t=0, w=0);
		## Hierarchical AUC (average and per.class) computed by PerfMeas package
		AUC.hier <- AUROC.single.over.classes(ann.test, S.test);
		## Hierarchical PxR at fixed recall levels (average and per.class) computed by PerfMeas package
		PXR.hier <- precision.at.multiple.recall.level.over.classes(ann.test, S.test);
		## Computing Hierarchical Examples-Measures 
		FMM.hier <- find.best.f(ann.test, S.test, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=TRUE);
		## Hierarchical PRC (average and per.class) computed by precrec package
		PRC.hier <- AUPRC.single.over.classes(ann.test, S.test); 
		## storing the hierarchical matrix
		S.hier <- S;
		rm(S);
	}else{
		folds <- do.unstratified.cv.data(S.training, kk=kk, seed=seed);
			
		## Let's start k-fold crossing validation for choosing best threshold maximizing on the basis of F.max value...
		for(k in 1:kk){
			training <- S.training[folds[[k]],];
			
			top.Fmax <- 0;
			best.Fmaxt <- 0;
			best.Fmaxw <- 0;
			for(t in threshold){
				for(w in weight){	
					pred.training <- TPR.DAG(training, g, root=root, positives=positives, bottomup=bottomup, w=w, t=t);
					target.training <- ann[rownames(pred.training),colnames(pred.training)];
					training.avg.meas <- find.best.f(target.training, pred.training, n.round=n.round, f.criterion=f.criterion, 
						verbose=FALSE, b.per.example=FALSE); 
					training.Fmax <- training.avg.meas[4];	## F.max maximization...
					if(training.Fmax > top.Fmax){
						top.Fmax <- training.Fmax;
						best.Fmaxt <- t;
						best.Fmaxw <- w;
						if(bottomup=="threshold" || bottomup=="tau"){
							cat("training fold:",k, "better Fmax.avg found:",top.Fmax, "best threshold:",best.Fmaxt, sep="\t", "\n");
						}else if(bottomup=="weighted.threshold.free"){
							cat("training fold:",k, "better Fmax.avg found:",top.Fmax, "best weight:", best.Fmaxw, sep="\t", "\n");
						}
						else{
							cat("training fold:",k, "better Fmax.avg found:",top.Fmax, "best threshold:",best.Fmaxt, 
								"best weight:", best.Fmaxw, sep="\t", "\n");
						}
					}
				}
			}
		}
		S.test <- TPR.DAG(S.test, g, root=root, positives=positives, bottomup=bottomup, t=best.Fmaxt, w=best.Fmaxw);
					
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
		rm(S, S.test, S.training, training);
	}
	## Storing Results #########
	## Storing Results #########
	if(positives=="children" && bottomup=="threshold.free"){
		meth.name <- "tprTF";
	}else if(positives=="children" && bottomup=="threshold"){
		meth.name <- "tprT";
	}else if(positives=="children" && bottomup=="weighted.threshold.free"){
		meth.name <- "tprW";
	}else if(positives=="children" && bottomup=="weighted.threshold"){
		meth.name <- "tprWT";
	}
	if(positives=="descendants" && bottomup=="threshold.free"){
		meth.name <- "descensTF";
	}else if(positives=="descendants" && bottomup=="threshold"){
		meth.name <- "descensT";
	}else if(positives=="descendants" && bottomup=="weighted.threshold.free"){
		meth.name <- "descensW";
	}else if(positives=="descendants" && bottomup=="weighted.threshold"){
		meth.name <- "descensWT";
	}else if(positives=="descendants" && bottomup=="tau"){
		meth.name <- "descensTAU";
	}
	if(norm){
		save(S.hier, file=paste0(hierScore.dir, flat.file, ".hierScores.",meth.name,".holdout.rda"), compress=TRUE);
		save(AUC.flat, AUC.hier, file=paste0(perf.dir, "AUC.", flat.file, ".hierScores.",meth.name,".holdout.rda"), compress=TRUE);
		save(PXR.flat, PXR.hier, file=paste0(perf.dir, "PXR.", flat.file, ".hierScores.",meth.name,".holdout.rda"), compress=TRUE);
		save(FMM.flat, FMM.hier, file=paste0(perf.dir, "FMM.", flat.file, ".hierScores.",meth.name,".holdout.rda"), compress=TRUE);
		save(PRC.flat, PRC.hier, file=paste0(perf.dir, "PRC.", flat.file, ".hierScores.",meth.name,".holdout.rda"), compress=TRUE);
	}else{
		save(S.hier, file=paste0(hierScore.dir, norm.type,".", flat.file, ".hierScores.",meth.name,".holdout.rda"), compress=TRUE);
		save(AUC.flat, AUC.hier, file=paste0(perf.dir, "AUC.", norm.type, ".", flat.file, ".hierScores.",meth.name,".holdout.rda"), compress=TRUE);
		save(PXR.flat, PXR.hier, file=paste0(perf.dir, "PXR.", norm.type, ".", flat.file, ".hierScores.",meth.name,".holdout.rda"), compress=TRUE);
		save(FMM.flat, FMM.hier, file=paste0(perf.dir, "FMM.", norm.type, ".", flat.file, ".hierScores.",meth.name,".holdout.rda"), compress=TRUE);
		save(PRC.flat, PRC.hier, file=paste0(perf.dir, "PRC.", norm.type, ".", flat.file, ".hierScores.",meth.name,".holdout.rda"), compress=TRUE);
	}
}

