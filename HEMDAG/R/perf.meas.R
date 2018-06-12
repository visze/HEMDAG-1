##***********************##
## TERM-CENTRIC MEASURES ##
##***********************##

#' @title AUPRC single class 
#' @description High-level function to compute the Area under the Precision Recall Curve (AUPRC) just for a single class through \pkg{precrec} package
#' @param target vector of the true labels (0 negative, 1 positive examples)
#' @param pred numeric vector of the values of the predicted labels (scores)
#' @return a numeric value corresponding to the AUPRC for the considered class 
#' @seealso \code{\link{AUPRC.single.over.classes}}
#' @export
#' @examples
#' data(labels);
#' data(scores);
#' data(graph);
#' root <- root.node(g);
#' L <- L[,-which(colnames(L)==root)];
#' S <- S[,-which(colnames(S)==root)];
#' PRC <- AUPRC.single.class(L[,3],S[,3]);
AUPRC.single.class <- function(target, pred){
	if (length(pred)!=length(target))
		stop("AUPRC.single.class: length of true and predicted labels does not match.");
	if (any((target!=0) & (target!=1)))
		stop("AUPRC.single.class: labels variable must take values 0 or 1");
	if ((all(target==0)) || (all(target==1))) ## degenerate case when all labels are equal 
		return (0);
	if(sum(target)==0){
		PRC <- 0;
		return(PRC)
	}
	else{
		res <- evalmod(scores=pred, labels=target);
		aucs <- auc(res);
		prc <- subset(aucs, curvetypes == "PRC");
		PRC <- prc$aucs;
		return(PRC);
	}
}

#' @title AUROC single class 
#' @description High-level function to compute the Area under the ROC Curve (AUPRC) just for a single class through \pkg{precrec} package
#' @param target vector of the true labels (0 negative, 1 positive examples)
#' @param pred numeric vector of the values of the predicted labels (scores)
#' @return a numeric value corresponding to the AUROC for the considered class 
#' @seealso \code{\link{AUROC.single.over.classes}}
#' @export
#' @examples
#' data(labels);
#' data(scores);
#' data(graph);
#' root <- root.node(g);
#' L <- L[,-which(colnames(L)==root)];
#' S <- S[,-which(colnames(S)==root)];
#' AUC <- AUROC.single.class(L[,3],S[,3]);
AUROC.single.class <- function(target, pred){
	if (length(pred)!=length(target))
		stop("AUROC.single.class: length of true and predicted labels does not match.");
	if (any((target!=0) & (target!=1)))
		stop("AUROC.single.class: labels variable must take values 0 or 1");
	if(sum(target)==0){
		AUC <- 0.5;
		return(AUC)
	}
	else{
		res <- evalmod(scores=pred, labels=target);
		aucs <- auc(res);
		roc <- subset(aucs, curvetypes == "ROC");
		AUC <- roc$aucs;
		return(AUC);
	}
}

#' @title AUPRC over classes
#' @description High-level function to compute the Area under the Precision Recall Curve (AUPRC) across a set of classes through \pkg{precrec} package
#' @param target matrix with the target multilabels: rows correspond to examples and columns to classes. 
#' \eqn{target[i,j]=1} if example \eqn{i} belongs to class \eqn{j}, \eqn{target[i,j]=0} otherwise.
#' @param pred a numeric matrix with predicted values (scores): rows correspond to examples and columns to classes.
#' @return a list with two elements:
#' \enumerate{
#' \item average: the average AUPRC across classes;        
#' \item per.class: a named vector with AUPRC for each class. Names correspond to classes
#' }
#' @seealso \code{\link{AUPRC.single.class}}
#' @export
#' @examples
#' data(labels);
#' data(scores);
#' data(graph);
#' root <- root.node(g);
#' L <- L[,-which(colnames(L)==root)];
#' S <- S[,-which(colnames(S)==root)];
#' PRC <- AUPRC.single.over.classes(L,S);
AUPRC.single.over.classes <- function(target, pred){
	n.examples <- nrow(target);
	n.classes <- ncol(target);
	if ((n.examples!=nrow(pred)) || (n.classes!=ncol(pred)))
		stop ("AUPRC.single.over.classes: number of rows or columns do not match between target and predicted classes");
	   
	## if there are classes with zero annotations, we remove them..
	target.names <- colnames(target);
	class.ann <- apply(target,2,sum);
	class.noann <- which(class.ann==0);
	check <- length(class.noann)!=0;
	check.degen <- length(class.noann)!=ncol(target); ## degenerate case when all the classes have zero annotation
	if(check & check.degen){
		target <- target[,-class.noann];
		pred <- pred[,-class.noann];
	}
	## degenerate case when target and pred are vector: just one class have an annotation. May happen in cross validation..
	if(!is.matrix(target)){
		target <- as.matrix(target);
		selected <- which(class.ann==1)
		colnames(target) <- names(selected);
	}
	if(!is.matrix(pred)){
		pred <- as.matrix(pred);
		selected <- which(class.ann==1)
		colnames(pred) <- names(selected);
	}

	## compute PRC considering only those class with non-zero annotations
	PRC.class <- rep(0,ncol(pred));
	names(PRC.class) <- colnames(pred);
	for(i in 1:ncol(pred)){
		PRC.class[i] <- AUPRC.single.class(target[,i],pred[,i]);
	}
	
	## if there are classes with zero annotations, set the prc of those classes to zero and restore the start classes order 
	if(check & check.degen){
		PRC.class <- PRC.class[target.names];
		PRC.class[is.na(PRC.class)] <- 0;
		names(PRC.class) <- target.names; 
	}

	#saving PRC result in the same format of package PerfMeas
	PRC.mean <- mean(PRC.class);
	PRC.res <- list(average=PRC.mean, per.class=PRC.class); 
	return(PRC.res);
}

#' @title AUROC over classes
#' @description High-level function to compute the Area under the ROC Curve (AUROC) for a set of classes through \pkg{precrec} package
#' @param target matrix with the target multilabels: rows correspond to examples and columns to classes. 
#' \eqn{target[i,j]=1} if example \eqn{i} belongs to class \eqn{j}, \eqn{target[i,j]=0} otherwise.
#' @param pred a numeric matrix with predicted values (scores): rows correspond to examples and columns to classes.
#' @return a list with two elements:
#' \enumerate{
#' \item average: the average AUROC across classes;        
#' \item per.class: a named vector with AUROC for each class. Names correspond to classes
#' }
#' @seealso \code{\link{AUROC.single.class}}
#' @export
#' @examples
#' data(labels);
#' data(scores);
#' data(graph);
#' root <- root.node(g);
#' L <- L[,-which(colnames(L)==root)];
#' S <- S[,-which(colnames(S)==root)];
#' AUC <- AUROC.single.over.classes(L,S);
AUROC.single.over.classes <- function(target, pred){
	n.examples <- nrow(target);
	n.classes <- ncol(target);
	if ((n.examples!=nrow(pred)) || (n.classes!=ncol(pred)))
		stop ("AUROC.single.over.classes: number of rows or columns do not match between target and predicted classes");
	
	## if there are classes with zero annotations, we remove them..
	target.names <- colnames(target);
	class.ann <- apply(target,2,sum);
	class.noann <- which(class.ann==0);
	check <- length(class.noann)!=0;
	check.degen <- length(class.noann)!=ncol(target); ## degenerate case when all the classes have zero annotation
	if(check & check.degen){
		target <- target[,-class.noann];
		pred <- pred[,-class.noann];
	}
	## degenerate case when target and pred are vector: just one class have an annotation. May happen in cross validation..
	if(!is.matrix(target)){
		target <- as.matrix(target);
		selected <- which(class.ann==1)
		colnames(target) <- names(selected);
	}
	if(!is.matrix(pred)){
		pred <- as.matrix(pred);
		selected <- which(class.ann==1)
		colnames(pred) <- names(selected);
	}
	
	## compute AUC considering only those class with non-zero annotations
	AUC.class <- rep(0,ncol(pred));
	names(AUC.class) <- colnames(pred);
	for(i in 1:ncol(pred)){
		AUC.class[i] <- AUROC.single.class(target[,i],pred[,i]); 
	}

	## if there are classes with zero annotations, set the AUC of those classes to zero and restore the start classes order 
	if(check & check.degen){
		AUC.class <- AUC.class[target.names];
		AUC.class[is.na(AUC.class)] <- 0.5;
		names(AUC.class) <- target.names; 
	}
		
	#saving AUC result in the same format of package PerfMeas
	AUC.mean <- mean(AUC.class);
	AUC.res <- list(average=AUC.mean, per.class=AUC.class); 
	return(AUC.res);
}

##************************************************************##
## Functions to compute Kiritchenko-like multi-label F-scores ##
##************************************************************##

#' @name Multilabel.F.measure
#' @aliases F.measure.multilabel
#' @title Multilabel F-measure 
#' @description Method for computing Precision, Recall, Specificity, Accuracy and F-measure for multiclass multilabel classification
#' @param target matrix with the target multilabels: rows correspond to examples and columns to classes.
#' \eqn{target[i,j]=1} if example \eqn{i} belongs to class \eqn{j}, \eqn{target[i,j]=0} otherwise
#' @param predicted a numeric matrix with predicted values (scores): rows correspond to examples and columns to classes
#' @param b.per.example boolean. 
#' \itemize{
#'	\item \code{TRUE}: results are returned for each example;
#'	\item \code{FALSE}: only the average results are returned
#' }
#' @return Two different outputs respect to the input paramenter \code{b.per.example}:
#' \itemize{
#'	\item \code{b.per.example==FALSE}: a list with a single element average. A named vector with average precision (P), recall (R), 
#' 	specificity (S), F-measure (F), average F-measure (avF) and Accuracy (A) across examples. F is the F-measure computed as the 
#' 	harmonic mean between the average precision and recall; av.F is the F-measure computed as the average across examples.
#' 	\item \code{b.per.example==FALSE}: a list with two elements:
#' 		\enumerate{
#' 			\item average: a named vector with average precision (P), recall (R), specificity (S), F-measure (F), average F-measure (avF) 
#' 			and Accuracy (A) across examples; 
#' 			\item per.example: a named matrix with the Precision (P), Recall (R), Specificity (S), Accuracy (A), F-measure (F) and 
#' 			av.F-measure (av.F) for each example. Row names correspond to examples, column names correspond respectively to Precision (P), Recall (R), 
#' 			Specificity (S), Accuracy (A), F-measure (F) and av.F-measure (av.F)
#' 		}
#' }
#' @examples
#' data(labels);
#' data(scores);
#' data(graph);
#' root <- root.node(g);
#' L <- L[,-which(colnames(L)==root)];
#' S <- S[,-which(colnames(S)==root)];
#' S[S>0.7] <- 1;
#' S[S<0.7] <- 0;
#' FMM <- F.measure.multilabel(L,S);
#' @export
#' @docType methods
setGeneric("F.measure.multilabel", 
	function(target, predicted, b.per.example=FALSE) standardGeneric("F.measure.multilabel"));

#' @rdname Multilabel.F.measure
setMethod("F.measure.multilabel", signature(target="matrix", predicted="matrix"),
	function(target, predicted, b.per.example=FALSE) { 
		n.examples <- nrow(target);
		n.classes <- ncol(target);
		if ((n.examples!=nrow(predicted)) || (n.classes!=ncol(predicted)))
			stop ("F.measure.multilabel: number of rows or columns do not match between target and predicted classes");
	  
		z <- target + predicted;
		TP <- apply(z, 1, function(x){
			return(sum(x==2));
			});
		TN <- apply(z, 1, function(x){
			return(sum(x==0));
		});
		z <- predicted - target;
		FP <- apply(z, 1, function(x){
			return(sum(x==1));
		});
		FN <- apply(z, 1, function(x){
			return(sum(x== -1));
		});
		rm(z);
		n <- sum(TP)+sum(TN)+sum(FN)+sum(FP);
		if( n != (n.examples*n.classes)){ 
			cat("n = ", n, "\n n.examples = ", n.examples, "\n n.classes = ", n.classes, "\n");
			cat (" sum(TP) = ", sum(TP), "\n sum(TN) = ", sum(TN), "\n sum(FN) = ", sum(FN), "\n sum(FP) = ", sum(FP), "\n");
			warning("F.measure.multilabel: Something went wrong in F-measure)");
		}
	   
		P <- TP+FP;
		P[which(P==0)] <- 1;  # to avoid division by 0 in precision
	   
		sum.TP.FN <- TP+FN;
		sum.TN.FP <- TN+FP;
	   
		sum.TP.FN[which(sum.TP.FN==0)] <- 1;  # to avoid division by 0 in recall
		sum.TN.FP[which(sum.TN.FP==0)] <- 1;  # to avoid division by 0 in specificity
		  
		precision <- TP/P;
		recall <- TP/sum.TP.FN;
		specificity <- TN/sum.TN.FP;
	   
		prec.rec <- precision+recall;
		prec.rec[which(prec.rec==0)] <- 1;  # to avoid division by 0 for f.measure
		f.measure <- (2*precision*recall)/prec.rec;
		accuracy <- (TP+TN)/n.classes;
	   
		av.precision <- sum(precision)/n.examples; 
		av.recall <- sum(recall)/n.examples; 
		av.specificity <- sum(specificity)/n.examples; 
		av.prec.rec <- av.precision+av.recall;
		if(av.prec.rec == 0){
			av.prec.rec <- 1;
		}
		overall.av.f.measure <- (2*av.precision*av.recall)/av.prec.rec;
		av.f.measure <- sum(f.measure)/n.examples; 
		av.accuracy  <- sum(accuracy)/n.examples; 
	   
		average <- c(av.precision, av.recall, av.specificity, overall.av.f.measure, av.f.measure,av.accuracy);
		names(average) <- c("P", "R", "S", "F", "avF", "A");
	   
	   if(b.per.example){
			per.example <- cbind(precision, recall, specificity, f.measure, accuracy);
			colnames(per.example) <- c("P", "R", "S", "F","A");
			return (list(average=average, per.example=per.example))
		}else{
			return (list(average=average));
		}
	} 
)

#' @title Best hierarchical F-score 
#' @description Function to select the best hierarchical F-score by choosing an appropriate threshold in the scores
#' @details All the examples having no positive annotations are discarded. The predicted scores matrix (\code{pred}) is rounded 
#' according to parameter \code{n.round} and all the values of \code{pred} are divided by \code{max(pred)}.
#' Then all the thresholds corresponding to all the different values included in \code{pred} are attempted, and the threshold 
#' leading to the maximum F-measure is selected.
#' @param target matrix with the target multilabels: rows correspond to examples and columns to classes.
#' \eqn{target[i,j]=1} if example \eqn{i} belongs to class \eqn{j}, \eqn{target[i,j]=0} otherwise
#' @param pred a numeric matrix with predicted values (scores): rows correspond to examples and columns to classes
#' @param n.round number of rounding digits to be applied to pred (\code{default=3})
#' @param f.criterion character. Type of F-measure to be used to select the best F-score. There are two possibilities:
#' \enumerate{
#'	\item \code{F} (def.) corresponds to the harmonic mean between the average precision and recall;
#'	\item \code{avF} corresponds to the per-example F-score averaged across all the examples.
#' }
#' @param verbose boolean. If \code{TRUE} (def.) the number of iterations are printed on stdout
#' @param b.per.example boolean. 
#' \itemize{
#'	\item \code{TRUE}: results are returned for each example;
#'	\item \code{FALSE}: only the average results are returned
#' }
#' @return Two different outputs respect to the input paramenter \code{b.per.example}:
#' \itemize{
#'	\item \code{b.per.example==FALSE}: a list with a single element average. A named vector with 7 elements relative to the best result in terms 
#' 	of the F.measure: Precision (P), Recall (R), Specificity (S), F.measure (F), av.F.measure (av.F), Accuracy (A) and the best selected Threshold (T). 
#' 	F is the F-measure computed as the harmonic mean between the average precision and recall; av.F is the F-measure computed as the average across 
#' examples and T is the best selected threshold;
#' 	\item \code{b.per.example==FALSE}: a list with two elements:
#' 		\enumerate{
#' 			\item average: a named vector with with 7 elements relative to the best result in terms of the F.measure: Precision (P), Recall (R), 
#'			Specificity (S), F.measure (F), av.F.measure (av.F), Accuracy (A) and the best selected Threshold (T). 
#' 			\item per.example: a named matrix with the Precision (P), Recall (R), Specificity (S), Accuracy (A), F-measure (F),	av.F-measure (av.F)
#' 			and the best selected Threhold (T) for each example. Row names correspond to examples, column names correspond respectively 
#'			to Precision (P), Recall (R), Specificity (S), Accuracy (A), F-measure (F), av.F-measure (av.F) and the best selected Threhold (T).
#' 		}
#' }
#' @export
#' @examples
#' data(graph);
#' data(labels);
#' data(scores);
#' root <- root.node(g);
#' L <- L[,-which(colnames(L)==root)];
#' S <- S[,-which(colnames(S)==root)];
#' FMM <- find.best.f(L,S,n.round=3, f.criterion ="F", verbose=TRUE, b.per.example=TRUE);
find.best.f <- function(target, pred, n.round=3, f.criterion ="F", verbose=TRUE, b.per.example=FALSE){
	x <- apply(target,1,sum);
	selected <- which(x>0);
	##  degenerate case when target is a full-zero matrix (all genes without annotations)
	if(length(selected)==0){
		selected <- which(x==0);
	}
	target <- target[selected,];
	## degenerate case when target is a vector (just one annotated gene)
	if(!is.matrix(target)){
		target <- t(as.matrix(target));
		rownames(target) <- names(selected);
	}
	pred <- pred[selected,];
	## degenerate case when pred is a vector (just one annotated gene)
	if(!is.matrix(pred)){
		pred <- t(as.matrix(pred));
		rownames(pred) <- names(selected);
	}
	pred <- pred/max(pred);
	pred <- round(pred,n.round);
	n.examples <- nrow(pred);
	n.classes <- ncol(pred);

	thresh <- unique(as.numeric(pred));
	thresh <- sort(thresh);
	best.res <- best <- best.thresh <- 0;
	i=0;
	for (t in thresh){
		pred.labels <- matrix(numeric(n.examples*n.classes), nrow=n.examples);
		pred.labels[pred>=t] <-1;
		res <- F.measure.multilabel(target, pred.labels, b.per.example);
		if(res$average[f.criterion] > best){
			best <- res$average[f.criterion];
			best.res <- res;  
			best.thresh <- t;
		}
		i <- i+1;
		if (i%%100 == 0  && verbose){
			cat("iteration ", i,  "\n");
		}
	}
	## degenerate case when target is a full-zero matrix: by.def F-score is zero
	if(!is.list(best.res)){
		best.res <- res;
	}
	if(b.per.example){
		best.res$average <- c(best.res$average, best.thresh);
		names(best.res$average)[7] <- "T"; 
		return(best.res);
	}else{
		best.res <- c(best.res$average, best.thresh);
		names(best.res)[7] <- "T";
		return(best.res);
	}
}
