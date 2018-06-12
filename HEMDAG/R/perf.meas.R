##***************##
## AUORC & AUROC ##
##***************##
#' @name AUPRC 
#' @aliases AUPRC.single.class
#' @aliases AUPRC.single.over.classes
#' @title AUPRC measures
#' @description Function to compute Area under the Precision Recall Curve (AUPRC) through \pkg{precrec} package
#' @details The AUPRC (for a single class or for a set of classes) is computed either one-shot or averaged across stratified folds.\cr
#' \code{AUPRC.single.class} computes the AUPRC just for a given class.\cr
#' \code{AUPR.single.over.classes} computes the AUPRC for a set of classes, including their average values across all the classes.\cr
#' For all those classes having zero annotations, the AUPRC is set to 0. These classes are discarded in the computing of the AUPRC   
#' averaged across classes, both when the AUPRC is computed one-shot or averaged across stratified folds.
#' @param folds number of folds on which computing the AUPRC. If \code{folds=NULL} (\code{def.}), the one-shot AUPRC is computed, 
#' otherwise the AUPRC is computed averaged across folds.
#' @param seed initialization seed for the random generator to create folds. Set \code{seed} only if \code{folds}\eqn{\neq}\code{NULL}.
#' If \code{seed}\eqn{\neq}\code{NULL} and \code{folds=NULL}, the AUPRC averaged across folds is computed without seed initialization.
#' @param target vector of the true labels (0 negative, 1 positive examples)
#' @param labels matrix with the target multilabels: rows correspond to examples and columns to classes. 
#' \eqn{target[i,j]=1} if example \eqn{i} belongs to class \eqn{j}, \eqn{target[i,j]=0} otherwise.
#' @param pred a numeric vector of the values of the predicted labels (scores)
#' @param predicted a numeric matrix with predicted values (scores): rows correspond to examples and columns to classes.
#' @return \code{AUPRC.single.class} returns a numeric value corresponding to the AUPRC for the considered class;\cr
#' \code{AUPR.single.over.classes} returns a list with two elements:
#' \enumerate{
#' 	\item average: the average AUPRC across classes;     
#' 	\item per.class: a named vector with AUPRC for each class. Names correspond to classes.
#' }
#' @export
#' @examples
#' data(labels);
#' data(scores);
#' data(graph);
#' root <- root.node(g);
#' L <- L[,-which(colnames(L)==root)];
#' S <- S[,-which(colnames(S)==root)];
#' PRC.single.class <- AUPRC.single.class(L[,3], S[,3], folds=5, seed=23);
#' PRC.over.classes <- AUPRC.single.over.classes(L, S, folds=5, seed=23);
AUPRC.single.class <- function(target, pred, folds=NULL, seed=NULL){
	if(is.matrix(target) || is.matrix(pred))
		stop("target or pred must be a vector", call.=FALSE);
	if(length(pred)!=length(target))
		stop("AUPRC.single.class: length of true and predicted labels does not match", call.=FALSE);
	if(any((target!=0) & (target!=1)))
		stop("AUPRC.single.class: labels variable must take values 0 or 1", call.=FALSE);
	if(is.null(folds) && !is.null(seed))
		seed <- NULL;
	if(!is.null(folds) && is.null(seed))
		warning("AUPRC.single.class: folds are generated without seed initialization", call.=FALSE);
	## degenerate case when all labels are equals
	if((all(target==0)) || (all(target==1))){
		if(!is.null(folds)){
			PRC.avg <- 0;
			PRC.fold <- rep(0,folds);
			PRC <- list(average=PRC.avg, across.fold=PRC.fold);
		}else{
			PRC <- 0;
			names(PRC) <- "one.shoot";
		}
		return(PRC);
	}
	if(is.null(folds) && !is.null(seed))
		seed <- NULL;
	## compute PRC averaged across folds
	if(!is.null(folds)){
		indices <- 1:length(target);
		positives <- which(target==1);
		foldIndex <- do.stratified.cv.data.single.class(indices, positives, kk=folds, seed=seed);
		testIndex <- mapply(c, foldIndex$fold.positives, foldIndex$fold.negatives, SIMPLIFY=FALSE); ## index of examples used for test set..
		fold.check <- unlist(lapply(testIndex,length));
		if(any(fold.check==0))
			stop("AUPRC.single.class: Number of folds selected too high: some folds have no examples. Please reduce the number of folds", call.=FALSE);
		PRC.fold <- numeric(folds);
		for(k in 1:folds){
			target.test <- target[testIndex[[k]]];
			pred.test <- pred[testIndex[[k]]];
			if(sum(target.test) > 0){
				if(all(target.test==0) || all(target.test==1)){ ## degenerate case when all labels in the k fold are equals 
					PRC.fold[k] <- 0;
				}else{
					res <- evalmod(scores=pred.test, labels=target.test);
					aucs <- auc(res);
					prc <- subset(aucs, curvetypes=="PRC");
					PRC.fold[k] <- prc$aucs;
				}
			}else{
				PRC.fold[k] <- 0;
			}
		}
		if(all(PRC.fold==0)){
			PRC.avg <- 0; ## degenerate case in which for each fold all the PRC computed by precrec are 0 
		}else{
			classIndex <- which(PRC.fold!=0);
			classPerfs <- PRC.fold[classIndex];
			PRC.avg <- mean(classPerfs);
		}
		PRC <- list(average=PRC.avg, across.fold=PRC.fold);
		return(PRC);
	}
	## compute PRC one-shoot
	res <- evalmod(scores=pred, labels=target);
	aucs <- auc(res);
	prc <- subset(aucs, curvetypes=="PRC");
	PRC <- prc$aucs;
	names(PRC) <- "one.shoot";
	return(PRC);
}

#' @rdname AUPRC
#' @export 
AUPRC.single.over.classes <- function(labels, predicted, folds=NULL, seed=NULL){
	if(!is.matrix(labels) || !is.matrix(predicted))
		stop("AUPRC.single.over.classes: labels or predicted must be a matrix", call.=FALSE);
	n.examples <- nrow(labels);
	n.classes <- ncol(labels);
	if((n.examples!=nrow(predicted)) || (n.classes!=ncol(predicted)))
		stop ("AUPRC.single.over.classes: number of rows or columns do not match between labels and predicted classes", call.=FALSE);
	if(any((labels!=0) & (labels!=1)))
		stop("AUPRC.single.over.classes: labels variable must take values 0 or 1", call.=FALSE);
	if(is.null(folds) && !is.null(seed))
		seed <- NULL;
	## AUPRC averaged across folds over classes
	if(!is.null(folds)){
		PRC.class <- rep(0,ncol(predicted));
		names(PRC.class) <- colnames(predicted);
		for(i in 1:ncol(predicted)){
			PRC.class[i] <- AUPRC.single.class(labels[,i], predicted[,i], folds=folds, seed=seed)$average;
		}
		if(all(PRC.class==0)){
			PRC.avg <- 0;
		}else{
			classIndex <- which(PRC.class!=0);
			classPerfs <- PRC.class[classIndex];
			PRC.avg <- mean(classPerfs);
		}
		PRC.res <- list(average=PRC.avg, per.class=PRC.class); 
		return(PRC.res);
	}
	## AUPRC one-shot over classes
	## if there are classes with zero annotations, we remove them..
	labels.names <- colnames(labels);
	class.ann <- apply(labels,2,sum);
	class.noann <- which(class.ann==0);
	check <- length(class.noann)!=0;
	check.degen <- length(class.noann)!=ncol(labels); ## degenerate case when all the classes have zero annotations
	if(check & check.degen){
		labels <- labels[,-class.noann];
		predicted <- predicted[,-class.noann];
	}
	## degenerate case when labels and predicted are vector: just one class has an annotation. Might happen in cross validation..
	if(!is.matrix(labels)){
		labels <- as.matrix(labels);
		selected <- which(class.ann==1);
		colnames(labels) <- names(selected);
	}
	if(!is.matrix(predicted)){
		predicted <- as.matrix(predicted);
		selected <- which(class.ann==1)
		colnames(predicted) <- names(selected);
	}
	## compute PRC considering only those class with non-zero annotations
	PRC.class <- rep(0,ncol(predicted));
	names(PRC.class) <- colnames(predicted);
	for(i in 1:ncol(predicted)){
		PRC.class[i] <- AUPRC.single.class(labels[,i], predicted[,i], folds=NULL, seed=NULL);
	}	
	## if there are classes with zero annotations, set the prc of those classes to zero and restore the start classes order 
	if(check & check.degen){
		PRC.class <- PRC.class[labels.names];
		PRC.class[is.na(PRC.class)] <- 0;
		names(PRC.class) <- labels.names; 
	}
	if(all(PRC.class==0)){
			PRC.avg <- 0;
	}else{
		classIndex <- which(PRC.class!=0);
		classPerfs <- PRC.class[classIndex];
		PRC.avg <- mean(classPerfs);
	}
	PRC.res <- list(average=PRC.avg, per.class=PRC.class); 
	return(PRC.res);
}

#' @name AUROC 
#' @aliases AUROC.single.class
#' @aliases AUROC.single.over.classes
#' @title AUROC measures
#' @description Function to compute the Area under the ROC Curve through \pkg{precrec} package
#' @details The AUROC (for a single class or for a set of classes) is computed either one-shot or averaged across stratified folds.\cr
#' \code{AUROC.single.class} computes the AUROC just for a given class.\cr
#' \code{AUROC.single.over.classes} computes the AUROC for a set of classes, including their average values across all the classes.\cr
#' For all those classes having zero annotations, the AUROC is set to 0.5. These classes are included in the computing of the AUROC 
#' averaged across classes, both when the AUROC is computed one-shot or averaged across stratified folds.
#' The AUROC is set to 0.5 to all those classes having zero annotations. 
#' @param folds number of folds on which computing the AUROC. If \code{folds=NULL} (\code{def.}), the one-shot AUROC is computed, 
#' otherwise the AUROC is computed averaged across folds.
#' @param seed initialization seed for the random generator to create folds. Set \code{seed} only if \code{folds}\eqn{\neq}\code{NULL}.
#' If \code{seed}\eqn{\neq}\code{NULL} and \code{folds=NULL}, the AUROC averaged across folds is computed without seed initialization.
#' @param target vector of the true labels (0 negative, 1 positive examples)
#' @param labels matrix with the target multilabels: rows correspond to examples and columns to classes. 
#' \eqn{target[i,j]=1} if example \eqn{i} belongs to class \eqn{j}, \eqn{target[i,j]=0} otherwise.
#' @param pred a numeric vector of the values of the predicted labels (scores)
#' @param predicted a numeric matrix with predicted values (scores): rows correspond to examples and columns to classes.
#' @return \code{AUROC.single.class} returns a numeric value corresponding to the AUROC for the considered class;\cr
#' \code{AUPR.single.over.classes} returns a list with two elements:
#' \enumerate{
#' 	\item average: the average AUROC across classes;        
#' 	\item per.class: a named vector with AUROC for each class. Names correspond to classes.
#' }
#' @export
#' @examples
#' data(labels);
#' data(scores);
#' data(graph);
#' root <- root.node(g);
#' L <- L[,-which(colnames(L)==root)];
#' S <- S[,-which(colnames(S)==root)];
#' AUC.single.class <- AUROC.single.class(L[,3], S[,3], folds=5, seed=23);
#' AUC.over.classes <- AUROC.single.over.classes(L, S, folds=5, seed=23);
AUROC.single.class <- function(target, pred, folds=NULL, seed=NULL){
	if(is.matrix(target) || is.matrix(pred))
		stop("AUROC.single.class: target or pred must be a vector", call.=FALSE);
	if(length(pred)!=length(target))
		stop("AUROC.single.class: length of true and predicted labels does not match.", call.=FALSE);
	if(any((target!=0) & (target!=1)))
		stop("AUROC.single.class: labels variable must take values 0 or 1", call.=FALSE);
	if(is.null(folds) && !is.null(seed))
		seed <- NULL;
	if(!is.null(folds) && is.null(seed))
		warning("AUROC.single.class: folds are generated without seed initialization", call.=FALSE);
	## degenerate case when all labels are equals
	if((all(target==0)) || (all(target==1))){
		if(!is.null(folds)){
			AUC.avg <- 0.5;
			AUC.fold <- rep(0.5,folds);
			AUC <- list(average=AUC.avg, across.fold=AUC.fold);
		}else{
			AUC <- 0.5;
			names(AUC) <- "one.shoot";
		}
		return(AUC);
	}
	if(is.null(folds) && !is.null(seed))
		seed <- NULL;
	## compute AUC averaged across folds
	if(!is.null(folds)){
		indices <- 1:length(target);
		positives <- which(target==1);
		foldIndex <- do.stratified.cv.data.single.class(indices, positives, kk=folds, seed=seed);
		testIndex <- mapply(c, foldIndex$fold.positives, foldIndex$fold.negatives, SIMPLIFY=FALSE); ## index of examples used for test set..
		fold.check <- unlist(lapply(testIndex,length));
		if(any(fold.check==0))
			stop("AUROC.single.class: Number of folds selected too high: some folds have no examples. Please reduce the number of folds", call.=FALSE);
		AUC.fold <- numeric(folds);
		for(k in 1:folds){
			target.test <- target[testIndex[[k]]];
			pred.test <- pred[testIndex[[k]]];
			if(sum(target.test) > 0){
				if(all(target.test==0) || all(target.test==1)){ ## degenerate case when all labels in the k fold are equals 
					AUC.fold[k] <- 0.5;
				}else{
					res <- evalmod(scores=pred.test, labels=target.test);
					aucs <- auc(res);
					auc <- subset(aucs, curvetypes=="ROC");
					AUC.fold[k] <- auc$aucs;
				}
			}else{				
				AUC.fold[k] <- 0.5;
			}
		}
		AUC.avg <- mean(AUC.fold);
		AUC <- list(average=AUC.avg, across.fold=AUC.fold);
		return(AUC);
	}
	## compute AUC one-shoot
	res <- evalmod(scores=pred, labels=target);
	aucs <- auc(res);
	prc <- subset(aucs, curvetypes=="ROC");
	AUC <- prc$aucs;
	names(AUC) <- "one.shoot";
	return(AUC);
}

#' @rdname AUROC
#' @export 
AUROC.single.over.classes <- function(labels, predicted, folds=NULL, seed=NULL){
	if(!is.matrix(labels) || !is.matrix(predicted))
		stop("AUROC.single.over.classes: labels or predicted must be a matrix", call.=FALSE);
	n.examples <- nrow(labels);
	n.classes <- ncol(labels);
	if((n.examples!=nrow(predicted)) || (n.classes!=ncol(predicted)))
		stop ("AUROC.single.over.classes: number of rows or columns do not match between labels and predicted classes", call.=FALSE);
	if(any((labels!=0) & (labels!=1)))
		stop("AUROC.single.over.classes: labels variable must take values 0 or 1", call.=FALSE);
	if(is.null(folds) && !is.null(seed))
		seed <- NULL;
	## AUROC averaged across folds over classes	
	if(!is.null(folds)){
		AUC.class <- rep(0,ncol(predicted));
		names(AUC.class) <- colnames(predicted);
		for(i in 1:ncol(predicted)){
			AUC.class[i] <- AUROC.single.class(labels[,i], predicted[,i], folds=folds, seed=seed)$average; 
		}
		AUC.avg <- mean(AUC.class);
		AUC.res <- list(average=AUC.avg, per.class=AUC.class); 
		return(AUC.res);
	}
	## if there are classes with zero annotations, we remove them..
	labels.names <- colnames(labels);
	class.ann <- apply(labels,2,sum);
	class.noann <- which(class.ann==0);
	check <- length(class.noann)!=0;
	check.degen <- length(class.noann)!=ncol(labels); ## degenerate case when all the classes have zero annotation
	if(check & check.degen){
		labels <- labels[,-class.noann];
		predicted <- predicted[,-class.noann];
	}
	## degenerate case when labels and predicted are vector: just one class have an annotation. May happen in cross validation..
	if(!is.matrix(labels)){
		labels <- as.matrix(labels);
		selected <- which(class.ann==1)
		colnames(labels) <- names(selected);
	}
	if(!is.matrix(predicted)){
		predicted <- as.matrix(predicted);
		selected <- which(class.ann==1)
		colnames(predicted) <- names(selected);
	}
	## compute AUC considering only those class with non-zero annotations
	AUC.class <- rep(0,ncol(predicted));
	names(AUC.class) <- colnames(predicted);
	for(i in 1:ncol(predicted)){
		AUC.class[i] <- AUROC.single.class(labels[,i],predicted[,i], folds=NULL, seed=NULL); 
	}
	## if there are classes with zero annotations, set the AUC of those classes to zero and restore the start classes order 
	if(check & check.degen){
		AUC.class <- AUC.class[labels.names];
		AUC.class[is.na(AUC.class)] <- 0.5;
		names(AUC.class) <- labels.names; 
	}
	## saving AUC result in the same format of package PerfMeas
	AUC.avg <- mean(AUC.class);
	AUC.res <- list(average=AUC.avg, per.class=AUC.class); 
	return(AUC.res);
}

##************************************************************##
## Functions to compute Kiritchenko-like multilabel F-scores ##
##************************************************************##
#' @name Multilabel.F.measure
#' @aliases F.measure.multilabel
#' @title Multilabel F-measure 
#' @description Method for computing Precision, Recall, Specificity, Accuracy and F-measure for multiclass and multilabel classification
#' @param target matrix with the target multilabels: rows correspond to examples and columns to classes.
#' \eqn{target[i,j]=1} if example \eqn{i} belongs to class \eqn{j}, \eqn{target[i,j]=0} otherwise
#' @param predicted a numeric matrix with discrete predicted values: rows correspond to examples and columns to classes.
#' \eqn{predicted[i,j]=1} if example \eqn{i} is predicted belonging to class \eqn{j}, \eqn{target[i,j]=0} otherwise
#' @param b.per.example boolean. 
#' \itemize{
#'	\item \code{TRUE}: results are returned for each example;
#'	\item \code{FALSE}: only the average results are returned
#' }
#' @return Two different outputs respect to the input parameter \code{b.per.example}:
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
	function(target, predicted, b.per.example=FALSE){
		n.examples <- nrow(target);
		n.classes <- ncol(target);
		if((n.examples!=nrow(predicted)) || (n.classes!=ncol(predicted)))
			stop ("F.measure.multilabel: number of rows or columns do not match between target and predicted classes", call.=FALSE);
		if(any((target!=0) & (target!=1)) || any((predicted!=0) & (predicted!=1)))
			stop("F.measure.multilabel: target and predicted variables must take values 0 or 1", call.=FALSE);
				
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
			warning("F.measure.multilabel: Something went wrong in F-measure", call.=FALSE);
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
#' @param pred a numeric matrix with continuous predicted values (scores): rows correspond to examples and columns to classes
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
#' @return Two different outputs respect to the input parameter \code{b.per.example}:
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
#' 			and the best selected Threshold (T) for each example. Row names correspond to examples, column names correspond respectively 
#'			to Precision (P), Recall (R), Specificity (S), Accuracy (A), F-measure (F), av.F-measure (av.F) and the best selected Threshold (T).
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
#' FMM <- find.best.f(L, S, n.round=3, f.criterion="F", verbose=TRUE, b.per.example=TRUE);
find.best.f <- function(target, pred, n.round=3, f.criterion="F", verbose=TRUE, b.per.example=FALSE){
	if(!is.matrix(target) || !is.matrix(pred))
		stop("find.best.f: target or pred must be a matrix", call.=FALSE);
	n.examples <- nrow(target);
	n.classes <- ncol(target);
	if((n.examples!=nrow(pred)) || (n.classes!=ncol(pred)))
		stop("find.best.f: number of rows or columns do not match between target and predicted classes", call.=FALSE);
	if(any((target!=0) & (target!=1)))
		stop("find.best.f: labels variable must take values 0 or 1", call.=FALSE);
	
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
	i <- 0;
	for(t in thresh){
		pred.labels <- matrix(numeric(n.examples*n.classes), nrow=n.examples);
		pred.labels[pred>=t] <-1;
		res <- F.measure.multilabel(target, pred.labels, b.per.example);
		if(res$average[f.criterion] > best){
			best <- res$average[f.criterion];
			best.res <- res;  
			best.thresh <- t;
		}
		i <- i+1;
		if(i%%100 == 0  && verbose){
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

#' @name FMM
#' @title Compute Multilabel F-measure
#' @description Function to compute the best hierarchical F-score either one-shot or averaged across folds
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
#' @param folds number of folds on which computing the AUROC. If \code{folds=NULL} (\code{def.}), the one-shot AUROC is computed, 
#' otherwise the AUROC is computed averaged across folds.
#' @param seed initialization seed for the random generator to create folds. Set \code{seed} only if \code{folds}\eqn{\neq}\code{NULL}.
#' If \code{seed}\eqn{\neq}\code{NULL} and \code{folds=NULL}, the AUROC averaged across folds is computed without seed initialization.
#' @return Two different outputs respect to the input parameter \code{b.per.example}:
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
#' 			and the best selected Threshold (T) for each example. Row names correspond to examples, column names correspond respectively 
#'			to Precision (P), Recall (R), Specificity (S), Accuracy (A), F-measure (F), av.F-measure (av.F) and the best selected Threshold (T).
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
#' FMM <- compute.Fmeasure.multilabel(L, S, n.round=3, f.criterion="F", verbose=TRUE, 
#' b.per.example=TRUE, folds=5, seed=23);
compute.Fmeasure.multilabel <- function(target, pred, n.round=3, f.criterion="F", verbose=TRUE, b.per.example=FALSE, folds=NULL, seed=NULL){
	if(!is.matrix(target) || !is.matrix(pred))
		stop("compute.Fmeasure.multilabel: target or pred must be a matrix", call.=FALSE);
	n.examples <- nrow(target);
	n.classes <- ncol(target);
	if((n.examples!=nrow(pred)) || (n.classes!=ncol(pred)))
		stop("compute.Fmeasure.multilabel: number of rows or columns do not match between target and predicted classes", call.=FALSE);
	if(any((target!=0) & (target!=1)))
		stop("compute.Fmeasure.multilabel: labels variable must take values 0 or 1", call.=FALSE);
	if(is.null(folds) && !is.null(seed))
		seed <- NULL;
	if(!is.null(folds) && is.null(seed))
		warning("compute.Fmeasure.multilabel: folds are generated without seed initialization", call.=FALSE);	
	## FMM averaged across folds 
	if(!is.null(folds)){
		testIndex <- do.unstratified.cv.data(pred, kk=folds, seed=seed);
		avg.res.list <- list();
		res.per.example <- c();
		for(k in 1:folds){
			fold.res <- find.best.f(target[testIndex[[k]],], pred[testIndex[[k]],], n.round=n.round, f.criterion=f.criterion, 
				verbose=verbose, b.per.example=b.per.example);
			avg.res.list[[k]] <- fold.res$average;
			res.per.example <- rbind(res.per.example, fold.res$per.example);
		}
		F.cv <- Reduce("+", avg.res.list)/folds;
		F.meas <-  apply(res.per.example,2,mean);
		names(F.meas)[4] <- "avF";
		## degenerate case when both precision and recall are zero..
		if((F.meas[["P"]] && F.meas[["R"]]) == 0){ 	## sum(res.per.example[,"F"])==0
			F.max <- 0;
		}else{
			F.max <- 2*(F.meas[["P"]] * F.meas[["R"]])/((F.meas[["P"]] + F.meas[["R"]]));
		}
		names(F.max) <- "F";
		FMM.avg <- append(F.meas, F.max, after=3);
		FMM.avg <- append(FMM.avg, F.cv["T"]);
		res <- list(average=FMM.avg, per.example=res.per.example);
		if(b.per.example){
			return(res);
		}else{
			res <- res$average;
			return(res);
		}		
	}
	## FMM one-shot
	res <- find.best.f(target=target, pred=pred, n.round=n.round, f.criterion=f.criterion, 
		verbose=verbose, b.per.example=b.per.example);
	return(res);
}

#' @name PXR
#' @title Precision at fixed Recall level
#' @description Function to compute the Precision at fixed Recall levels (PXR) over classes through \pkg{PerfMeas} package
#' @details The PXR across all classes is computed either one-shot or averaged across stratified folds
#' @details Function to compute the precision at fixed recall levels across all classes either one-shot or averaged across stratified folds
#' @param target matrix with the target multilabels: rows correspond to examples and columns to classes. 
#' \eqn{target[i,j]=1} if example \eqn{i} belongs to class \eqn{j}, \eqn{target[i,j]=0} otherwise.
#' @param pred a numeric matrix with predicted values (scores): rows correspond to examples and columns to classes.
#' @param rec.levels a vector with the desired recall levels (\code{def:} \code{from:0.1}, \code{to:0.9}, \code{by:0.1})
#' @param folds number of folds on which computing the AUROC. If \code{folds=NULL} (\code{def.}), the one-shot AUROC is computed, 
#' otherwise the AUROC is computed averaged across folds.
#' @param seed initialization seed for the random generator to create folds. Set \code{seed} only if \code{folds}\eqn{\neq}\code{NULL}.
#' If \code{seed}\eqn{\neq}\code{NULL} and \code{folds=NULL}, the AUROC averaged across folds is computed without seed initialization.
#' @return a list with two elements:
#' \enumerate{
#' \item average: a vector with the the average precisions at different recall levels across classes
#' \item per.class: a matrix with the precisions at different recall levels: rows are classes, columns precisions at different recall level
#' }
#' @export
#' @examples
#' data(labels);
#' data(scores);
#' data(graph);
#' root <- root.node(g);
#' L <- L[,-which(colnames(L)==root)];
#' S <- S[,-which(colnames(S)==root)];
#' rec.levels <- seq(from=0.1, to=1, by=0.1);
#' PXR <- PXR.at.multiple.recall.levels.over.classes(L, S, rec.levels=rec.levels, folds=5, seed=23);
PXR.at.multiple.recall.levels.over.classes <- function(target, pred, rec.levels=seq(from=0.1, to=1, by=0.1), folds=NULL, seed=NULL){
	if(!is.matrix(target) || !is.matrix(pred))
		stop("PXR.at.multiple.recall.levels.over.classes: target or pred must be a matrix", call.=FALSE);	
	n.examples <- nrow(target);
	n.classes <- ncol(target);
	if((n.examples!=nrow(pred)) || (n.classes!=ncol(pred)))
		stop("PXR.at.multiple.recall.levels.over.classes: number of rows or columns do not match between target and predicted classes", call.=FALSE);
	if(any((target!=0) & (target!=1)))
		stop("PXR.at.multiple.recall.levels.over.classes: labels variable must take values 0 or 1", call.=FALSE);
	if(is.null(folds) && !is.null(seed))
		seed <- NULL;
	if(!is.null(folds) && is.null(seed))
		warning("PXR.at.multiple.recall.levels.over.classes: folds are generated without seed initialization", call.=FALSE);	
	## PXR averaged across folds
	if(!is.null(folds)){
		PXR <- c();
		for(i in 1:ncol(target)){
			y.target <- target[,i];
			indices <- 1:length(y.target);
			positives <- which(y.target==1);
			foldIndex <- do.stratified.cv.data.single.class(indices, positives, kk=folds, seed=seed);
			testIndex <- mapply(c, foldIndex$fold.positives, foldIndex$fold.negatives, SIMPLIFY=FALSE);
			fold.pxr <- list();
			for(k in 1:folds){
				fold.pxr[[k]] <- precision.at.multiple.recall.level(pred[testIndex[[k]]], target[testIndex[[k]]], rec.levels=rec.levels)$precision;
			}
			PXR <- rbind(PXR, Reduce("+", fold.pxr)/folds);
			rownames(PXR)[i] <- colnames(target)[i];
		}
		avgPXR <- apply(PXR,2,mean);
		res <- list(avgPXR=avgPXR, PXR=PXR);
		return(res);
	}
	## PXR one-shot
	res <- precision.at.multiple.recall.level.over.classes(target=target, predicted=pred, rec.levels=rec.levels);
	return(res);
}

