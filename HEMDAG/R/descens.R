##*********##
## DESCENS ##
##*********##
#' @name DESCENS
#' @aliases descens.threshold
#' @aliases descens.threshold.free
#' @aliases descens.weighted.threshold.free
#' @aliases descens.weighted.threshold
#' @aliases descens.tau
#' @seealso \code{\link{TPR-DAG}}
#' @title DESCENS variants 
#' @description The novelty of DESCENS with respect to TPR-DAG algorithm consists in considering the contribution of all the descendants 
#' of each node instead of only that of its children, since with the TPR-DAG algorithm the contribution of the descendants of a given node 
#' decays exponentially with their distance from the node itself, thus reducing the impact of the predictions made at the most specific levels of the ontology. 
#' On the contrary DESCENS predictions are more influenced by the information embedded in the most specific terms of the taxonomy (e.g. leaf nodes), 
#' thus putting more emphasis on the terms that most characterize the gene under study. 
#' @details The \emph{vanilla} DESCENS adopts a per-level bottom-up traversal of the DAG to correct the flat predictions \eqn{\hat{y}_i}:
#' \deqn{
#' 	\bar{y}_i := \frac{1}{1 + |\Delta_i|} (\hat{y}_i + \sum_{j \in \Delta_i} \bar{y}_j)
#' }
#' where \eqn{\Delta_i} are the positive descendants of \eqn{i}.
#' Different strategies to select the positive descendants \eqn{\Delta_i} can be applied:
#' \enumerate{
#' 	\item \strong{Threshold-Free} strategy: as positive descendants we choose those nodes that achieve a score higher than that of their ancestor node \eqn{i}:
#' 	\deqn{
#' 		\Delta_i := \{ j \in descendats(i) | \bar{y}_j > \hat{y}_i \}
#' 	}
#' 	\item \strong{Threshold} strategy: the positive descendants are selected on the basis of a threshold that can ben selected in two different ways:
#' 	\enumerate{
#' 		\item for each node a constant threshold \eqn{\bar{t}} is a priori selected:
#'		\deqn{
#'			\phi_i := \{ j \in descendats(i) | \bar{y}_j > \bar{t} \}
#'		}
#' 		For instance if the predictions represent probabilities it could be meaningful to a priori select \eqn{\bar{t}=0.5}.
#' 		\item the threshold is selected to maximize some performance metric \eqn{\mathcal{M}} estimated on the training data, as for instance
#' 		the F-score or the AUPRC. In other words the threshold is selected to maximize some measure of accuracy of the predictions 
#' 		\eqn{\mathcal{M}(j,t)} on the training data for the class \eqn{j} with respect to the threshold \eqn{t}. 
#' 		The corresponding set of positives \eqn{\forall i \in V} is:
#' 		\deqn{
#' 			\phi_i := \{ j \in descendants(i) | \bar{y}_j > t_j^*,  t_j^* = \arg \max_{t} \mathcal{M}(j,t) \}
#' 		}
#' 		For instance \eqn{t_j^*} can be selected from a set of \eqn{t \in (0,1)} through internal cross-validation techniques.
#'	}
#' }
#' The weighted DESCENS variants can be simply designed by adding a weight \eqn{w \in [0,1]} to balance the contribution between 
#' the prediction of the classifier associated with the node \eqn{i} and that of its positive descendants:
#' \deqn{
#' 	\bar{y}_i := w \hat{y}_i + \frac{(1 - w)}{|\Delta_i|} \sum_{j \in \phi_i} \bar{y}_j
#' }
#' The DESCENS-\eqn{\tau} variants balances the contribution between the positives children of a node \eqn{i} and that of
#' its positives descendants excluding the children by adding a weight \eqn{\tau \in [0,1]}:
#' \deqn{
#' \bar{y}_i := \frac{\tau}{ 1 +|\phi_i|} ( \hat{y}_i + \sum_{j \in \phi_i} \bar{y}_j ) + \frac{1-\tau}{1+|\delta_i|} ( \hat{y}_i + \sum_{j\in \delta_i} \bar{y}_j )
#' }
#' where \eqn{\phi_i} are the positive children of \eqn{i} and \eqn{\delta_i=\Delta_i \setminus \phi_i} the descendants of \eqn{i} without its children. 
#' If \eqn{\tau=1} we consider only the contribution of the positive children of \eqn{i}; if \eqn{\tau=0} only the descendants that are not
#' children contribute to the score, while for intermediate values of \eqn{\tau} we can balance the contribution of \eqn{\phi_i} and 
#' \eqn{\delta_i} positive nodes.
#' @param S a named flat scores matrix with examples on rows and classes on columns
#' @param g a graph of class \code{graphNEL}. It represents the hierarchy of the classes
#' @param root name of the class that it is on the top-level of the hierarchy (\code{def. root="00"})
#' @param t threshold for the choice of the positive descendants (\code{def. t=0.5}); whereas in the \code{descens.tau} variant 
#' the parameter \code{t} balances the contribution between the positives children of a node \eqn{i} and that of its
#' positives descendants excluding the positives children
#' @param w weight to balance between the contribution of the node \eqn{i} and that of its positive descendants
#' @return a named matrix with the scores of the classes corrected according to the DESCENS algorithm.
#' @examples
#' data(graph);
#' data(scores);
#' data(labels);
#' root <- root.node(g);
#' S.descensTF <- descens.threshold.free(S,g,root);
#' S.descensT <- descens.threshold(S,g,root,t=0.5);
#' S.descensW <- descens.weighted.threshold.free(S,g,root,w=0.5);
#' S.descensWT <- descens.weighted.threshold(S,g,root,w=0.5, t=0.5);
#' S.descensTAU <- descens.tau(S,g,root, t=0.5);

#' @rdname DESCENS
#' @export 
descens.threshold <- function(S, g, root="00", t=0.5){
	if(!(root %in% colnames(S))){
		max.score <- max(S);
		z <- rep(max.score,nrow(S));
		S <- cbind(z,S);
		colnames(S)[1] <- root;
	}
	## check consistency between nodes of g and classes of S
	class.check <- ncol(S)!=numNodes(g);
	if(class.check)
		stop("DESCENS: the number of nodes of the graph and the number of classes of the flat scores matrix does not match", call.=FALSE);
	## compute graph levels
	levels <- graph.levels(g,root);
	# bottom-up visit
	desc.bup <- build.descendants.bottom.up(g,levels);
	nodes <- names(desc.bup);
		for(i in 1:length(desc.bup)){
			if(length(desc.bup[[i]])!=1){
				node.curr <- nodes[i];
				parent <- S[,names(desc.bup[i])];
				tmp <- setdiff(desc.bup[[i]],node.curr);
				# if(i<=10) cat(node.curr, length(tmp), length(desc.bup[[i]]), "****************", sep="\n");
				desc <- as.matrix(S[,tmp]);
				# colnames(desc) <- desc.bup[[i]]
				for(j in 1:length(parent)){
					desc.set <- desc[j,] > t;    # positive descendants selection
					desc.pos <- desc[j,][desc.set];
					parent[j] <- (parent[j] + sum(desc.pos))/(1+length(desc.pos));   # flat scores correction
				}
				S[,names(desc.bup[i])] <- parent;
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

#' @rdname DESCENS
#' @export 
descens.threshold.free <- function(S, g, root="00"){
	if(!(root %in% colnames(S))){
		max.score <- max(S);
		z <- rep(max.score,nrow(S));
		S <- cbind(z,S);
		colnames(S)[1] <- root;
	}
	## check consistency between nodes of g and classes of S
	class.check <- ncol(S)!=numNodes(g);
	if(class.check)
		stop("DESCENS: the number of nodes of the graph and the number of classes of the flat scores matrix does not match", call.=FALSE);
	## compute graph levels
	levels <- graph.levels(g,root);
	# bottom-up visit
	desc.bup <- build.descendants.bottom.up(g,levels);
	nodes <- names(desc.bup);
		for(i in 1:length(desc.bup)){
			if(length(desc.bup[[i]])!=1){
				node.curr <- nodes[i];
				parent <- S[,names(desc.bup[i])];
				tmp <- setdiff(desc.bup[[i]],node.curr);
				# if(i<=10) cat(node.curr, length(tmp), length(desc.bup[[i]]), "****************", sep="\n");
				desc <- as.matrix(S[,tmp]);
				# colnames(desc) <- desc.bup[[i]]
				for(j in 1:length(parent)){
					desc.set <- desc[j,] > parent[j];	# positive descendants selection
					desc.pos <- desc[j,][desc.set];
					parent[j] <- (parent[j] + sum(desc.pos))/(1+length(desc.pos));   # flat scores correction
				}
				S[,names(desc.bup[i])] <- parent;
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

#' @rdname DESCENS
#' @export 
descens.weighted.threshold.free <- function(S, g, root="00", w=0.5){
	if(!(root %in% colnames(S))){
		max.score <- max(S);
		z <- rep(max.score,nrow(S));
		S <- cbind(z,S);
		colnames(S)[1] <- root;
	}
	## check consistency between nodes of g and classes of S
	class.check <- ncol(S)!=numNodes(g);
	if(class.check)
		stop("DESCENS: the number of nodes of the graph and the number of classes of the flat scores matrix does not match", call.=FALSE);
	## compute graph levels
	levels <- graph.levels(g,root);
	# bottom-up visit
	desc.bup <- build.descendants.bottom.up(g,levels);
	nodes <- names(desc.bup);
		for(i in 1:length(desc.bup)){
			if(length(desc.bup[[i]])!=1){
				node.curr <- nodes[i];
				parent <- S[,names(desc.bup[i])];
				tmp <- setdiff(desc.bup[[i]],node.curr);
				# if(i<=10) cat(node.curr, length(tmp), length(desc.bup[[i]]), "****************", sep="\n");
				desc <- as.matrix(S[,tmp]);
				# colnames(desc) <- desc.bup[[i]]
				for(j in 1:length(parent)){
					desc.set <- desc[j,] > parent[j];
					desc.pos <- desc[j,][desc.set];
					if(length(desc.pos)!=0){
						parent[j] <- w*parent[j] + (1-w)*sum(desc.pos)/length(desc.pos);  # flat scores correction
					}
				}
				S[,names(desc.bup[i])] <- parent;
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

#' @rdname DESCENS
#' @export 
descens.weighted.threshold <- function(S, g, root="00", t=0.5, w=0.5){
	if(!(root %in% colnames(S))){
		max.score <- max(S);
		z <- rep(max.score,nrow(S));
		S <- cbind(z,S);
		colnames(S)[1] <- root;
	}
	## check consistency between nodes of g and classes of S
	class.check <- ncol(S)!=numNodes(g);
	if(class.check)
		stop("DESCENS: the number of nodes of the graph and the number of classes of the flat scores matrix does not match", call.=FALSE);
	## compute graph levels
	levels <- graph.levels(g,root);
	# bottom-up visit
	desc.bup <- build.descendants.bottom.up(g,levels);
	nodes <- names(desc.bup);
		for(i in 1:length(desc.bup)){
			if(length(desc.bup[[i]])!=1){
				node.curr <- nodes[i];
				parent <- S[,names(desc.bup[i])];
				tmp <- setdiff(desc.bup[[i]],node.curr);
				# if(i<=10) cat(node.curr, length(tmp), length(desc.bup[[i]]), "****************", sep="\n");
				desc <- as.matrix(S[,tmp]);
				# colnames(desc) <- desc.bup[[i]]
				for(j in 1:length(parent)){
					desc.set <- desc[j,] > t;
					desc.pos <- desc[j,][desc.set];
					if(length(desc.pos)!=0){
						parent[j] <- w*parent[j] + (1-w)*sum(desc.pos)/length(desc.pos);  # flat scores correction
					}
				}
				S[,names(desc.bup[i])] <- parent;
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

#' @rdname DESCENS
#' @export 
descens.tau <- function(S, g, root="00", t=0.5){
	if(!(root %in% colnames(S))){
		max.score <- max(S);
		z <- rep(max.score,nrow(S));
		S <- cbind(z,S);
		colnames(S)[1] <- root;
	}
	## check consistency between nodes of g and classes of S
	class.check <- ncol(S)!=numNodes(g);
	if(class.check)
		stop("DESCENS: the number of nodes of the graph and the number of classes of the flat scores matrix does not match", call.=FALSE);
	## compute graph levels
	levels <- graph.levels(g,root);
	# bottom-up visit
	chd.bup <- get.children.bottom.up(g,levels);
	desc.bup <- build.descendants.bottom.up(g,levels);
	nodes <- names(desc.bup);
	for(i in 1:length(desc.bup)){
		if(length(desc.bup[[i]])!=1){ 		# in the desc list is included also the node itself
			node.curr <- nodes[i];
			parent <- S[,names(desc.bup[i])];
			tmp <- setdiff(desc.bup[[i]], node.curr);		# descendants
			delta <- setdiff(tmp, chd.bup[[i]]);  			# descendants without children 
			children <- as.matrix(S[,chd.bup[[i]]]);		# genes considering children node 
			desc <-  as.matrix(S[,delta]);					# genes considering descendants nodes without children
			for(j in 1:length(parent)){
				desc.set <- desc[j,] > parent[j];			# positive descendants (without children) selection
				desc.pos <- desc[j,][desc.set];
				child.set <- children[j,] > parent[j];  	# positive children selection
				child.pos <- children[j,][child.set];
				parent[j] <- t * ((parent[j] + sum(child.pos))/(1+length(child.pos))) + (1- t) * ((parent[j] + sum(desc.pos))/(1+length(desc.pos)));
			}
			S[,names(desc.bup[i])] <- parent;
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
