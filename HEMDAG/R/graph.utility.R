##**************************************************##
## Utility functions to process and analyze graphs	##
##**************************************************##

#' @title Build Graph Levels 
#' @description This function groups a set of nodes in according to their maximum depth in the graph. It first inverts the weights 
#' of the graph and then applies the Bellman Ford algorithm to find the shortest path, achieving in this way the longest path
#' @param g an object of class \code{graphNEL} 
#' @param root name of the root node (\code{def. root="00"}) 
#' @return a list of the nodes grouped w.r.t. the distance from the root: the first element of the list corresponds to the root node (level 0),
#' the second to nodes at maximum distance 1 (level 1), the third to the node at maximum distance 3 (level 2) and so on.
#' @export 
#' @examples
#' data(graph);
#' root <- root.node(g);
#' lev <- graph.levels(g, root=root);
graph.levels <- function(g, root="00"){
	if(sum(nodes(g) %in% root)==0) 
		stop("root node not found in g. Insert the root node");
	if(root.node(g)!=root) 
		stop("root is not the right root node of g. Use the function root.node(g) to find the root node of g");
	ed <- edges(g);
	ew <- edgeWeights(g);
	for(i in 1:length(ed)){
		l <- length(ew[[i]]);
		if(l!=0)
			ew[[i]][1:l] <- -1;
	}
	edL <- vector(mode="list", length=length(ed));
	names(edL) <- names(ed);
	for(i in 1:length(ed)){
		edL[[i]] <- list(edges=ed[[i]], weights=ew[[i]]);
	}
	G <- graphNEL(nodes=nodes(g), edgeL=edL, edgemode="directed");  
	depth.G <- bellman.ford.sp(G,root)$distance;
	depth.G <- -depth.G  
	levels <- vector(mode="list", length=max(depth.G)+1);
	names(levels) <- paste(rep("level", max(depth.G)+1), 0:max(depth.G), sep="_");
	for(i in 1:(max(depth.G)+1)){
		levels[[i]] <- names(which(depth.G==i-1));  
	}
	return(levels);
}

#' @title Flip Graph
#' @description Compute a directed graph with edges in the opposite direction
#' @param g a \code{graphNEL} directed graph
#' @return a graph (as an object of class \code{graphNEL}) with edges in the opposite direction w.r.t. g
#' @export
#' @examples
#' data(graph);
#' g.flipped <- compute.flipped.graph(g);
compute.flipped.graph <- function(g){
	ed <- edges(g);
	ndL <- vector(mode="list", length=length(ed));
	names(ndL) <- names(ed);
	for(i in 1:length(ed)){
		children <- ed[[i]];
		parent   <- names(ed[i]);
		if(length(children)!=0){
			for(j in 1:length(children)){
				ndL[[children[j]]] <- c(ndL[[children[j]]],parent); 
			}
		}
	}
	for (i in 1:length(ndL)){
		ndL[[i]] <- list(edges=ndL[[i]]);
	}
	og <- graphNEL(nodes=nodes(g), edgeL=ndL, edgemode="directed");
	return(og);
}

#' @name parents
#' @aliases get.parents
#' @aliases get.parents.top.down
#' @aliases get.parents.bottom.up
#' @aliases get.parents.topological.sorting
#' @title Build parents 
#' @description Compute the parents for each node of a graph
#' @param g a graph of class \code{graphNEL}. It represents the hierarchy of the classes
#' @param root name of the root node (\code{def. root="00"})
#' @param levels a list of character vectors. Each component represents a graph level and the elements of any 
#' component correspond to nodes. The level 0 coincides with the root node.
#' @seealso \code{\link{graph.levels}}
#' @examples
#' data(graph);
#' root <- root.node(g)
#' parents <- get.parents(g, root=root);
#' lev <- graph.levels(g, root=root);
#' parents.tod <- get.parents.top.down(g, lev, root=root);
#' parents.bup <- get.parents.bottom.up(g, lev, root=root);
#' parents.tsort <- get.parents.topological.sorting(g, root=root);

#' @rdname parents
#' @return \code{get.parents} returns a named list of character vectors. Each component corresponds to a node \eqn{x} of the graph (i.e. child node) 
#' and its vector is the set of its parents (the root node is not included)
#' @export
get.parents <- function(g, root="00"){
	if(sum(nodes(g) %in% root)==0) 
		stop("root node not found in g. Insert the root node");
	if(root.node(g)!=root) 
		stop("root is not the right root node of g. Use the function root.node(g) to find the root node of g");
	nd <- nodes(g)
	ndL <- vector(mode="list", length=length(nd));
	names(ndL) <- nd;
	ed <- edges(g);
	for(i in 1:length(ed)){
		children <- ed[[i]];
		parent   <- names(ed[i]);
		if(length(children)!=0){
			for(j in 1:length(children)){
				ndL[[children[j]]] <- c(ndL[[children[j]]],parent); 
			}
		}
	}
	ndL <- ndL[-which(names(ndL)==root)]; 
	return(ndL);
}

#' @rdname parents
#' @return \code{get.parents.top.down} returns a named list of character vectors. Each component corresponds to a node 
#' \eqn{x} of the graph (i.e. child node) and its vector is the set of its parents. 
#' The nodes order follows the levels of the graph from root (excluded) to leaves.
#' @export
get.parents.top.down <- function(g,levels, root="00"){
	if(sum(nodes(g) %in% root)==0) 
		stop("root node not found in g. Insert the root node");
	if(root.node(g)!=root) 
		stop("root is not the right root node of g. Use the function root.node(g) to find the root node of g");	ord.nd <- unlist(levels); 
	ndL <- vector(mode="list", length=length(ord.nd));
	names(ndL) <- ord.nd;
	ed <- edges(g);
	for(i in 1:length(ed)){
		children <- ed[[i]];
		parent   <- names(ed[i]);
		if(length(children)!=0){
			for(j in 1:length(children)){
				ndL[[children[j]]] <- c(ndL[[children[j]]],parent); 
			}
		}
	}
	ndL <- ndL[-which(names(ndL)==root)]; 
	return(ndL);
}

#' @rdname parents
#' @return \code{get.parents.bottom.up} returns a named list of character vectors. Each component corresponds to a node \eqn{x} of the 
#' graph (i.e. child node) and its vector isthe set of its parents. The nodes are ordered from leaves to root (excluded).
#' @export
get.parents.bottom.up <- function(g,levels, root="00"){
	if(sum(nodes(g) %in% root)==0) 
		stop("root node not found in g. Insert the root node");
	if(root.node(g)!=root) 
		stop("root is not the right root node of g. Use the function root.node(g) to find the root node of g");
	flip.ord.nd <- rev(unlist(levels)); 
	ndL <- vector(mode="list", length=length(flip.ord.nd));
	names(ndL) <- flip.ord.nd;
	ed <- edges(g);
	for(i in 1:length(ed)){
		children <- ed[[i]];
		parent   <- names(ed[i]);
		if(length(children)!=0){
			for(j in 1:length(children)){
				ndL[[children[j]]] <- c(ndL[[children[j]]],parent); 
			}
		}
	}
	ndL <- ndL[-which(names(ndL)==root)]; 
	return(ndL);
}

#' @rdname parents
#' @return \code{get.parents.topological.sorting} a named list of character vectors. Each component corresponds to a 
#' node \eqn{x} of the graph (i.e. child node) and its vector is the set of its parents. The nodes are ordered according to a 
#' topological sorting, i.e. parents node come before children node.
#' @export
get.parents.topological.sorting <- function(g, root="00"){
	if(sum(nodes(g) %in% root)==0) 
		stop("root node not found in g. Insert the root node");
	if(root.node(g)!=root) 
		stop("root is not the right root node of g. Use the function root.node(g) to find the root node of g");	ord.nd <- tsort(g); 
	ndL <- vector(mode="list", length=length(ord.nd));
	names(ndL) <- ord.nd;
	ed <- edges(g);
	for(i in 1:length(ed)){
		children <- ed[[i]];
		parent   <- names(ed[i]);
		if(length(children)!=0){
			for(j in 1:length(children)){
				ndL[[children[j]]] <- c(ndL[[children[j]]],parent); 
			}
		}
	}
	ndL <- ndL[-which(names(ndL)==root)]; 
	return(ndL);
}

#' @name descendants
#' @aliases build.descendants
#' @aliases build.descendants.per.level 
#' @aliases build.descendants.bottom.up
#' @title Build descendants 
#' @description Compute the descendants for each node of a graph
#' @param g a graph of class \code{graphNEL}. It represents the hierarchy of the classes
#' @param levels a list of character vectors. Each component represents a graph level and the elements of any 
#' component correspond to nodes. The level 0 coincides with the root node.
#' @seealso \code{\link{graph.levels}}
#' @examples
#' data(graph);
#' root <- root.node(g);
#' desc <- build.descendants(g);
#' lev <- graph.levels(g, root=root);
#' desc.tod <- build.descendants.per.level(g,lev);
#' desc.bup <- build.descendants.bottom.up(g,lev);

#' @rdname descendants
#' @return \code{build.descendants} returns a named list of vectors. 
#' Each component corresponds to a node \eqn{x} of the graph, and its vector is the set of its descendants including also \eqn{x}.
#' @export
build.descendants <- function(g){
	name.nodes <- nodes(g);
	g2 <- transitive.closure(g);
	desc <- edges(g2);
	for(x in name.nodes){
		desc[[x]] <- c(desc[[x]],x);
	}
	return(desc);
}

#' @rdname descendants
#' @return \code{build.descendants.per.level} returns a named list of vectors. 
#' Each component corresponds to a node \eqn{x} of the graph and its vector is the set of its descendants including also \eqn{x}.
#' The nodes are ordered from root (included) to leaves.
#' @export
build.descendants.per.level <- function(g,levels){
	ord.nd <- unlist(levels);
	g2 <- transitive.closure(g);
	desc <- edges(g2)[ord.nd];
	for(x in ord.nd){
		desc[[x]] <- c(desc[[x]],x);
	}
	return(desc);
}

#' @rdname descendants
#' @return \code{build.descendants.bottom.up} returns a named list of vectors. 
#' Each component corresponds to a node \eqn{x} of the graph and its vector is the set of its descendants including also \eqn{x}.
#' The nodes are ordered from leaves to root (included).
#' @export
build.descendants.bottom.up <- function(g,levels) {
	flip.ord.nd <- rev(unlist(levels));
	g2 <- transitive.closure(g);
	desc <- edges(g2)[flip.ord.nd];
	for(x in flip.ord.nd)
		desc[[x]] <- c(desc[[x]],x);
	return(desc);
}

#' @name children
#' @aliases build.children
#' @aliases get.children.top.down
#' @aliases get.children.top.down
#' @title Build children 
#' @description Compute the children for each node of a graph
#' @param g a graph of class \code{graphNEL}. It represents the hierarchy of the classes
#' @param levels a list of character vectors. Each component represents a graph level and the elements of any 
#' component correspond to nodes. The level 0 coincides with the root node. 
#' @seealso \code{\link{graph.levels}}
#' @examples
#' data(graph);
#' root <- root.node(g);
#' children <- build.children(g);
#' lev <- graph.levels(g, root=root);
#' children.tod <- get.children.top.down(g,lev);
#' children.bup <- get.children.bottom.up(g,lev);

#' @rdname children
#' @return \code{build.children} returns a named list of vectors. Each component corresponds to a node \eqn{x} of the graph and its 
#' vector is the set of its children
#' @export
build.children <- function(g){
	return(edges(g));
}

#' @rdname children
#' @return \code{get.children.top.down} returns a named list of character vectors. Each component corresponds to a node \eqn{x} 
#' of the graph (i.e. parent node) and its vector is the set of its children. The nodes are ordered from root (included) to leaves.
#' @export
get.children.top.down <- function(g,levels){
	child <- build.children(g)
	nd <- c();
	for(i in 1:length(levels)){
		level.nodes <- levels[[i]];
		nd <- append(nd,child[level.nodes]);
	}
	return(nd);
}

#' @rdname children
#' @return \code{get.children.bottom.up} returns a named list of character vectors. Each component corresponds to a node \eqn{x} 
#' of the graph (i.e. parent node) and its vector is the set of its children. The nodes are ordered from leaves (included) to root.
#' @export
get.children.bottom.up <- function(g,levels){
	flip.ord.nd <- rev(unlist(levels));
	ed <- edges(g);  
	nd <- ed[flip.ord.nd];
	return(nd);
}

#' @name ancestors
#' @aliases build.ancestors
#' @aliases build.ancestors.per.level
#' @aliases build.ancestors.bottom.up
#' @title Build ancestors 
#' @description Compute the ancestors for each node of a graph
#' @param g a graph of class \code{graphNEL}. It represents the hierarchy of the classes
#' @param levels a list of character vectors. Each component represents a graph level and the elements of any 
#' component correspond to nodes. The level 0 coincides with the root node.
#' @seealso \code{\link{graph.levels}}
#' @examples
#' data(graph);
#' root <- root.node(g);
#' anc <- build.ancestors(g);
#' lev <- graph.levels(g, root=root);
#' anc.tod <-build.ancestors.per.level(g,lev);
#' anc.bup <- build.ancestors.bottom.up(g,lev);

#' @rdname ancestors
#' @return \code{build.ancestos} returns a named list of vectors. Each component corresponds to a node \eqn{x} of the graph and its vector 
#' is the set of its ancestors including also \eqn{x}.
#' @export
build.ancestors <- function(g){
	og <- compute.flipped.graph(g);
	names.nodes <- nodes(og);
	og2 <- transitive.closure(og);
	anc <- edges(og2);
	for(x in names.nodes){
		anc[[x]] <- c(anc[[x]],x);
	}
	return(anc);
}

#' @rdname ancestors
#' @return \code{build.ancestors.per.level} returns a named list of vectors. Each component corresponds to a node \eqn{x} 
#' of the graph and its vector is the set of its ancestors including also \eqn{x}. The nodes are ordered from root (included) to leaves.
#' @export
build.ancestors.per.level <- function(g,levels){
	og <- compute.flipped.graph(g);
	ord.nd <- unlist(levels);
	og2 <- transitive.closure(og);
	anc <- edges(og2)[ord.nd];
	for(x in ord.nd)
		anc[[x]] <- c(anc[[x]],x);
	return(anc);
}

#' @rdname ancestors
#' @return \code{build.ancestors.bottom.up} a named list of vectors. Each component corresponds to a node \eqn{x} of the 
#' graph and its vector is the set of its ancestors including also \eqn{x}. The nodes are ordered from leaves to root (included).
#' @export
build.ancestors.bottom.up <- function(g,levels){
	og <- compute.flipped.graph(g);
	flip.ord.nd <- rev(unlist(levels));
	og2 <- transitive.closure(og);
	anc <- edges(og2)[flip.ord.nd];
	for(x in flip.ord.nd)
		anc[[x]] <- c(anc[[x]],x);
	return(anc);
}

#' @title Root node
#' @description Find the root node of a directed graph
#' @param g a graph of class \code{graphNEL}. It represents the hierarchy of the classes
#' @return name of the root node
#' @export
#' @examples
#' data(graph);
#' root <- root.node(g);
root.node <- function(g){
	d <- degree(g);
	root <- names(which(d$inDegree==0));
	return(root);
}

#' @title Leaves
#' @description Find the leaves of a directed graph
#' @param g a graph of class \code{graphNEL}. It represents the hierarchy of the classes
#' @return a vector with the names of the leaves of g
#' @export
#' @examples
#' data(graph);
#' leaves <- find.leaves(g);
find.leaves <- function(g){
	d <- degree(g);
	leaves <- names(which(d$outDegree==0));
	return(leaves);
}

#' @title Distances from leaves
#' @description This function returns the minimum distance of each node from one of the leaves of the graph
#' @param g a graph of class \code{graphNEL}. It represents the hierarchy of the classes
#' @return a named vector. The names are the names of the nodes of the graph \code{g}, and their values represent the distance from the leaves.
#' A value equal to 0 is assigned to the leaves, 1 to nodes with distance 1 from a leaf and so on
#' @export
#' @examples
#' data(graph);
#' dist.leaves <- distances.from.leaves(g);
distances.from.leaves <- function(g){
	leaves <- find.leaves(g);
	n.leaves <- length(leaves);
	og <- compute.flipped.graph(g);
	og <- addNode("root", og)
	#   for (x in leaves)
	#    og = addEdge("root", x, og, 1);
	og <- addEdge(rep("root",n.leaves), leaves, og, rep(1,n.leaves));
	dist <- acc(og,"root")[[1]]-1;
	return(dist);
}

#' @title Constraints Matrix
#' @description This function returns a matrix with two columns and as many rows as there are edges.
#' The entries of the first columns are the index of the node the edge cames from (i.e. children nodes), 
#' the entries of the second columns indicate the index of node the edge is to (i.e. parents nodes). 
#' Referring to a DAG this matrix defines a partial order. 
#' @param g a graph of class \code{graphNELL}. It represents the hierarchy of the classes
#' @return a constraints matrix w.r.t the graph \code{g}
#' @export
#' @examples
#' data(graph);
#' m <- constraints.matrix(g);
constraints.matrix <- function(g){
	eM <- edgeMatrix(g);
	eM <- cbind(eM[2,],eM[1,]);
	nd <- nodes(g);
	dimnames(eM) <- list(nd[eM[,2]], c("child","parent"))
	return(eM);
}

#' @title Weighted Adjacency Matrix
#' @description Construct a Weighted Adjacency Matrix (wadj matrix) of a graph
#' @param file name of the plain text file to be read (\code{def. edges}). The format of the file is a sequence of rows. 
#' Each row corresponds to an edge represented through a pair of vertices separated by blanks and the weight of the edges.\cr
#' For instance: \code{nodeX nodeY score}
#' @param compressed boolean value:
#' \itemize{
#'	\item TRUE (def.): the input file must be in a .gz compressed format;
#'	\item FALSE: the input file must be in a plain text format;
#' }
#' @param nodename boolean value:
#' \itemize{
#' 	\item TRUE (def.): the names of nodes are gene symbol (i.e. characters);
#' 	\item FALSE: the names of the nodes are entrez gene ID (i.e. integer numbers);
#' }
#' @details The input paramenter nodename sorts the row names of the wadj matrix in increasing order if they are integer number or 
#' in alphabetic order if they are characters.
#' @return a named symmetric weighted adjacency matrix of the graph
#' @export
#' @examples
#' edges <- system.file("extdata/edges.txt", package="HEMDAG");
#' W <- weighted.adjacency.matrix(file=edges, compressed=FALSE, nodename=TRUE);
weighted.adjacency.matrix <- function(file="edges.txt", compressed=TRUE, nodename=TRUE){
	if(compressed){
		m <- read.table(gzfile(file), colClasses="character", stringsAsFactors=FALSE);
	}else{
		m <- as.matrix(read.table(file, colClasses="character", stringsAsFactors=FALSE));
	}
	if(nodename){
		nodes <- sort(unique(as.vector(as.matrix(m[,1:2])))); ##NB:df must be converted as matrix to make as.vector workig..
	}else{
		nodes <- as.character(sort(as.numeric(unique(as.vector(m[,1:2]))))); 
	}
	n.nodes <- length(nodes);
	# building the adjacency matrix
	W <- matrix(0, nrow=n.nodes, ncol=n.nodes);
	dimnames(W) <- list(nodes,nodes);
	W[cbind(m[,1], m[,2])] <- as.numeric(m[,3]);
	W[cbind(m[,2], m[,1])] <- as.numeric(m[,3]);
	return(W);
}

#' @title Tupla Matrix
#' @description Trasform a Weighted Adjacency Matrix (wadj matrix) of a graph in a tupla, i.e. as a sequences of rows separated by 
#' blank and the weight of the edges, e.g \code{nodeX nodeY score}
#' @param m a weighetd adjacency matrix of the graph. Rows and columns are examples. It must be a square named matrix.
#' @param compressed boolean value:
#' \itemize{
#'	\item TRUE (def.): the output file will be a \code{gz} compressed format;
#'	\item FALSE: the output file will be a plain \code{text} format;
#' }
#' @param output.file name of the file of the  to be written
#' @details Only the \emph{non-zero} interactions are kept, while the \emph{zero} interactions are discarded. 
#' In other words in the \code{output.file} are reported only those nodes having a weight different from zero 
#' @return if \code{compressed=TRUE} the weighted adjacency matrix as tupla is stored in a compressed \code{gz}, 
#' otherwise (\code{compressed=FALSE}) it is stored in a plain \code{text} file. 
#' @export
#' @examples
#' \dontrun{
#' data(wadj);
#' tupla.matrix(W,compressed=TRUE, output.file="tupla.wadj.gz");
#' tupla.matrix(W,compressed=FALSE, output.file="tupla.wadj.txt");}
tupla.matrix <- function(m, compressed=TRUE, output.file="net.file.gz"){
	im <- which(m!=0, arr.ind=TRUE);
	rows <- rownames(im);
	colrep.names <- intersect(colnames(m), rownames(im));
	colrep.times <- table(im[,2])
	cols <- rep(colrep.names, times=colrep.times);
	df <- data.frame(row=rows, col=cols, score=m[im]);
	if(compressed){
		write.table(df, file=gzfile(output.file), quote=FALSE, row.names=FALSE, col.names=FALSE);
	}else{
		write.table(df, file=output.file, quote=FALSE, row.names=FALSE, col.names=FALSE);
	}
}

#' @title HPO specific annotations matrix
#' @description Construct the labels matrix of the most specific HPO terms
#' @details The input plain text file representing the most specific associations gene-HPO term can be obtained by forking the GitHub repositoty
#' \href{https://github.com/gecko515/HPOparser}{HPOparser}, a collection of \code{Perl} subroutines to parse the HPO OBO file 
#' and the HPO annotations file.
#' @param file text file representing the most specific associations gene-HPO term \cr (\code{def: "gene2pheno.txt"}). 
#' The file must be written as sequence of rows. Each row represents a gene and all its
#' associations with abnormal phenotype tab separated, \cr \emph{e.g.: gene_1 <tab> phen1 <tab> ... phen_N}. \cr
#' See \strong{Details} section to know more information about how to obtain this file.
#' @param genename boolean value: 
#' \itemize{
#' 	\item TRUE (def.): the names of genes are \emph{gene symbol} (i.e. characters);
#' 	\item FALSE: the names of gene are entrez \emph{gene ID} (i.e. integer numbers);
#' }
#' @return the annotation matrix of the most specific annotations (0/1): rows are genes and columns are HPO terms.
#' Let's denote \eqn{M} the labels matrix. If \eqn{M[i,j]=1}, means that the gene \eqn{i} is annotated with the class \eqn{j}, otherwise \eqn{M[i,j]=0}.
#' @export
#' @examples
#' gene2pheno <- system.file("extdata/gene2pheno.txt", package="HEMDAG");
#' spec.ann <- specific.annotation.matrix(gene2pheno, genename=TRUE);
specific.annotation.matrix <- function(file="gene2pheno.txt", genename="TRUE"){
	line <- readLines(file);
	tmp <- strsplit(line, split="[ +\t]");

	gene.names <- c();
	for(i in 1:length(tmp)) gene.names <- c(gene.names,tmp[[i]][1]);

	ann.list <- list();
	for(i in 1:length(tmp)) ann.list[[i]] <- tmp[[i]][-1];
	names(ann.list) <- gene.names;

	hpoID <- unique(unlist(ann.list));
 
	n.genes <- length(gene.names);
	n.hpoID <- length(hpoID);
	m <- matrix(integer(n.genes * n.hpoID), nrow=n.genes);
	rownames(m) <- gene.names;
	colnames(m) <- hpoID;

	for (i in gene.names){
		spec.ann <- ann.list[[i]]; 
		m[i, spec.ann] <- 1;  
	}

	if(genename){
		m <- m[sort(rownames(m)),sort(colnames(m))];
	}else{
		rname <- as.character(sort(as.numeric(gene.names)));
		m <- m[rname, sort(colnames(m))];
	}

	return(m);
}

#' @title Specific annotations list
#' @description Construct a list of the most specific annotations starting from the table of the most specific annotations
#' @param ann annotation matrix (0/1). Rows are examples and columns are most specific terms. It must be a named matrix. 
#' @return a named list, where the names of each component correspond to an examples (genes) and the elements of each component 
#' are the most specific classes associated to that genes
#' @seealso \code{\link{specific.annotation.matrix}}
#' @export
#' @examples
#' data(labels);
#' spec.list <- specific.annotation.list(L);
specific.annotation.list <- function(ann){
 ann.list <- apply(ann, 1, function(gene){
		terms <- which(gene==1);
		return(names(gene[terms]));
	});
	return(ann.list);
}

#' @title Transitive closure of annotations 
#' @description Performs the transitive closure of the annotations using ancestors and the most specific annotation table.
#' The annotations are propagated from bottom to top, enriching the most specific annotations table.
#' The rows of the matrix correspond to the genes of the most specific annotation table and the columns to the HPO terms/classes
#' @param ann.spec the annotation matrix of the most specific annotations (0/1): rows are genes and columns are HPO terms.
#' @param anc list of the ancestors of the ontology. 
#' @return an annotation table T: rows correspond to genes and columns to HPO terms. \eqn{T[i,j]=1} means that gene \eqn{i} is annotated for the term \eqn{j},
#' \eqn{T[i,j]=0} means that gene \eqn{i} is not annotated for the term \eqn{j}.
#' @seealso \code{\link{specific.annotation.matrix}}, \code{\link{build.ancestors}}
#' @export
#' @examples
#' data(graph);
#' data(labels);
#' anc <- build.ancestors(g);
#' tca <- transitive.closure.annotations(L, anc);
transitive.closure.annotations <- function(ann.spec, anc){
	## costructiion of annotation list
	ann.list <- specific.annotation.list(ann.spec);

	## cotruction the full empty annotation matrix
	entrezIDs <- rownames(ann.spec);
	n.genes <- length(entrezIDs);
	HPOIDs <- names(anc);
	n.HPOID <- length(anc);
	hpo.ann <- matrix(numeric(n.HPOID * n.genes), nrow=n.genes, ncol=n.HPOID);	#empty label matrix
	dimnames(hpo.ann) <- list(entrezIDs,HPOIDs);
	
	## fill the full empty annotation matrix with the most specific annotation 	
	hpo.spec.term <- colnames(ann.spec); # the most specific hpo terms
	# might happen that there are same HPO IDs that are classified as "obsolete" in obo file, but that still exist in the annotation file 
	hpo.spec.term.sel <- HPOIDs[HPOIDs  %in% hpo.spec.term]; # removing obsolete HPO terms...
	hpo.ann[entrezIDs,hpo.spec.term.sel] <- ann.spec[,hpo.spec.term.sel];	

	## transitive closure: annotation propagation from the most specific nodes to all its ancestors
	for (i in entrezIDs){
		spec.ann <- ann.list[[i]];
		all.anc <- lapply(spec.ann, function(x) return(anc[[x]]));
		all.anc <- unique(unlist(all.anc));
		hpo.ann[i, all.anc] <- 1;  # setting the annotations derived by transitive closure
	}

	## remove HPO empty terms 
	hpo.ann <- hpo.ann[,colSums(hpo.ann)!=0];
	return(hpo.ann);
}

#' @title Full annotations matrix
#' @description Construct a full annotations table using ancestors and the most specific annotations table w.r.t. a given weighted adjacency matrix (wadj). 
#' The rows of the full annotations matrix correspond to all the examples of the given weighted adjacency matrix and the columns to the class/terms.
#' The transitive closure of the annotations is performed. 
#' @details The examples present in the annotation matrix (\code{ann.spec}) but not in the adjacency weighted matrix (\code{W}) are purged.
#' @param W symmetric adjacency weighted matrix of the graph 
#' @param anc list of the ancestors of the ontology. 
#' @param ann.spec the annotation matrix of the most specific annotations (0/1): rows are genes and columns are classes.
#' @return a full annotation table T, that is a matrix in which the transitive closure of annotations was performed. 
#' Rows correspond to genes of the weighted adjiacency matrix and columns to terms. 
#' \eqn{T[i,j]=1} means that gene \eqn{i} is annotated for the term \eqn{j}, \eqn{T[i,j]=0} means that gene \eqn{i} is not annotated for the term \eqn{j}.
#' @seealso \code{\link{weighted.adjacency.matrix}}, \code{\link{build.ancestors}}, \cr
#' \code{\link{specific.annotation.matrix}}, \code{\link{transitive.closure.annotations}}
#' @export
#' @examples
#' data(wadj);
#' data(graph);
#' data(labels);
#' anc <- build.ancestors(g);
#' full.ann <- full.annotation.matrix(W, anc, L);
full.annotation.matrix <- function(W, anc, ann.spec){
	## construction of annotation list
	ann.list <- specific.annotation.list(ann.spec);

	## construction the full empty annotation matrix
	entrezIDs <- rownames(W);
	n.genes <- length(entrezIDs);
	HPOIDs <- names(anc);
	n.HPOID <- length(anc);
	hpo.ann <- matrix(numeric(n.HPOID * n.genes), nrow=n.genes, ncol=n.HPOID);	#empty label matrix
	dimnames(hpo.ann) <- list(entrezIDs,HPOIDs);
	
	## fill the full empty annotation matrix with the most specific annotation 
	entrezIDs2hpo <- rownames(ann.spec);								# all genes that are associated with hpo terms
	entrezIDs.sel <- entrezIDs[entrezIDs %in% entrezIDs2hpo];			# genes 2 hpo terms 2 entrez id of wadj
	hpo.spec.term <- colnames(ann.spec);								# the most specific hpo terms
	#might happen that there are same HPO IDs that are classified as "obsolete" in obo file, but that still exist in the annotation file (e.g. build 1233)
	hpo.spec.term.sel <- HPOIDs[HPOIDs  %in% hpo.spec.term]; # removing obsolete HPO terms...
	hpo.ann[entrezIDs.sel,hpo.spec.term.sel] <- ann.spec[entrezIDs.sel,hpo.spec.term.sel];	# setting the most specific annotations

	## transitive closure: annotation propagation from the most specific nodes to all its ancestors
	for (i in entrezIDs){
		spec.ann <- ann.list[[i]];
		all.anc <- lapply(spec.ann, function(x) return(anc[[x]]));
		all.anc <- unique(unlist(all.anc));
		hpo.ann[i, all.anc] <- 1;  # setting the annotations derived by transitive closure
	}

	## remove HPO empty terms 
	hpo.ann <- hpo.ann[,colSums(hpo.ann)!=0];
	return(hpo.ann);
}

#' @title Do full annotations matrix
#' @description High-level function to obtain a full annotation matrix, that is a matrix in which the transitive closure of annotations was performed, 
#' respect to a given weighted adiacency matrix
#' @param anc.file.name name of the file containg the list for each node the list of all its ancestor (without \code{rda} extension)
#' @param anc.dir relative path to directory where the ancestor file is stored 
#' @param net.file name of the file containing the weighted adjiacency matrix of the graph (without \code{rda} extension)
#' @param net.dir relative path to directory where the weighted adjiacency matrix is stored 
#' @param ann.file.name name of the file containing the matrix of the most specific annotations (without \code{rda} extension)
#' @param ann.dir relative path to directory where the matrix of the most specific annotation is stored 
#' @param output.name name of the output file without rda extension (without \code{rda} extension)
#' @param output.dir relative path to directory where the output file must be stored   
#' @return a full annotation matrix T, that is a matrix in which the transitive closure of annotations was performed.
#' Rows correspond to genes of the input weighted adjiacency matrix and columns to terms. 
#' \eqn{T[i,j]=1} means that gene \eqn{i} is annotated for the term \eqn{j}, \eqn{T[i,j]=0} means that gene \eqn{i} is not annotated for the term \eqn{j}.
#' @seealso \code{\link{full.annotation.matrix}}
#' @export
#' @examples
#' data(graph);
#' data(labels);
#' data(wadj);
#' if (!dir.exists("data")){
#' 	dir.create("data");
#' }
#' if (!dir.exists("results")){
#' 	dir.create("results");
#' }
#' anc <- build.ancestors(g);
#' save(anc,file="data/ancestors.rda");
#' save(g,file="data/graph.rda");
#' save(L,file="data/labels.rda");
#' save(W,file="data/wadj.rda");
#' anc.dir <- net.dir <- ann.dir <- "data/";
#' output.dir <- "results/";
#' anc.file.name <- "ancestors";
#' net.file <- "wadj";
#' ann.file.name <- "labels";
#' output.name <- "full.ann.matrix";
#' Do.full.annotation.matrix(anc.file.name=anc.file.name, anc.dir=anc.dir, net.file=net.file, 
#'	net.dir=net.dir, ann.file.name=ann.file.name, ann.dir=ann.dir, output.name=output.name, 
#' 	output.dir=output.dir);
Do.full.annotation.matrix <- function(anc.file.name=anc.file.name, anc.dir=anc.dir, net.file=net.file, net.dir=net.dir, 
	ann.file.name=ann.file.name, ann.dir=ann.dir, output.name=output.name, output.dir=output.dir){
	## loading list of ancestors
	anc.path <- paste0(anc.dir, anc.file.name, ".rda");
	anc.name <- load(anc.path);
	anc <- eval(parse(text=anc.name));  

	## loading wadj matrix
	net.path <- paste0(net.dir, net.file, ".rda");
	net.name <- load(net.path);
	W <- eval(parse(text=net.name));  

	## loading the specific annotation matrix
	ann.path <- paste0(ann.dir, ann.file.name, ".rda");
	ann.name <- load(ann.path);
	ann.spec <- eval(parse(text=ann.name));

	## costruction of full HPO annotation matrix
	ann <- full.annotation.matrix(W=W, anc=anc, ann.spec=ann.spec);
	
	## saving labels matrix
	ann.file <- paste0(output.dir, output.name, ".rda");
	save(ann, file=ann.file, compress=TRUE);
}

#' @title Build submatrix
#' @title Build an annotation matrix with only those terms having more than n annotations.
#' @description Terms having less than n annotations are pruned. Terms having exactly n annotations are discarded as well.
#' @param hpo.ann the annotations matrix (0/1). Rows are examples and columns are classes 
#' @param n integer number of annotations to be pruned
#' @return Matrix of annotations having only those terms with more than n annotations
#' @export
#' @examples
#' data(labels);
#' subm <- do.submatrix(L,5);
do.submatrix <- function(hpo.ann,n){
	hpo.ann.sel <- hpo.ann[,colSums(hpo.ann)>n];
	return(hpo.ann.sel);
}

#' @title Build subgraph 
#' @description This function returns a subgraph with only the supplied nodes and any edges between them
#' @param nd a vector with the nodes for which the subgraph must be built
#' @param g a graph of class \code{graphNEL}. It represents the hierarchy of the classes
#' @param edgemode can be "directed" or "undirected"
#' @return a subgraph with only the supplied nodes 
#' @export
#' @examples
#' data(graph);
#' anc <- build.ancestors(g);
#' nd <- anc[["HP:0001371"]];
#' subg <- do.subgraph(nd, g, edgemode="directed");
do.subgraph <- function(nd, g, edgemode="directed"){
	ed <- edges(g);
	ed.sel <- ed[nd];

	ndL <- vector(mode="list", length=length(ed.sel));
	names(ndL) <- names(ed.sel);

	for(i in 1:length(ed.sel)){ 
		parent   <- names(ed.sel[i]);
		children <- ed.sel[[i]];
		if(length(children!=0)){ 
			children.map <- children[children %in% nd]
			ndL[[i]] <- append(ndL[[i]],children.map);
		} 
	}

	for (i in 1:length(ndL)){
		ndL[[i]] <- list(edges=ndL[[i]]);
	}

	G <- graphNEL(nodes=nd, edgeL=ndL, edgemode=edgemode);
	return(G);
}

#' @title Annotation matrix checker
#' @description This function assess the integrity of an annotation table in which a transitive closure of annotations was performed
#' @param anc list of the ancestors of the ontology. 
#' @param ann.spec the annotation matrix of the most specific annotations (0/1): rows are genes and columns are terms. 
#' @param hpo.ann the full annotations matrix (0/1), that is the matrix in which the transitive closure of the annotation was performed.
#' Rows are examples and columns are classes. 
#' @return If the transitive closure of the annotations is well performed "OK" is returned, otherwise a message error is printed on the stdout
#' @seealso \code{\link{build.ancestors}}, \code{\link{transitive.closure.annotations}}, \code{\link{full.annotation.matrix}}
#' @export
#' @examples
#' data(graph);
#' data(labels);
#' anc <- build.ancestors(g);
#' tca <- transitive.closure.annotations(L, anc);
#' check.annotation.matrix.integrity(anc, L, tca);
check.annotation.matrix.integrity <- function(anc, ann.spec, hpo.ann){
	## construction of annotation list
	ann.list <- specific.annotation.list(ann.spec);
	genes <- rownames(hpo.ann);

	check <- c();
	for (i in genes){
		spec.ann <- which(hpo.ann[i,]==1);
		len.ann <- length(spec.ann);
		all.anc <- lapply(ann.list[[i]], function(x) return(anc[[x]]));
		all.anc <- unique(unlist(all.anc));
		len.anc <- length(all.anc);
		cmp <- len.anc == len.ann;
		if(cmp==TRUE){
			check <- c(check,"OK");
		} else {
			check <- c(check,"NOTOK");
		}
	}
	names(check) <- genes;

	violated <- any(check!="OK");
	if(violated){
		n <- names(check)[check=="NOTOK"];
		cat("check.annotation.matrix: NOT_OK. Transitive closure NOT RESPECTED", "\n");
	}else{
		cat("check.annotation.matrix: OK", "\n");	
	}
}

#' @title DAG checker
#' @description This function assess the integrity of a DAG
#' @param g a graph of class \code{graphNEL}. It represents the hierarchy of the classes
#' @param root name of the class that is on the top-level of the hierarchy (def:"00")
#' @return If there are nodes not accessible from the root "OK" is printed, 
#' otherwise a message error and the list of the not accessible nodes is printed on the stdout
#' @export
#' @examples
#' data(graph);
#' root <- root.node(g);
#' check.DAG.integrity(g, root=root);
check.DAG.integrity <- function(g, root="00"){
	if(sum(nodes(g) %in% root)==0) 
		stop("root node not found in g. Insert the root node");
	if(root.node(g)!=root) 
		stop("the supplied root node is not the right root node of g. Use the function root.node(g) to find the root node of g");
  all.nodes <- nodes(g);
  acc.nodes <- names(acc(g,root)[[1]]);
  if((length(all.nodes) - length(acc.nodes)) > 1) {
		n <- setdiff(all.nodes,c(acc.nodes,root));
		cat("check.GO.integrity: not all nodes accessible from root", "\n");
		cat("Nodes not accessible from root: \n");
		cat(n,"\n");
	}else{ 
		cat("OK \n")
	};
}

#' @name hierarchical.checkers
#' @aliases check.hierarchy
#' @aliases check.hierarchy.single.sample
#' @title Hierarchical constraints cheker
#' @description Check if the true path rule is violated or not. In other words this function checks if the score of a 
#' parent or an ancestor node is always larger or equal than that of its children or descendants nodes
#' @param y.hier vector of scores relative to a single example. This must be a named numeric vector
#' @param S.hier the matrix with the scores of the classes corrected in according to hierarchy. This must be a named matrix: rows are
#' examples and colums are classes
#' @param g a graph of class \code{graphNEL}. It represents the hierarchy of the classes
#' @param root name of the class that it is the top-level (root) of the hierarchy (\code{def:00})
#' @examples
#' data(graph);
#' data(scores);
#' root <- root.node(g);
#' S.hier <- htd(S,g,root);
#' S.hier.single.example <- S.hier[sample(ncol(S.hier),1),];
#' check.hierarchy.single.sample(S.hier.single.example, g, root=root);
#' check.hierarchy(S.hier, g, root);

#' @return return a list of 3 elements:
#' \itemize{
#'  \item Status: 
#'	 \itemize{
#'	  \item OK if none hierarchical constraints have bee broken;
#'	  \item NOTOK if there is at least one hierarchical constraints broken;
#'   }
#' 	\item Hierarchy_Constraints_Broken:
#'	 \itemize{
#'	  \item TRUE: example did not respect the hierarchical constraints; 
#'	  \item FALSE: example broke the hierarchical constraints;
#'   }
#'  \item Hierarchy_costraints_satisfied: how many terms satisfied the hierarchical constraint
#' }
#' @export
check.hierarchy.single.sample <- function(y.hier,g, root="00"){
	if(!(root %in% names(y.hier))){
		max.score <- max(y.hier);
		y.hier <- c(max.score,y.hier);
		names(y.hier)[1] <- root;
	}
	par <- get.parents(g,root);
	v <- c()
	for(i in 1:length(par)){
		child <- y.hier[names(par[i])];
		parents <- y.hier[par[[i]]]
		x <- parents >= child   
		y <- any(x==0)  
		v <- append(v,y)
	}
	names(v) <- names(par)
	violated <- any(v==TRUE);
	if(violated)
		Status = "NOTOK"
	else
		Status = "OK";
	h <- as.factor(v);
	k <- summary(h);
	l <- list(Status=Status, hierarchy.constraints.broken=v, hierarchy.costraints.satisfied=k);
	return(l);
}

#' @rdname hierarchical.checkers
#' @export
check.hierarchy <- function(S.hier,g, root="00"){
	if(!(root %in% colnames(S.hier))){
		max.score <- max(S.hier);
		z <- rep(max.score,nrow(S.hier));
		S.hier <- cbind(z,S.hier);
		colnames(S.hier)[1] <- root;
	}
	par <- get.parents(g,root);
	v <- c()
	for(i in 1:length(par)){
		child <- S.hier[,names(par[i])];
		parents <- S.hier[,par[[i]]]
		x <- parents >= child   
		y <- any(x==0)   
		v <- append(v,y)
	}
	names(v) <- names(par)
	violated <- any(v==TRUE);
	if(violated)
		Status = "NOTOK"
	else
		Status = "OK";
	h <- as.factor(v);
	k <- summary(h);
	l <- list(Status=Status, hierarchy.constraints.broken=v, hierarchy.costraints.satisfied=k);
	return(l);
}

#' @title Unstratified cross-validation
#' @description This function splits a dataset in k-fold in an unstratified way (that is a fold may not have an equal amount of positive and 
#' negative examples). This function is used to perform k-fold cross-validation experiments in a hierarchical correction contest where 
#' splitting dataset in a stratified way is not needed. 
#' @param S matrix of the flat scores. It must be a named matrix, where rows are example (e.g. genes) and columns are classes/terms (e.g. HPO terms)
#' @param kk number of folds in which to split the dataset (\code{def. k=5})
#' @param seed seed for the random generator. If \code{NULL} (def.) no initialization is performed
#' @return a list with \eqn{k=kk} components (folds). Each component of the list is a character vector contains the names of the examples.
#' @export
#' @examples
#' data(scores);
# folds <- do.unstratified.cv.data(S, kk=5, seed=23);
do.unstratified.cv.data <- function(S, kk=5, seed=NULL){
	set.seed(seed);
	examples <- rownames(S);
	n <- nrow(S);
	size <- c();
	folds <- vector(mode="list", length=kk)
	names(folds) <- paste0(rep("fold",kk), 1:kk)
	for (k in 1:kk) {
		first <- ((k - 1) * n) %/% kk
		last <- (k * n) %/% kk
		size <- last-first;
		x	<- sample(examples,size);
		folds[[k]] <- x;
		examples <- setdiff(examples,x);
	}
	return(folds);
}

#' @name stratified.cross.validation
#' @aliases do.stratified.cv.data.single.class
#' @aliases do.stratified.cv.data.over.classes
#' @title Stratified cross validation
#' @description Generate data for the stratified cross-validation 
#' @param labels labels matrix. Rows are genes and columns are classes. Let's denote \eqn{M} the labels matrix. 
#' If \eqn{M[i,j]=1}, means that the gene \eqn{i} is annotated with the class \eqn{j}, otherwise \eqn{M[i,j]=0}.
#' @param examples indices or names of the examples. Can be either a vector of integers or a vector of names. 
#' @param positives vector of integers or vector of names. The indices (or names) refer to the indices (or names) of 'positive' examples	
#' @param kk number of folds (\code{def=5})
#' @param seed seed of the random generator (\code{def=NULL}). If is set to \code{NULL} no initialization is performed
#' @examples
#' data(labels);
#' examples.index <- 1:nrow(L);
#' examples.name <- rownames(L);
#' positives <- which(L[,3]==1);
#' x <- do.stratified.cv.data.single.class(examples.index, positives, kk=5, seed=23);
#' y <- do.stratified.cv.data.single.class(examples.name, positives, kk=5, seed=23);
#' z <- do.stratified.cv.data.over.classes(L, examples.index, kk=5, seed=23);
#' k <- do.stratified.cv.data.over.classes(L, examples.name, kk=5, seed=23);

#' @rdname stratified.cross.validation
#' @return \code{do.stratified.cv.data.single.class} returns a list with 2 two component:
#' \itemize{
#'  \item fold.non.positives: a list with \eqn{k} components. Each component is a vector with the indices (or names) of the non-positive elements. 
#' 	Indices (or names) refer to row numbers (or names) of a data matrix.
#'	 \item fold.positives: a list with \eqn{k} components. Each component is a vector with the indices (or names) of the positive elements. 
#' 	Indices (or names) refer to row numbers (or names) of a data matrix.	 
#' }
#' @export
do.stratified.cv.data.single.class <- function(examples, positives, kk=5, seed=NULL){
	set.seed(seed);

	if(is.numeric(examples) && length(names(positives))!=0){
		positives <- unname(positives);
	}

	if(is.character(examples) && length(names(positives))!=0){
		positives <- names(positives);
	}

	negatives <- setdiff(examples,positives); 	
	n <- length(positives);		
	m <- length(negatives);		
	set.pos <- list();
	set.neg <- list();

	for (k in 1:kk) {
		#fold positives 
		last.pos <- (k * n) %/% kk;
		first.pos  <- ((k - 1) * n) %/% kk;
		size.pos <-  last.pos - first.pos;				
		subset.pos <- 	sample(positives, size.pos);			
		set.pos[[k]] <- subset.pos;						
		positives <- setdiff(positives, subset.pos);	
		
		#fold non positives
		last.neg <- (k * m) %/% kk;
		first.neg  <- ((k - 1) * m) %/% kk;
		size.neg <-  last.neg - first.neg;				
		subset.neg <- sample(negatives, size.neg);					
		set.neg[[k]] <- subset.neg;						
		negatives <- setdiff(negatives, subset.neg);	
	}
	return(list(fold.positives=set.pos, fold.negatives=set.neg));
}

#' @rdname stratified.cross.validation
#' @return \code{do.stratified.cv.data.over.classes} returns a list with \eqn{n} components, where \eqn{n} is the number of classes of the labels matrix. 
#' Each component \eqn{n} is in turn a list with \eqn{k} elements, where \eqn{k} is the number of folds. 
#' Each fold contains an equal amount of examples positives and negatives.
#' @export
do.stratified.cv.data.over.classes <- function(labels, examples, kk=5, seed=NULL){
	set.seed(seed);
	folds <- list();
	for(class in colnames(labels)){
		folds[[class]] <- list();
		positives <- which(labels[,class]==1);
		strfold <- do.stratified.cv.data.single.class(examples,positives, kk=kk, seed=seed);
		for(k in 1:kk){
			folds[[class]][[k]] <- list();
			names(folds[[class]])[k] <- paste0("fold",k);
			folds[[class]][[k]] <- append(strfold$fold.positives[[k]], strfold$fold.negatives[[k]]);
		}
	}
	return(folds);
}
