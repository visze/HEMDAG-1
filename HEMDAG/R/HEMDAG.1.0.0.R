#' @import  graph					   
#' @import  RBGL				   
#' @import  precrec				   
#' @import  PerfMeas			   
#' @import  preprocessCore
#' @import methods
#' @importFrom utils read.table write.table

## Quiet concerns of R CMD check. 
## Avoid this warning: no visible binding for global variable ‘curvetypes’
if(getRversion() >= "2.15.1") utils::globalVariables(c("curvetypes"))

##**************##
## IO FUNCTIONS ##
##**************##

#' @title Parse an HPO OBO file
#' @description Read an HPO OBO file (\href{http://human-phenotype-ontology.github.io/}{HPO}) and write 
#' the edges of the DAG on a plain text file.The format of the file is a sequence of
#' rows and each row corresponds to an edge represented through a pair of vertices separated by blanks
#' @param file an HPO OBO file
#' @param output.file name of the file of the edges to be written
#' @return a text file representing the edges in the format: source  destination (i.e. one row for each edge) 
#' @export
#' @examples
#' \dontrun{
#' hpobo <- "http://purl.obolibrary.org/obo/hp.obo";
#' do.edges.from.HPO.obo(file=hpobo, output.file="hp.edge");}
do.edges.from.HPO.obo <- function(file="hp.obo", output.file="edge.file"){
	line <- readLines(file);
	n.lines <- length(line);
	m <- matrix(character(1000000*2), ncol=2);
	colnames(m)=c("source", "destination");
	i=1;
	j=0; # number of edges;
	#browser();
	while(i<=n.lines){
		while((i<=n.lines) && (line[i]!="[Term]")){
			i <- i + 1;
		}
		if(i>=n.lines){break();}
		i <- i + 1; # id
		destination <- strsplit(line[i], split="[ +\t]")[[1]][2];
		while( (line[i]!="") && (strsplit(line[i], split="[ +\t]")[[1]][1]!="is_a:") ){ # checking first is_a entry
			i <- i + 1;
		}  
		if (line[i] == ""){next();}  # we are at the end of the record and is_a has been found
		source <- strsplit(line[i], split="[ +\t]")[[1]][2];
		j <- j + 1;
		i <- i + 1;
		m[j,]<-c(source,destination);
		while( (line[i]!="") && (strsplit(line[i], split="[ +\t]")[[1]][1]=="is_a:") ){# checking successive is_a entry
			source <- strsplit(line[i], split="[ +\t]")[[1]][2];
			i <- i + 1;
			j <- j + 1;
			m[j,]<-c(source,destination);
		} 	  
	}
	m <- m[1:j,];
	write.table(m, file=output.file, quote=FALSE, row.names=FALSE, col.names=FALSE);
}

#' @title Write a directed graph on file
#' @description An object of class \code{graphNEL} is read and the graph is written on a plain text file as sequence of rows 
#' @param g a graph of class \code{graphNEL}
#' @param file name of the file to be written
#' @return a plain text file representing the graph. Each row corresponds to an edge represented through a pair of vertices separated by blanks
#' @export
#' @examples
#' \dontrun{
#' data(graph);
#' write.graph(g, file="graph.edges.txt");}
write.graph <- function(g, file="graph.txt"){  
	num.edges <- length(unlist(edges(g)));
	num.v <- numNodes(g);
	m <- matrix(character(num.edges*2), ncol=2);
	res <- edges(g);
	count=0;
	node1 <- names(res);
	for (i in 1:num.v) {
		x <- res[[i]];
	len.x <- length(x);
	if (len.x!=0)
		for (j in 1:len.x) {
			count <- count + 1;
		m[count,] <- c(node1[i],x[j]);
		}
	}
	write.table(m, file=file, quote=FALSE, row.names=FALSE, col.names=FALSE);
}

#' @title Read a directed graph from a file
#' @description A directed graph is read from a file and a \code{graphNEL} object is built
#' @param file name of the file to be read. The format of the file is a sequence of rows and each row corresponds 
#' to an edge represented through a pair of vertices separated by blanks
#' @return an object of class \code{graphNEL}
#' @export
#' @examples
#' ed <- system.file("extdata/graph.edges.txt", package= "HEMDAG");
#' g <- read.graph(file=ed);
read.graph <- function(file="graph.txt"){  
	m <- as.matrix(read.table(file, colClasses="character"));
	thenodes<-sort(unique(as.vector(m))); # nodes
	n.nodes <- length(thenodes);
	n.edges <- nrow(m);
	# building the graph
	edL <- vector("list", length=n.nodes);
	names(edL) <- thenodes;
	for(i in 1:n.nodes)
		edL[[i]]<-list(edges=NULL);
	g <- graphNEL(nodes=thenodes, edgeL=edL, edgemode="directed");
	g <- addEdge(m[1:n.edges,1], m[1:n.edges,2], g, rep(1,n.edges));
	return(g);
}

#' @title Read an undirected graph from a file
#' @description The graph is read from a file and a \code{graphNEL} object is built. The format of the input file is a sequence of rows. 
#' Each row corresponds to an edge represented through a pair of vertices separated by blanks, and the weight of the edge.
#' @param file name of the file to be read
#' @return a graph of class \code{graphNEL}
#' @export
#' @examples
#' edges <- system.file("extdata/edges.txt" ,package="HEMDAG");
#' g <- read.undirected.graph(file=edges);
read.undirected.graph <- function(file="graph.txt") {  
	m <- as.matrix(read.table(file, colClasses="character"));
	thenodes<-sort(unique(as.vector(m[,1:2]))); # nodes
	n.nodes <- length(thenodes);
	n.edges <- nrow(m);
	# building the graph
	edL <- vector("list", length=n.nodes);
	names(edL) <- thenodes;
	for(i in 1:n.nodes)
		edL[[i]]<-list(edges=NULL);
	g <- graphNEL(nodes=thenodes, edgeL=edL, edgemode="undirected");
	g <- addEdge(m[1:n.edges,1], m[1:n.edges,2], g, as.numeric(m[1:n.edges,3]));
	return(g);
}

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
#' @param g a graph of class \code{graphNEL}. It represents the hierarchy of the classes.
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
#' @param g a graph of class \code{graphNEL}. It represents the hierarchy of the classes.
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
#' @param g a graph of class \code{graphNEL}. It represents the hierarchy of the classes.
#' @return a named vector. The names are the names of the nodes of the graph g, and their values represent the distance from the leaves.
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

#' @title Constraints matrix
#' @description This function returns a matrix with two columns and as many rows as there are edges.
#' The entries of the first columns are the index of the node the edge cames from (i.e. children nodes), 
#' the entries of the second columns indicate the index of node the edge is to (i.e. parents nodes). 
#' Referring to a DAG this matrix defines a partial order. 
#' @param g a graph of class \code{graphNEL}L. Represent the hierarchy of the class
#' @return a constraints matrix w.r.t the graph g
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
#' Each row corresponds to an edge represented through a pair of vertices separated by blanks and the weight of the edge.\cr
#' For instance: nodeX nodeY score
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

#' @title HPO specific annotations matrix
#' @description Construct the labels matrix of the most specific HPO terms
#' @param file text file representing the most specific associations gene-HPO term (\code{def: "gene2pheno.txt"}). 
#' The file must be written as sequence of rows. Each row represents a gene and all its
#' associations with abnormal phenotype tab separated, \cr \emph{e.g.: gene_1 <tab> phen1 <tab> ... phen_N}
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
		return(names(gene[terms])) 
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
#' @seealso \code{\link{weighted.adjacency.matrix}}, \code{\link{build.ancestors}}, \code{\link{specific.annotation.matrix}},
#' \code{\link{transitive.closure.annotations}}
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
#' @param g a graph of class \code{graphNEL}. It represents the hierarchy of the classes.
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
#' @param g a graph of class \code{graphNEL}. It represents the hierarchy of the classes.
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
#' @param S.hier the matrix with the scores of the classes corrected in according to hierarchy
#' @param g a graph of class \code{graphNEL}. It represents the hierarchy of the classes
#' @param root name of the class that it is the top-level (root) of the hierarchy (\code{def:00})
#' @seealso \code{\link{htd}}
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
#' @param norm.type can assume two character values:
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
	}else{
		## Quantile Normalization 
		## NOTE: normalize.quantiles function returns a unnamed matrix. colnames are essential for hier.corr..
		S.norm <- normalize.quantiles(S);
		dimnames(S.norm) <- list(rownames(S),colnames(S));
		S <- S.norm;
		rm(S.norm);	
	}

	## Storing results
	save(S, file=paste0(flat.norm.dir, norm.type, ".", flat.file, ".rda"));
}


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
		stop("AUPRC.single.class: lengths of true and predicted labels do not match.");
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
		stop("AUROC.single.class: lengths of true and predicted labels do not match.");
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

##*********##
## HTD-DAG ##
##*********##

#' @title HTD-DAG
#' @description Implementetion of a top-down procedure to correct the scores in the hierarchy according to the 
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
#' @param g a graph of class graphNEL. It represents the hierarchy of the classes
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
#' Furthermore it is possible to add a weight \eqn{w \in [0,1]} to balance between the contribution of the node \eqn{i} and that of its positive children
#' \eqn{\phi}, through their convex combination, obtaining in this way the weighted TPR-DAG version:
#' \deqn{
#' 	\bar{y}_i := w \hat{y}_i + \frac{(1 - w)}{|\phi_i|} \sum_{j \in \phi_i} \bar{y}_j
#' }
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
#' TPR.TF <- tpr.threshold.free(S,g,root);
#' TPR.T <- tpr.threshold(S,g,root,t=0.5);
#' TPR.W <- tpr.weighted.threshold.free(S,g,root,w=0.5);
#' TPR.WT <- tpr.weighted.threshold(S,g,root,w=0.5, t=0.5);

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
## DO HTD-DAG ##
##************##

#' @title HTD-DAG vanilla
#' @description High level function to compute hierarchical correction according to HTD-DAG algorithm
#' @param norm boolean value: 
#' \itemize{
#' \item \code{TRUE} (def.): the flat scores matrix has been already normalized in according to a normalization method;	
#' \item \code{FALSE}: the flat scores matrix has not been normalized yet. See the parameter \code{norm.type} for which normalization can be applied.
#' }
#' @param norm.type three values: 
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
#' for each class and example considered. It stored in \code{hierScore.dir} directory.
#' \item \code{PCM} (Protein Centric Measures) average and per-example: compute \code{F-score} measure by \code{find.best.f} function. 
#' It stored in \code{perf.dir} directory.
#' \item \code{PRC} (area under Precision-Recall Curve) average and per.class: compute \code{PRC} by \pkg{precrec} package. 
#' It stored in \code{perf.dir} directory.
#' \item \code{AUC} (Area Under ROC Curve) average and per-class: compute \code{AUC} by \pkg{precrec} package. 
#' It stored in \code{perf.dir} directory.
#' \item \code{PxR} (Precision at fixed Recall levels) average and per classes: compute \code{PxR}  by \pkg{PerfMeas} package. 
#' It stored in \code{perf.dir} directory.
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
#' Do.HTD(norm=FALSE, norm.type= "MaxNorm", flat.file=flat.file, ann.file=ann.file, 
#' dag.file=dag.file, flat.dir=flat.dir, ann.dir=ann.dir, dag.dir=dag.dir, 
#' flat.norm.dir=flat.norm.dir, n.round=3, f.criterion ="F", hierScore.dir=hierScore.dir, 
#' perf.dir=perf.dir);
Do.HTD <- function	(norm=TRUE, norm.type= "NONE", flat.file=flat.file, ann.file=ann.file, dag.file=dag.file, flat.dir=flat.dir, 
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
			dag.dir=dag.dir,flat.norm.dir=flat.norm.dir);
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
	AUC.htd <- AUROC.single.over.classes(ann, S); gc();

	## Hierarchical PxR at fixed recall levels 
	PXR.htd <- precision.at.multiple.recall.level.over.classes(ann, S); gc();

	## Computing Hierarchical Examples-Measures 
	FMM.htd <- find.best.f(ann, S, n.round=n.round, f.criterion =f.criterion, verbose=FALSE, b.per.example=TRUE);

	## Hierarchical PRC (average and per.class) computed by precrec package
	PRC.htd <- AUPRC.single.over.classes(ann, S); 

	## storing the hierarchical matrix
	S.htd <- S;
	rm(S); gc();

	## Storing Results #########
	if(norm){
		save(S.htd, file=paste0(hierScore.dir, flat.file, ".hierScores.htd.rda"), compress=TRUE);
		save(AUC.flat, AUC.htd, file=paste0(perf.dir, "AUC.", flat.file, ".hierScores.htd.rda"), compress=TRUE);
		save(PXR.flat, PXR.htd, file=paste0(perf.dir, "PXR.", flat.file, ".hierScores.htd.rda"), compress=TRUE);
		save(FMM.flat, FMM.htd, file=paste0(perf.dir, "PCM.", flat.file, ".hierScores.htd.rda"), compress=TRUE);
		save(PRC.flat, PRC.htd, file=paste0(perf.dir, "PRC.", flat.file, ".hierScores.htd.rda"), compress=TRUE);
	}else{
		save(S.htd, file=paste0(hierScore.dir, norm.type,".", flat.file, ".hierScores.htd.rda"), compress=TRUE);	
		save(AUC.flat, AUC.htd, file=paste0(perf.dir, "AUC.", norm.type,".", flat.file, ".hierScores.htd.rda"), compress=TRUE);	
		save(PXR.flat, PXR.htd, file=paste0(perf.dir, "PXR.", norm.type,".", flat.file, ".hierScores.htd.rda"), compress=TRUE);	
		save(FMM.flat, FMM.htd, file=paste0(perf.dir, "PCM.", norm.type,".", flat.file, ".hierScores.htd.rda"), compress=TRUE);
		save(PRC.flat, PRC.htd, file=paste0(perf.dir, "PRC.", norm.type,".", flat.file, ".hierScores.htd.rda"), compress=TRUE);
	}
}

#' @title HTD-DAG holdout
#' @description High level function to correct the scores with a hierarchy according to HTD-DAG algorithm applying a classical holdout procedure
#' @param norm boolean value: 
#' \itemize{
#' \item \code{TRUE} (def.): the flat scores matrix has been already normalized in according to a normalization method;	
#' \item \code{FALSE}: the flat scores matrix has not been normalized yet. See the parameter \code{norm} for which normalization can be applied.
#' }
#' @param norm.type three values: 
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
#' 	\item \code{hierarchical scores matrix}: a matrix with examples on rows and classes on columns representing the computed hierarchical scores 
#' 	for each class and example considered. It stored in \code{hierScore.dir} directory.
#' 	\item \code{PCM} (Protein Centric Measures) average and per-example: compute \code{F-score} measure by \code{find.best.f} function. 
#' 	It stored in \code{perf.dir} directory.
#' 	\item \code{PRC} (area under Precision-Recall Curve) average and per.class: compute \code{PRC} by \pkg{precrec} package. 
#' 	It stored in \code{perf.dir} directory.
#' 	\item \code{AUC} (Area Under ROC Curve) average and per-class: compute \code{AUC} by \pkg{precrec} package. 
#' 	It stored in \code{perf.dir} directory.
#' 	\item \code{PxR} (Precision at fixed Recall levels) average and per classes: compute \code{PxR}  by \pkg{PerfMeas} package. 
#' 	It stored in \code{perf.dir} directory.
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
#' Do.HTD.holdout(norm=FALSE, norm.type= "MaxNorm", flat.file=flat.file, ann.file=ann.file, 
#' dag.file=dag.file, ind.test.set=ind.test.set, ind.dir=ind.dir, flat.dir=flat.dir, 
#' ann.dir=ann.dir, dag.dir=dag.dir, flat.norm.dir=flat.norm.dir, n.round=3, f.criterion ="F", 
#' hierScore.dir=hierScore.dir, perf.dir=perf.dir);
Do.HTD.holdout <- function(norm=TRUE, norm.type= "NONE", flat.file=flat.file, ann.file=ann.file, dag.file=dag.file, 
	ind.test.set=ind.test.set, ind.dir=ind.dir, flat.dir=flat.dir, ann.dir=ann.dir, dag.dir=dag.dir, flat.norm.dir=NULL, 
	n.round=3, f.criterion ="F", hierScore.dir="hierScore.dir/", perf.dir="perf.dir/"){

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
			dag.dir=dag.dir,flat.norm.dir=flat.norm.dir);
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
	AUC.htd <- AUROC.single.over.classes(ann, S); 
		
	## Hierarchical PxR at fixed recall levels 
	PXR.htd <- precision.at.multiple.recall.level.over.classes(ann, S); 

	## Computing Hierarchical Examples-Measures 
	FMM.htd <- find.best.f(ann, S, n.round=n.round, f.criterion =f.criterion, verbose=FALSE, b.per.example=TRUE);	

	## Hierarchical PRC (average and per.class) computed by precrec package
	PRC.htd <- AUPRC.single.over.classes(ann, S);  

	## storing the hierarchical matrix
	S.htd <- S;
	rm(S); gc();
	
	## Storing Results #########
	if(norm){
		save(S.htd, file=paste0(hierScore.dir, flat.file, ".hierScores.htd.rda"), compress=TRUE);
		save(AUC.flat, AUC.htd, file=paste0(perf.dir, "AUC.", flat.file, ".hierScores.htd.rda"), compress=TRUE);
		save(PXR.flat, PXR.htd, file=paste0(perf.dir, "PXR.", flat.file, ".hierScores.htd.rda"), compress=TRUE);
		save(FMM.flat, FMM.htd, file=paste0(perf.dir, "PCM.", flat.file, ".hierScores.htd.rda"), compress=TRUE);
		save(PRC.flat, PRC.htd, file=paste0(perf.dir, "PRC.", flat.file, ".hierScores.htd.rda"), compress=TRUE);
	}else{
		save(S.htd, file=paste0(hierScore.dir, norm.type,".", flat.file, ".hierScores.htd.rda"), compress=TRUE);	
		save(AUC.flat, AUC.htd, file=paste0(perf.dir, "AUC.", norm.type,".", flat.file, ".hierScores.htd.rda"), compress=TRUE);	
		save(PXR.flat, PXR.htd, file=paste0(perf.dir, "PXR.", norm.type,".", flat.file, ".hierScores.htd.rda"), compress=TRUE);	
		save(FMM.flat, FMM.htd, file=paste0(perf.dir, "PCM.", norm.type,".", flat.file, ".hierScores.htd.rda"), compress=TRUE);
		save(PRC.flat, PRC.htd, file=paste0(perf.dir, "PRC.", norm.type,".", flat.file, ".hierScores.htd.rda"), compress=TRUE);
	}
}

##************##
## DO TPR-DAG ##
##************##

#' @name TPR.DAG.CV
#' @aliases Do.tpr.threshold.cv
#' @aliases Do.tpr.weighted.threshold.free.cv
#' @aliases Do.tpr.weighted.threshold.cv
#' @title TPR-DAG cross-validation 
#' @description High level function to correct the scores with a hierarchy according to one of the TPR-DAG variants performing cross-validation experiments
#' @details These high level functions perform a classical cross-validation procedure to find the best threshold maximizing on F-score measure.
#' @param threshold range of threshold values to be tested in order to find the best threshold (def: from:0.1, to:0.9, by:0.1 step).
#' The denser the range is, the higher the probability to find the best theshold is, but obviously the execution time will be increasing.
#' @param weight range of weight values to be tested in order to find the best weight (def: from:0.1, to:0.9, by:0.1 step).
#' The denser the range is, the higher the probability to find the best weight is, but obviously the execution time will be increasing.
#' @param kk number of folds of the cross validation (def: kk=5);
#' @param seed intialization seed for the random generator to create folds (def:0). If \code{NULL} (def.) no initialization is performed
#' @param norm boolean value: 
#' \itemize{
#' \item \code{TRUE} (def.): the flat scores matrix has been already normalized in according to a normalization method;	
#' \item \code{FALSE}: the flat scores matrix has not been normalized yet. See the parameter \code{norm.type} for which normalization can be applied.
#' }
#' @param norm.type three values: 
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
#' for each class and example considered. It stored in \code{hierScore.dir} directory.
#' \item \code{PCM} (Protein Centric Measures) average and per-example: compute \code{F-score} measure by \code{find.best.f} function. 
#' It stored in \code{perf.dir} directory.
#' \item \code{PRC} (area under Precision-Recall Curve) average and per.class: compute \code{PRC} by \pkg{precrec} package. 
#' It stored in \code{perf.dir} directory.
#' \item \code{AUC} (Area Under ROC Curve) average and per-class: compute \code{AUC} by \pkg{precrec} package. 
#' It stored in \code{perf.dir} directory.
#' \item \code{PxR} (Precision at fixed Recall levels) average and per classes: compute \code{PxR}  by \pkg{PerfMeas} package. 
#' It stored in \code{perf.dir} directory.
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
			dag.dir=dag.dir,flat.norm.dir=flat.norm.dir);
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
	S.tpr <- c();

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
		S.tpr <- rbind(S.tpr, pred.test);
	}
	## put the rows (i.e. genes) of assembled k sub-matrix in the same order of the full beginning matrix
	S.tpr <- S.tpr[rownames(S),];

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
	FMM.tpr <- list(average=FMM.avg, per.example=FMM.per.example);

	## averaging AUC (average and per.class) across k testing sets 
	AUC.average.tpr.over.test <- Reduce("+", AUC.average.test)/kk;
	AUC.class.tpr.over.test <- Reduce("+", AUC.class.test)/kk;
	AUC.tpr <- list(average=AUC.average.tpr.over.test, per.class=AUC.class.tpr.over.test);

	## averaging PxR (average and per.class) across k testing sets 
	PXR.average.tpr.over.test <- Reduce("+", PXR.average.test)/kk;
	PXR.class.tpr.over.test <- Reduce("+", PXR.class.test)/kk;
	PXR.tpr <- list(average=PXR.average.tpr.over.test, per.class=PXR.class.tpr.over.test);

	## averaging PRC (average and per.class) across k testing sets 
	PRC.average.tpr.over.test <- Reduce("+", PRC.average.test)/kk;
	PRC.class.tpr.over.test <- Reduce("+", PRC.class.test)/kk;
	PRC.tpr <- list(average=PRC.average.tpr.over.test, per.class=PRC.class.tpr.over.test);

	## Storing Results #########
	if(norm){
		save(S.tpr, file=paste0(hierScore.dir, flat.file, ".hierScores.tprT.rda"), compress=TRUE);
		save(AUC.flat, AUC.tpr, file=paste0(perf.dir, "AUC.", flat.file, ".hierScores.tprT.rda"), compress=TRUE);
		save(PXR.flat, PXR.tpr, file=paste0(perf.dir, "PXR.", flat.file, ".hierScores.tprT.rda"), compress=TRUE);
		save(FMM.flat, FMM.tpr, file=paste0(perf.dir, "PCM.", flat.file, ".hierScores.tprT.rda"), compress=TRUE);
		save(PRC.flat, PRC.tpr, file=paste0(perf.dir, "PRC.", flat.file, ".hierScores.tprT.rda"), compress=TRUE);
	}else{
		save(S.tpr, file=paste0(hierScore.dir, norm.type,".", flat.file, ".hierScores.tprT.rda"), compress=TRUE);
		save(AUC.flat, AUC.tpr, file=paste0(perf.dir, "AUC.", norm.type, ".", flat.file, ".hierScores.tprT.rda"), compress=TRUE);
		save(PXR.flat, PXR.tpr, file=paste0(perf.dir, "PXR.", norm.type, ".", flat.file, ".hierScores.tprT.rda"), compress=TRUE);
		save(FMM.flat, FMM.tpr, file=paste0(perf.dir, "PCM.", norm.type, ".", flat.file, ".hierScores.tprT.rda"), compress=TRUE);
		save(PRC.flat, PRC.tpr, file=paste0(perf.dir, "PRC.", norm.type, ".", flat.file, ".hierScores.tprT.rda"), compress=TRUE);
	}
}

#' @name TPR.DAG.HOLDOUT
#' @aliases Do.tpr.threshold.holdout
#' @aliases Do.tpr.weighted.threshold.free.holdout
#' @aliases Do.tpr.weighted.threshold.holdout
#' @title TPR-DAG holdout 
#' @description High level function to correct the scores with a hierarchy according to one of the TPR-DAG variants performing a classical holdout procedure.
#' @details These high level functions perform a classical holdout procedure to find the best threshold maximizing on F.score measure.
#' @param threshold range of threshold values to be tested in order to find the best threshold (def: from:0.1, to:0.9, by:0.1 step).
#' The denser the range is, the higher the probability to find the best theshold is, but obviously the execution time will be increasing.
#' @param weight range of weight values to be tested in order to find the best weight (def: from:0.1, to:0.9, by:0.1 step).
#' The denser the range is, the higher the probability to find the best weight is, but obviously the execution time will be increasing.
#' @param kk number of folds of the cross validation (def: kk=5);
#' @param seed intialization seed for the random generator to create folds (def:0). If \code{NULL} (def.) no initialization is performed
#' @param norm boolean value: 
#' \itemize{
#' \item \code{TRUE} (def.): the flat scores matrix has been already normalized in according to a normalization method;	
#' \item \code{FALSE}: the flat scores matrix has not been normalized yet. See the parameter \code{norm.type} for which normalization can be applied.
#' }
#' @param norm.type three values: 
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
#' for each class and example considered. It stored in \code{hierScore.dir} directory.
#' \item \code{PCM} (Protein Centric Measures) average and per-example: compute \code{F-score} measure by \code{find.best.f} function. 
#' It stored in \code{perf.dir} directory.
#' \item \code{PRC} (area under Precision-Recall Curve) average and per.class: compute \code{PRC} by \pkg{precrec} package. 
#' It stored in \code{perf.dir} directory.
#' \item \code{AUC} (Area Under ROC Curve) average and per-class: compute \code{AUC} by \pkg{precrec} package. 
#' It stored in \code{perf.dir} directory.
#' \item \code{PxR} (Precision at fixed Recall levels) average and per classes: compute \code{PxR}  by \pkg{PerfMeas} package. 
#' It stored in \code{perf.dir} directory.
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
			dag.dir=dag.dir,flat.norm.dir=flat.norm.dir);
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
	AUC.tpr <- AUROC.single.over.classes(ann.test, S.test);

	## PxR at fixed recall levels (average and per.class) computed with PerfMeas package on the test set
	PXR.tpr <- precision.at.multiple.recall.level.over.classes(ann.test, S.test);

	## F.measure: Computing Hierarchical Examples-Measures 
	FMM.tpr <- find.best.f(ann.test, S.test, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=TRUE);	

	## Hierarchical PRC (average and per.class) computed by precrec package
	PRC.tpr <- AUPRC.single.over.classes(ann.test, S.test); 

	## storing the hierarchical matrix
	S.tpr <- S.test;
	rm(S.test, S.training); gc();

	## Storing Results #########
	if(norm){
		save(S.tpr, file=paste0(hierScore.dir, flat.file, ".hierScores.tprT.rda"), compress=TRUE);
		save(AUC.flat, AUC.tpr, file=paste0(perf.dir, "AUC.", flat.file, ".hierScores.tprT.rda"), compress=TRUE);
		save(PXR.flat, PXR.tpr, file=paste0(perf.dir, "PXR.", flat.file, ".hierScores.tprT.rda"), compress=TRUE);
		save(FMM.flat, FMM.tpr, file=paste0(perf.dir, "PCM.", flat.file, ".hierScores.tprT.rda"), compress=TRUE);
		save(PRC.flat, PRC.tpr, file=paste0(perf.dir, "PRC.", flat.file, ".hierScores.tprT.rda"), compress=TRUE);
	}else{
		save(S.tpr, file=paste0(hierScore.dir, norm.type,".", flat.file, ".hierScores.tprT.rda"), compress=TRUE);
		save(AUC.flat, AUC.tpr, file=paste0(perf.dir, "AUC.", norm.type, ".", flat.file, ".hierScores.tprT.rda"), compress=TRUE);
		save(PXR.flat, PXR.tpr, file=paste0(perf.dir, "PXR.", norm.type, ".", flat.file, ".hierScores.tprT.rda"), compress=TRUE);
		save(FMM.flat, FMM.tpr, file=paste0(perf.dir, "PCM.", norm.type, ".", flat.file, ".hierScores.tprT.rda"), compress=TRUE);
		save(PRC.flat, PRC.tpr, file=paste0(perf.dir, "PRC.", norm.type, ".", flat.file, ".hierScores.tprT.rda"), compress=TRUE);
	}
}

#' @rdname TPR.DAG.CV
#' @export
Do.tpr.weighted.threshold.free.cv <- function(weight=seq(from=0.1, to=1, by=0.1), kk=5, seed=NULL, norm=TRUE, norm.type= "NONE", 
	flat.file=flat.file, ann.file=ann.file, dag.file=dag.file, flat.dir=flat.dir, ann.dir=ann.dir,dag.dir=dag.dir, 
	flat.norm.dir=NULL, n.round=3, f.criterion ="F", hierScore.dir="hierScore.dir/", perf.dir="perf.dir/" ){
	
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
			dag.dir=dag.dir,flat.norm.dir=flat.norm.dir);
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
	S.tpr <- c();

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
		S.tpr <- rbind(S.tpr, pred.test);
	}
	## put the rows (i.e. genes) of assembled k sub-matrix in the same order of the full beginning matrix
	S.tpr <- S.tpr[rownames(S),];

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
	FMM.tpr <- list(average=FMM.avg, per.example=FMM.per.example);

	## averaging macro-AUC (average and per.class) across k testing sets 
	AUC.average.tpr.over.test <- Reduce("+", AUC.average.test)/kk;
	AUC.class.tpr.over.test <- Reduce("+", AUC.class.test)/kk;
	AUC.tpr <- list(average=AUC.average.tpr.over.test, per.class=AUC.class.tpr.over.test);

	## averaging PxR (average and per.class) across k testing sets 
	PXR.average.tpr.over.test <- Reduce("+", PXR.average.test)/kk;
	PXR.class.tpr.over.test <- Reduce("+", PXR.class.test)/kk;
	PXR.tpr <- list(average=PXR.average.tpr.over.test, per.class=PXR.class.tpr.over.test);

	## averaging PRC (average and per.class) across k testing sets 
	PRC.average.tpr.over.test <- Reduce("+", PRC.average.test)/kk;
	PRC.class.tpr.over.test <- Reduce("+", PRC.class.test)/kk;
	PRC.tpr <- list(average=PRC.average.tpr.over.test, per.class=PRC.class.tpr.over.test);

	## Storing Results #########
	if(norm){
		save(S.tpr, file=paste0(hierScore.dir, flat.file, ".hierScores.tprW.rda"), compress=TRUE);
		save(AUC.flat, AUC.tpr, file=paste0(perf.dir, "AUC.", flat.file, ".hierScores.tprW.rda"), compress=TRUE);
		save(PXR.flat, PXR.tpr, file=paste0(perf.dir, "PXR.", flat.file, ".hierScores.tprW.rda"), compress=TRUE);
		save(FMM.flat, FMM.tpr, file=paste0(perf.dir, "PCM.", flat.file, ".hierScores.tprW.rda"), compress=TRUE);
		save(PRC.flat, PRC.tpr, file=paste0(perf.dir, "PRC.", flat.file, ".hierScores.tprW.rda"), compress=TRUE);
	}else{
		save(S.tpr, file=paste0(hierScore.dir, norm.type,".", flat.file, ".hierScores.tprW.rda"), compress=TRUE);
		save(AUC.flat, AUC.tpr, file=paste0(perf.dir, "AUC.", norm.type, ".", flat.file, ".hierScores.tprW.rda"), compress=TRUE);
		save(PXR.flat, PXR.tpr, file=paste0(perf.dir, "PXR.", norm.type, ".", flat.file, ".hierScores.tprW.rda"), compress=TRUE);
		save(FMM.flat, FMM.tpr, file=paste0(perf.dir, "PCM.", norm.type, ".", flat.file, ".hierScores.tprW.rda"), compress=TRUE);
		save(PRC.flat, PRC.tpr, file=paste0(perf.dir, "PRC.", norm.type, ".", flat.file, ".hierScores.tprW.rda"), compress=TRUE);
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
			dag.dir=dag.dir,flat.norm.dir=flat.norm.dir);
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
	AUC.tpr <- AUROC.single.over.classes(ann.test, S.test);	

	## PxR at fixed recall levels (average and per.class) computed with PerfMeas package on the test set
	PXR.tpr <- precision.at.multiple.recall.level.over.classes(ann.test, S.test);

	## F.measure: Computing Hierarchical Examples-Measures 
	FMM.tpr <- find.best.f(ann.test, S.test, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=TRUE);

	## Hierarchical PRC (average and per.class) computed by precrec package
	PRC.tpr <- AUPRC.single.over.classes(ann.test, S.test); 

	## storing the hierarchical matrix
	S.tpr <- S.test;
	rm(S.test, S.training);

	## Storing Results #########
	if(norm){
		save(S.tpr, file=paste0(hierScore.dir, flat.file, ".hierScores.tprW.rda"), compress=TRUE);
		save(AUC.flat, AUC.tpr, file=paste0(perf.dir, "AUC.", flat.file, ".hierScores.tprW.rda"), compress=TRUE);
		save(PXR.flat, PXR.tpr, file=paste0(perf.dir, "PXR.", flat.file, ".hierScores.tprW.rda"), compress=TRUE);
		save(FMM.flat, FMM.tpr, file=paste0(perf.dir, "PCM.", flat.file, ".hierScores.tprW.rda"), compress=TRUE);
		save(PRC.flat, PRC.tpr, file=paste0(perf.dir, "PRC.", flat.file, ".hierScores.tprW.rda"), compress=TRUE);
	}else{
		save(S.tpr, file=paste0(hierScore.dir, norm.type,".", flat.file, ".hierScores.tprW.rda"), compress=TRUE);
		save(AUC.flat, AUC.tpr, file=paste0(perf.dir, "AUC.", norm.type, ".", flat.file, ".hierScores.tprW.rda"), compress=TRUE);
		save(PXR.flat, PXR.tpr, file=paste0(perf.dir, "PXR.", norm.type, ".", flat.file, ".hierScores.tprW.rda"), compress=TRUE);
		save(FMM.flat, FMM.tpr, file=paste0(perf.dir, "PCM.", norm.type, ".", flat.file, ".hierScores.tprW.rda"), compress=TRUE);
		save(PRC.flat, PRC.tpr, file=paste0(perf.dir, "PRC.", norm.type, ".", flat.file, ".hierScores.tprW.rda"), compress=TRUE);
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
			dag.dir=dag.dir,flat.norm.dir=flat.norm.dir);
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
	S.tpr <- c();

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
		S.tpr <- rbind(S.tpr, pred.test);
	}
	## put the rows (i.e. genes) of assembled k sub-matrix in the same order of the full beginning matrix
	S.tpr <- S.tpr[rownames(S),];

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
	FMM.tpr <- list(average=FMM.avg, per.example=FMM.per.example);

	## averaging macro-AUC (average and per.class) across k testing sets 
	AUC.average.tpr.over.test <- Reduce("+", AUC.average.test)/kk;
	AUC.class.tpr.over.test <- Reduce("+", AUC.class.test)/kk;
	AUC.tpr <- list(average=AUC.average.tpr.over.test, per.class=AUC.class.tpr.over.test);

	## averaging PxR (average and per.class) across k testing sets 
	PXR.average.tpr.over.test <- Reduce("+", PXR.average.test)/kk;
	PXR.class.tpr.over.test <- Reduce("+", PXR.class.test)/kk;
	PXR.tpr <- list(average=PXR.average.tpr.over.test, per.class=PXR.class.tpr.over.test);

	## averaging PRC (average and per.class) across k testing sets 
	PRC.average.tpr.over.test <- Reduce("+", PRC.average.test)/kk;
	PRC.class.tpr.over.test <- Reduce("+", PRC.class.test)/kk;
	PRC.tpr <- list(average=PRC.average.tpr.over.test, per.class=PRC.class.tpr.over.test);

	## Storing Results #########
	if(norm){
		save(S.tpr, file=paste0(hierScore.dir, flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
		save(AUC.flat, AUC.tpr, file=paste0(perf.dir, "AUC.", flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
		save(PXR.flat, PXR.tpr, file=paste0(perf.dir, "PXR.", flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
		save(FMM.flat, FMM.tpr, file=paste0(perf.dir, "PCM.", flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
		save(PRC.flat, PRC.tpr, file=paste0(perf.dir, "PRC.", flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
	}else{
		save(S.tpr, file=paste0(hierScore.dir, norm.type,".", flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
		save(AUC.flat, AUC.tpr, file=paste0(perf.dir, "AUC.", norm.type, ".", flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
		save(PXR.flat, PXR.tpr, file=paste0(perf.dir, "PXR.", norm.type, ".", flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
		save(FMM.flat, FMM.tpr, file=paste0(perf.dir, "PCM.", norm.type, ".", flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
		save(PRC.flat, PRC.tpr, file=paste0(perf.dir, "PRC.", norm.type, ".", flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
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
			dag.dir=dag.dir,flat.norm.dir=flat.norm.dir);
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
	AUC.tpr <- AUROC.single.over.classes(ann.test, S.test);

	## PxR at fixed recall levels (average and per.class) computed with PerfMeas package on the test set
	PXR.tpr <- precision.at.multiple.recall.level.over.classes(ann.test, S.test);

	## F.measure: Computing Hierarchical Examples-Measures 
	FMM.tpr <- find.best.f(ann.test, S.test, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=TRUE);	

	## Hierarchical PRC (average and per.class) computed by precrec package
	PRC.tpr <- AUPRC.single.over.classes(ann.test, S.test); 

	## storing the hierarchical matrix
	S.tpr <- S.test;
	rm(S.test, S.training);

	## Storing Results #########
	if(norm){
		save(S.tpr, file=paste0(hierScore.dir, flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
		save(AUC.flat, AUC.tpr, file=paste0(perf.dir, "AUC.", flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
		save(PXR.flat, PXR.tpr, file=paste0(perf.dir, "PXR.", flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
		save(FMM.flat, FMM.tpr, file=paste0(perf.dir, "PCM.", flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
		save(PRC.flat, PRC.tpr, file=paste0(perf.dir, "PRC.", flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
	}else{
		save(S.tpr, file=paste0(hierScore.dir, norm.type,".", flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
		save(AUC.flat, AUC.tpr, file=paste0(perf.dir, "AUC.", norm.type, ".", flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
		save(PXR.flat, PXR.tpr, file=paste0(perf.dir, "PXR.", norm.type, ".", flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
		save(FMM.flat, FMM.tpr, file=paste0(perf.dir, "PCM.", norm.type, ".", flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
		save(PRC.flat, PRC.tpr, file=paste0(perf.dir, "PRC.", norm.type, ".", flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
	}
}

#' @title TPR-DAG vanilla
#' @description High level function to compute hierarchical correction according to TPR-DAG vanilla algorithm
#' @param norm boolean value: 
#' \itemize{
#' \item \code{TRUE} (def.): the flat scores matrix has been already normalized in according to a normalization method;	
#' \item \code{FALSE}: the flat scores matrix has not been normalized yet. See the parameter \code{norm.type} for which normalization can be applied.
#' }
#' @param norm.type three values: 
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
#' for each class and example considered. It stored in \code{hierScore.dir} directory.
#' \item \code{PCM} (Protein Centric Measures) average and per-example: compute \code{F-score} measure by \code{find.best.f} function. 
#' It stored in \code{perf.dir} directory.
#' \item \code{PRC} (area under Precision-Recall Curve) average and per.class: compute \code{PRC} by \pkg{precrec} package. 
#' It stored in \code{perf.dir} directory.
#' \item \code{AUC} (Area Under ROC Curve) average and per-class: compute \code{AUC} by \pkg{precrec} package. 
#' It stored in \code{perf.dir} directory.
#' \item \code{PxR} (Precision at fixed Recall levels) average and per classes: compute \code{PxR}  by \pkg{PerfMeas} package. 
#' It stored in \code{perf.dir} directory.
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
			dag.dir=dag.dir,flat.norm.dir=flat.norm.dir);
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
	AUC.tpr <- AUROC.single.over.classes(ann, S);
	
	## Hierarchical PxR at fixed recall levels (average and per.class) computed by PerfMeas package
	PXR.tpr <- precision.at.multiple.recall.level.over.classes(ann, S);

	## Computing Hierarchical Examples-Measures 
	FMM.tpr <- find.best.f(ann, S, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=TRUE);

	## Hierarchical PRC (average and per.class) computed by precrec package
	PRC.tpr <- AUPRC.single.over.classes(ann, S); 

	## storing the hierarchical matrix
	S.tpr <- S;
	rm(S);

	## Storing Results #########
	if(norm){
		save(S.tpr, file=paste0(hierScore.dir, flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
		save(AUC.flat, AUC.tpr, file=paste0(perf.dir, "AUC.", flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
		save(PXR.flat, PXR.tpr, file=paste0(perf.dir, "PXR.", flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
		save(FMM.flat, FMM.tpr, file=paste0(perf.dir, "PCM.", flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
		save(PRC.flat, PRC.tpr, file=paste0(perf.dir, "PRC.", flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
	}else{
		save(S.tpr, file=paste0(hierScore.dir, norm.type,".", flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
		save(AUC.flat, AUC.tpr, file=paste0(perf.dir, "AUC.", norm.type, ".", flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
		save(PXR.flat, PXR.tpr, file=paste0(perf.dir, "PXR.", norm.type, ".", flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
		save(FMM.flat, FMM.tpr, file=paste0(perf.dir, "PCM.", norm.type, ".", flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
		save(PRC.flat, PRC.tpr, file=paste0(perf.dir, "PRC.", norm.type, ".", flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
	}
}

#' @title TPR-DAG vanilla holdout
#' @description High level function to correct the scores with a hierarchy according to TPR-DAG vanill algorithm applying a classical holdout procedure
#' @param norm boolean value: 
#' \itemize{
#' \item \code{TRUE} (def.): the flat scores matrix has been already normalized in according to a normalization method;	
#' \item \code{FALSE}: the flat scores matrix has not been normalized yet. See the parameter \code{norm} for which normalization can be applied.
#' }
#' @param norm.type three values: 
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
#' 	\item \code{hierarchical scores matrix}: a matrix with examples on rows and classes on columns representing the computed hierarchical scores 
#' 	for each class and example considered. It stored in \code{hierScore.dir} directory.
#' 	\item \code{PCM} (Protein Centric Measures) average and per-example: compute \code{F-score} measure by \code{find.best.f} function. 
#' 	It stored in \code{perf.dir} directory.
#' 	\item \code{PRC} (area under Precision-Recall Curve) average and per.class: compute \code{PRC} by \pkg{precrec} package. 
#' 	It stored in \code{perf.dir} directory.
#' 	\item \code{AUC} (Area Under ROC Curve) average and per-class: compute \code{AUC} by \pkg{precrec} package. 
#' 	It stored in \code{perf.dir} directory.
#' 	\item \code{PxR} (Precision at fixed Recall levels) average and per classes: compute \code{PxR}  by \pkg{PerfMeas} package. 
#' 	It stored in \code{perf.dir} directory.
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
			dag.dir=dag.dir,flat.norm.dir=flat.norm.dir);
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
	AUC.tpr <- AUROC.single.over.classes(ann, S);
	
	## Hierarchical PxR at fixed recall levels (average and per.class) computed by PerfMeas package
	PXR.tpr <- precision.at.multiple.recall.level.over.classes(ann, S);

	## Computing Hierarchical Examples-Measures 
	FMM.tpr <- find.best.f(ann, S, n.round=n.round, f.criterion =f.criterion, verbose=FALSE, b.per.example=TRUE);

	## Hierarchical PRC (average and per.class) computed by precrec package
	PRC.tpr <- AUPRC.single.over.classes(ann, S); 

	## storing the hierarchical matrix
	S.tpr <- S;
	rm(S);

	## Storing Results #########
	if(norm){
		save(S.tpr, file=paste0(hierScore.dir, flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
		save(AUC.flat, AUC.tpr, file=paste0(perf.dir, "AUC.", flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
		save(PXR.flat, PXR.tpr, file=paste0(perf.dir, "PXR.", flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
		save(FMM.flat, FMM.tpr, file=paste0(perf.dir, "PCM.", flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
		save(PRC.flat, PRC.tpr, file=paste0(perf.dir, "PRC.", flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
	}else{
		save(S.tpr, file=paste0(hierScore.dir, norm.type,".", flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
		save(AUC.flat, AUC.tpr, file=paste0(perf.dir, "AUC.", norm.type, ".", flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
		save(PXR.flat, PXR.tpr, file=paste0(perf.dir, "PXR.", norm.type, ".", flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
		save(FMM.flat, FMM.tpr, file=paste0(perf.dir, "PCM.", norm.type, ".", flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
		save(PRC.flat, PRC.tpr, file=paste0(perf.dir, "PRC.", norm.type, ".", flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
	}
}


## Attach short description of HEMDAG library when load it
.onAttach <- function(libname=.libPaths(), pkgname="HEMDAG"){
	packageStartupMessage("HEMDAG: Hierarchical Ensemble Methods for DAG-structured taxonomies\n")
}
