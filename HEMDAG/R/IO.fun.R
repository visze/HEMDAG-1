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
