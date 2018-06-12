\name{HEMDAG-package}
\alias{HEMDAG-package}
\alias{HEMDAG}
\docType{package}
\title{HEMDAG: Hierarchical Ensemble Methods for Directed Acyclic Graphs}

\description{
Implementation of Hierarchical Ensemble Methods for DAGs: HTD-DAG (Hierarchical Top Down) and TPR-DAG (True Path Rule).
}

\details{
HEMDAG can be used to enhance the predictions of virtually any off-the-shelf flat learning method, merely by taking into account the hierarchical nature 
of the classes of a bio-ontology. HEMDAG is specifically designed for exploiting the hierarchical relationships of DAG-structured taxonomies, 
such as the Human Phenotype Ontology (HPO) or the Gene Ontology (GO), but it can be also safely applied to tree-structured taxonomies (as FunCat), 
since obviously trees are DAGs. Finally HEMDAG scale nicely both in terms of the complexity of the taxonomy and in the cardinality of the examples.
}

\author{
\emph{Marco Notaro} and \emph{Giorgio Valentini}

\href{https://sites.google.com/site/anacletolaboratory/}{AnacletoLab} 

DI, Dipartimento di Informatica

Universita' degli Studi di Milano

Maintainer: 
\emph{Marco Notaro} 
\email{marco.notaro@unimi.it}
}

\references{
Marco Notaro, Max Schubach, Peter N. Robinson and Giorgio Valentini, 
\emph{Prediction of Human Phenotype Ontology terms by means of Hierarchical Ensemble methods},
BMC Bioinformatics 2017
}

\keyword{package}
