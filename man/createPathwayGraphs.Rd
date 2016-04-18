\name{createPathwayGraphs}
\alias{createPathwayGraphs}
\title{Convert KEGG Pathways to Gene-Gene Network Graphs.}
\usage{
createPathwayGraphs(org, pathways, edgeTypes, doubleEdges, choice, groupMode)
}
\arguments{
\item{org}{KEGG organism identifier.}

\item{pathways}{Vector of KEGG pathway identifiers.}

\item{edgeTypes}{Vector of edge types mappings.}

\item{doubleEdges}{Specify which edgeTypes should be considered 
bidirectional.}

\item{choice}{Create metabolic graph either by using relations or reactions
from KGML file ('reactions', 'relations')\cr}

\item{groupMode}{'expand' to consider each group member a node, or 
'collapse' to consider all components' genes as a node}
}
\value{
A list containing a list of compact adjacency matrices, a list of 
expanded adjacency matrices, and list detailing all nodes, edges and 
interaction types.
}
\description{
Convert KEGG Pathways to Gene-Gene Network Graphs.
}
\details{
KEGG pathways consist of nodes each one containing one or more genes. 
Thus, two kinds of adjacency matrices are created. The compact adjacency 
matrix retains the groupings and stores edge types between genes and genes,
genes and groups of genes or between group of genes. The expanded adjacency
matrix stores edge type information between individual genes. 


}
\examples{

# Download Insulin Signaling Pathway
pathways <- c('04915', '04917', '04930', '05031')
paths    <- downloadPathways(org='hsa', pathways=pathways)

# Create pathway graph
graphs   <- createPathwayGraphs(org='hsa', pathways=paths)


}
\references{
Li, C., Han, J., Yao, Q., Zou, C., Xu, Y., Zhang, C., ... & Li, X. (2013).
Subpathway-GM: identification of metabolic subpathways via joint power 
of interesting genes and metabolites and their topologies within pathways.
Nucleic acids research, 41(9), e101-e101.
}

