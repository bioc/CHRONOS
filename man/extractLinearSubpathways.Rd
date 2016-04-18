\name{extractLinearSubpathways}
\alias{extractLinearSubpathways}
\title{Linear subpathway extraction from pathway graphs}
\usage{
extractLinearSubpathways(graphs, pathways, a, b, filter, export, groupMode,
    verbose)
}
\arguments{
\item{graphs}{Pathway graphs as returned from 
\code{\link{createPathwayGraphs}}.}

\item{pathways}{The subset of pathways from whom subpathways are to be 
extracted. \cr
If missing, all pathway graphs are used.}

\item{a}{Minimum subpathway length.}

\item{b}{Maximum subpathway length.}

\item{filter}{Filter the subpaths with user genes (TRUE). \cr}

\item{export}{Exports subpaths in CHRONOS/extdata/Output/Subpaths/Linear/<org> 
folder. \cr 
Availiable formats are '.txt' and/or '.RData'.}

\item{groupMode}{Expand paralogues ('expand') or collapse them to a 
single entry ('collapse'). \cr}

\item{verbose}{Display informative messages (TRUE)

Requires previous execution of \code{\link{importExpressions}}.}
}
\value{
Returns a list consisting of
\itemize{
    \item{A matrix of linear subpathways  }{(subpaths)}
    \item{A list of processed pathway graphs adjacency matrices  }{(adjMats)}
    \item{A list of processed pathway genes and interactions between them  }
    {(lexicon)}
}
}
\description{
Linear subpathway extraction from pathway graphs
}
\details{
Subpath filtering supports the removal of subpaths that have at least 
one member not belonging to the set of user supplied genes. These genes 
are extracted from the user's mRNA expressions matrix.Thus, the execution  
of \code{\link{importExpressions}} is a prerequisite.

To extract linear subpathways from a pathway graph, all possible start 
and end nodes are considered. A start node has only outgoing edges while
an end node only has incoming edges. For each such pair, all linear  
subpathways are found by traversing the corresponding graph. Since the 
initial pathway graph's nodes contain one or more genes, resulting 
subpathways consist of bins of one or more genes. These subpaths are 
expanded to subpathways with one gene per bin in order to obtain usable 
subpathways.
}
\examples{

# Load pathway graphs from toy data
load(system.file('extdata', 'Examples//data.RData', package='CHRONOS'))

# Extract linear subpathways
linSubs <- extractLinearSubpathways(graphs=graphs)

}
