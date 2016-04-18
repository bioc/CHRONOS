\name{extractNonLinearSubpathways}
\alias{extractNonLinearSubpathways}
\title{Non linear subpathway extraction from pathway graphs}
\usage{
extractNonLinearSubpathways(graphs, pathways, a, b, k, filter, groupMode,
    export, verbose)
}
\arguments{
\item{graphs}{Pathway graphs as returned from 
\code{\link{createPathwayGraphs}}.}

\item{pathways}{The subset of pathways from whom subpathways are to be 
extracted. \cr
If missing, all pathway graphs are used.}

\item{a}{Minimum subpathway length.}

\item{b}{Maximum subpathway length.}

\item{k}{Clique size.}

\item{filter}{Filter the subpaths with user genes (TRUE). \cr}

\item{groupMode}{Expand paralogues ('expand') or collapse them to a single 
entry ('collapse'). \cr}

\item{export}{Exports subpaths in CHRONOS/extdata/Output/Subpaths/Non-Linear/
<org> folder. \cr
Availiable formats are '.txt' and/or '.RData'.}

\item{verbose}{Display informative messages (TRUE)

Requires previous execution of \code{\link{importExpressions}}.}
}
\value{
Returns a list consisting of
\itemize{
    \item{A matrix of linear subpathways  }{(subpaths)}
    \item{A list of processed pathway graphs adjacency matrices}{(adjMats)}
    \item{A list of processed pathway genes and interactions between them}
    {(lexicon)}
}

To extract non linear subpaths from a pathway graph, all interactions 
between nodes of belonging to k-cliques are found. The ones that correspond


To extract non linear subpaths from a pathway graph, all interactions 
between nodes of belonging to k-cliques are found. The ones that correspond
to actual interactions between genes make up the non linear subpath.
}
\description{
Non linear subpathway extraction from pathway graphs
}
\examples{

# Load pathway graphs from toy data
load(system.file('extdata', 'Examples//data.RData', package='CHRONOS'))

# Extract linear subpathways
nliSubs  <- extractNonLinearSubpathways(graphs=graphs)

}

