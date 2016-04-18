\name{pathwayMeasures}
\alias{pathwayMeasures}
\title{Pathway structural and functional aspects}
\usage{
pathwayMeasures(graphs)
}
\arguments{
\item{graphs}{Pathway graphs as returned from
\code{\link{createPathwayGraphs}}.}
}
\value{
Matrix with pathness, betweeness centrality and degree values 
for each gene in the pathway graphs at it's columns.

}
\description{
Pathway structural and functional aspects
}
\details{
Structural and functional aspects of a pathway are calculated in 
respect to all organism pathways.
}
\examples{

# Load pathway graphs from toy data
load(system.file('extdata', 'Examples//data.RData', package='CHRONOS'))

# Calculate pathway structural and functional aspects
measures <- pathwayMeasures(graphs)

}

