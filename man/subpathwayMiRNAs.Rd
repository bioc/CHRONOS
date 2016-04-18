\name{subpathwayMiRNAs}
\alias{subpathwayMiRNAs}
\title{Create a circulat plot of a subpathway and the miRNAs that target it.}
\usage{
subpathwayMiRNAs(summary, subIdx, timePoints)
}
\arguments{
\item{summary}{Output from \code{\link{scoreSubpathways}}}

\item{subIdx}{Subpathway index}

\item{timePoints}{Time points to include in visualization, default to all.}
}
\value{
.
}
\description{
Create a circulat plot of a subpathway and the miRNAs that target it.
}


\examples{

# Load scored subpathways from toy data
load(system.file('extdata', 'Examples//data.RData', package='CHRONOS'))
# Visualize one or more subpathways. 
subpathwayMiRNAs(summary=linSubsScored, subIdx=2)

}

