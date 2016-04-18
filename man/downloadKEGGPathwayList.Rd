\name{downloadKEGGPathwayList}
\alias{downloadKEGGPathwayList}
\title{Retrieve all availiable pathways for an organism.}
\usage{
downloadKEGGPathwayList(org)
}
\arguments{
\item{org}{KEGG organism identifier.}
}
\value{
Data frame of pathway ids and names.
}
\description{
Retrieve all availiable pathways for an organism.
}
\details{
.
}
\examples{

# Load extracted linear subpathways from toy data
load(system.file('extdata', 'Examples//data.RData', package='CHRONOS'))

# Retrieve all availiable hsa pathways
\dontrun{pathways <- downloadKEGGPathwayList(org='hsa')}

}
\references{
\itemize{
    \item{\url{http://www.genome.jp/kegg/pathway.html}}
}
}

