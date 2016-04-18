\name{subpathwayKEGGmap}
\alias{subpathwayKEGGmap}
\title{Create links to KEGG pathway map with highlighted subpathways.}
\usage{
subpathwayKEGGmap(subpathways, type, openInBrowser)
}
\arguments{
\item{subpathways}{Subpathways as returned by 
\code{\link{extractLinearSubpathways}} \cr
or \code{\link{extractNonLinearSubpathways}}}

\item{type}{Subpathway type (Linear, Non-Linear)}

\item{openInBrowser}{Open link in default browser.}
}
\value{
Vector of links of KEGG pathway maps.
}
\description{
Create links to KEGG pathway map with highlighted subpathways.
}

\examples{

# Load extracted linear subpathways from toy data
load(system.file('extdata', 'Examples//data.RData', package='CHRONOS'))

# Opening selected subpathways in default browser
subs <- linSubs$subpaths[1:3, ]

subpathwayKEGGmap(subpathways=subs, type='Linear', openInBrowser=FALSE)
}
