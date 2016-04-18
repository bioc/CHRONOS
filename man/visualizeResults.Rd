\name{visualizeResults}
\alias{visualizeResults}
\title{Visualize results in tabular form (txt, xlsx)}
\usage{
visualizeResults(summary, export, expand, colors, from, to)
}
\arguments{
\item{summary}{Evaluation results as returned from 
\code{\link{scoreSubpathways}}}

\item{export}{'.xlsx' exports a xlsx file and '.txt' a .txt file.}

\item{expand}{TRUE if each subpathway member and miRNA belongs to a single 
cell, FALSE if all subpathway members belong to one cell and
miRNAs to another cell.}

\item{colors}{The color scheme used in subScores heatmap.}

\item{from}{Primary annotation \code{\link{convertNomenclature}}.
Defaults to EntrezGene ID.}

\item{to}{Secondary annotation \code{\link{convertNomenclature}}}
}
\value{
A txt or a xlsx file in CHRONOS/extdata/Output/Scores/Linear/<org> \cr
or CHRONOS/extdata/Output/Scores/Non-Linear/<org>
}
\description{
Visualize results in tabular form (txt, xlsx)
}


\examples{

# Load scored subpathways from toy data
load(system.file('extdata', 'Examples//data.RData', package='CHRONOS'))

visualizeResults(linSubsScored, export='txt')

}

