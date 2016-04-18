\name{scoreSubpathways}
\alias{scoreSubpathways}
\title{Evaluate subpathways using an interacting scorng scheme (IS) 
for each time point.}
\usage{
scoreSubpathways(subpathways, filters, measures, parameters, miRNAinteractions)
}
\arguments{
\item{subpathways}{Subpaths as returned from 
\code{\link{extractLinearSubpathways}} \cr
and \code{\link{extractNonLinearSubpathways}}.}

\item{filters}{Named vector of filters used for subpathway evaluation. \cr
Values denote corresponding thresholds.
\tabular{llll}{
    \tab \tab pvalue \tab  Statistical evaluation \cr
    \tab \tab measures \tab Structural and functional evaluation\cr
    \tab \tab subScore \tab mRNA-mRNA interaction scoring\cr
    \tab \tab mirScore \tab miRNA-mRNA interaction scoring
}}

\item{measures}{Subpathway structural and functional aspects as returned 
    from \cr
\code{\link{pathwayMeasures}}.}

\item{parameters}{C,K,T parameters of scoring scheme.}
\item{miRNAinteractions}{An edgelist of miRNA-mRNA interactions used to  
    override downloaded interactions from miRecords.}
}
\value{
\tabular{llll}{
    \tab subpathways \tab  High ranking subpathways \cr
    \tab subScores \tab miRNA-subpathway scores\cr
    \tab mRNAScores \tab mRNA-mRNA scores for each subpathway and for each 
    time point\cr
    \tab miRNAsOverSubpathway \tab High ranking miRNAs hitting 
    each subpathway \cr
    \tab pValues\tab P-value of each subpathway\cr
    \tab filters\tab Filters used for the evaluation
}
}
\description{
Evaluate subpathways using an interacting scorng scheme (IS) 
for each time point.
}
\details{
...
}
\examples{

# Load extracted subpathways from toy data
load(system.file('extdata', 'Examples//data.RData', package='CHRONOS'))

# Import mRNA expressions
mRNAexpr <- importExpressions(data=mRNAexpr, type='mRNA', org='hsa')

# Score extracted linear subpathways
filters       <- c('subScore'=0.4)
linSubsScored <- scoreSubpathways(subpathways=linSubs, filters=filters)

}
\references{
Jethava, V., Bhattacharyya, C., Dubhashi, D., & Vemuri, G. N. (2011).  
Netgem: Network embedded temporal generative model for gene expression  
data. BMC bioinformatics, 12(1), 327.

Kim,Y. et al. (2011). Principal network analysis: identification of 
subnetworks representing major dynamics using gene expression data. 
Bioinformatics, 27(3), 391-398
}

