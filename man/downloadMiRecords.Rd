\name{downloadMiRecords}
\alias{downloadMiRecords}
\title{Download miRNA-mRNA interactions for an organism.}
\usage{
downloadMiRecords(org, pn, update, databases)
}
\arguments{
\item{org}{KEGG organism identifier.}

\item{pn}{Number of databases that verify miRNA-mRNA interactions.}

\item{update}{Download preprocessed data (update=FALSE) or new data from 
miRecords (update=TRUE).}

\item{databases}{Specify which miRNA-mRNA interaction databases 
will be used.}
}
\value{
Downloaded data is stored in 
    CHRONOS/extdata/Downloads/miRecords/<org>/miRNATargets.RData
}
\description{
Download miRNA-mRNA interactions for an organism.
}
\details{
miRecords is a resource for animal miRNA-target interactions. 
The Predicted Targets component of miRecords 
is an integration of predicted miRNA targets produced by 11 established 
miRNA target prediction tools, namely DIANA-microT, MicroInspector, 
miRanda, MirTarget2, miTarget, NBmiRTar, PicTar, PITA, RNA22, RNAhybrid, 
and TargetScan/TargertScanS.
}
\examples{

# Load extracted linear subpathways from toy data
load(system.file('extdata', 'Examples//data.RData', package='CHRONOS'))

\dontrun{downloadMiRecords(org='hsa', pn=5, update=FALSE, databases='All')}

}
\references{
\itemize{
    \item{\url{http://c1.accurascience.com/miRecords}}
}
}

