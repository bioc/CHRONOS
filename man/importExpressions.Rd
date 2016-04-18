\name{importExpressions}
\alias{importExpressions}
\title{Import gene and miRNA expressions from}
\usage{
importExpressions(data, type, sep, org, mRNAnomenclature)
}
\arguments{
\item{data}{Expressions data filename or matrix.}

\item{type}{Expressions data type. (or mRNA expressions, type=<nomenType>.
Availiable gene expression nomenclature can be found in
\code{\link{convertNomenclature}}.
For miRNA expressions, type='miRNA'.}

\item{sep}{File delimiter.}

\item{org}{KEGG organism identifier}

\item{mRNAnomenclature}{Nomenclature of user's mRNA expressions}
}
\value{
.
}
\description{
Import gene and miRNA expressions from
}
\details{
\itemize{
\item{ Import gene expressions data from CHRONOS/extdata/Input/<userFile>.txt
        or a supplied matrix. }
\item{ Import miRNA expressions data from CHRONOS/extdata/Input/<userFile>.txt
        or a supplied matrix.}
}
}
\examples{

# Example

load(system.file('extdata', 'Examples//data.RData', package='CHRONOS'))

importExpressions(data=mRNAexpr, type='mRNA', 
                mRNAnomenclature='entrezgene', sep='\\t', org='hsa')
importExpressions(data=miRNAexpr, type='miRNA', sep='\\t', org='hsa')

}

