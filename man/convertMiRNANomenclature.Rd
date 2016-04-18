\name{convertMiRNANomenclature}
\alias{convertMiRNANomenclature}
\title{Conform miRNA annotations to the ones currently used by miRecords.}
\usage{
convertMiRNANomenclature(org, miRNAs, update)
}
\arguments{
\item{org}{KEGG organism identifier.}

\item{miRNAs}{Vector of miRNAs identifiers.}

\item{update}{Update annotation mapper with latest annotation changes.}
}
\value{
.
}
\description{
Conform miRNA annotations to the ones currently used by miRecords.
}
\details{
Determine which miRNAs are incompatible with miRecords annotations
and retreive the suitable ones from www.mirbase.org. 

}
\examples{

data <- c('hsa-let-7g-5p', 'hsa-miR-154-5p', 'hsa-miR-376b-3p')

convertMiRNANomenclature(org='hsa', miRNAs=data)

}

