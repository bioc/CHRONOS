\name{CHRONOSrun}
\alias{CHRONOSrun}
\title{Default run of CHRONOS}
\usage{
CHRONOSrun(mRNAexp, mRNAlabel, miRNAexp, pathType, subType, measures,
    thresholds, org, export, verbose, miRNAinteractions)
}
\arguments{
\item{mRNAexp}{mRNA expressions filename located in CHRONOS/extdata/Input}

\item{mRNAlabel}{mRNA nomenclature 
(for supported types see \code{\link{convertNomenclature}})}

\item{miRNAexp}{miRNA expressions filename located in CHRONOS/extdata/Input}

\item{pathType}{Pathway type 
('Metabolic', 'Non-Metabolic', 'All' or vector of pathway ids)}

\item{subType}{Subpathway type ('Linear', 'Non-Linear', 'All')}

\item{measures}{Include subpathway structural and functional aspects 
('TRUE', 'FALSE')}

\item{thresholds}{Subscore, mirscore and p-value thresholds \cr
c('pvalue'=pvalue, 'subscore'=subscore, 'mirscore'=mirscore)}

\item{org}{KEGG organism identifier}

\item{export}{Export file type ('.xlsx', '.txt') }

\item{verbose}{Show informative messages (TRUE/FALSE).}

\item{miRNAinteractions}{Edgelist of miRNA-mRNA interactions.}
}
\value{
.
}
\description{
Default run of CHRONOS
}
\details{
\itemize{
\item{ Imports gene and miRNA expressions from 
        CHRONOS/extdata/Input/<mRNAexpFile>.txt \cr
        and CHRONOS/extdata/Input/<miRNAexpFile>.txt}
\item{ Downloads all availiable pathways for the specified organism 
        from KEGG.}
\item{ Creates pathway graphs from downloaded KGML files.}
\item{ Extracts linear subpathways from metabolic and non 
        metabolic graphs.}
\item{ Extracts non linear subpathways from metabolic and non metabolic 
        graphs.}
\item{ Downloads miRecords miRNA-mRNA interactions.}
\item{ Scores and evaluates (linear and non linear) subpathways to extract
        significant results.}
\item{ Organism identifier.}
\item{ Visualizes most the significant results ('.xlsx' or '.txt').}
\item{ Display informative messages (TRUE/FALSE).}
\item{ User-defined miRNA-mRNA interactions can be supplied in the form of 
        an edgelist with two columns. If no such information is availiable, 
        a missing or a NULL argument forces the use of default interactions 
        by using \code{\link{downloadMiRecords}}. }
}
}
\examples{

# Default run

load(system.file('extdata', 'Examples//data.RData', package='CHRONOS'))

res <- CHRONOSrun(  mRNAexp=mRNAexpr,
                    mRNAlabel='entrezgene',
                    miRNAexp=miRNAexpr,
                    pathType=c('04915', '04917', '04930', '05031'),
                    org='hsa',
                    subType='Linear',
                    thresholds=c('subScore'=0.4, 'mirScore'=0.4),
                    miRNAinteractions=miRNAinteractions)
}

