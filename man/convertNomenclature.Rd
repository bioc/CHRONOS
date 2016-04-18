\name{convertNomenclature}
\alias{convertNomenclature}
\title{Convert genes identifier nomenclature.}
\usage{
convertNomenclature(ids, org, from, to)
}
\arguments{
\item{ids}{Vector of gene identifiers}

\item{org}{KEGG organism identifier}

\item{from}{Initial identifier type}

\item{to}{A vector of final identifier types}
}
\value{
Vector of converted gene identifiers
}
\description{
Convert genes identifier nomenclature.
}
\details{
\tabular{ll}{
    EntrezGene ID \tab 'entrezgene'\cr
    Ensembl Gene ID \tab 'ensembl_gene_id'\cr
    Ensemble Transcript ID \tab'ensembl_transcript_id'\cr
    Ensemble Protein ID \tab 'ensembl_peptide_id'\cr
    HGNC ID \tab 'hgnc_id'\cr
    HGNC Symbol \tab 'hgnc_symbol'\cr
    HGNC Transcript name \tab 'hgnc_transcript_name'\cr
    Refseq mRNA ID \tab 'refseq_mrna'\cr
    Refseq Protein ID \tab 'refseq_peptide'\cr
    UniProt/Swissprot Accession \tab 'uniprot_swissprot_accession'\cr
    UniProt/Swissprot ID \tab 'uniprot_swissprot'\cr
    UniGene ID \tab 'unigene'\cr
    UniProt Genename ID \tab 'uniprot_genename'
}
}
\examples{
# Identifiers to be converted
ids <- c('5091', '5105')

# Convert to HGNC ID, Ensembl Gene ID and UniProt Genename ID
from <- 'entrezgene'
to <- c('hgnc_symbol','ensembl_gene_id', 'uniprot_genename')
\dontrun{res <- convertNomenclature(ids=ids, org='hsa', from=from, to=to)}

}

