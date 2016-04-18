
convertNomenclature  <- function(ids, org, from, to)
{
    
    tryCatch(
        .mRNAconversion(values=ids, 
                        org=org, 
                        filter=from, 
                        attributes=to),
        error=function(e) ' error ')
}

.mRNAconversion <- function(values, org, filter, attributes)
{
    # Find dataset
    datasets <- listDatasets(useMart("ensembl"))[, 'dataset']
    names(datasets) <- sapply(as.list(datasets), function(x) {substr(x,1,3)})
    dataset  <- unname(datasets[org])

    # Basic error checking
    if (org=='mmu') 
        { dataset <- 'mmusculus_gene_ensembl' }
    if (is.na(dataset) && !missing(org))     
        { cat('No dataset availiable for',org,'.\n'); return() }
    if (is.na(dataset) && missing(org))
        { cat('Please specify an organism.\n'); return() }
    if(missing(filter))
        { cat('Please specify a valid filter\n'); return() }
    if(missing(attributes))
        { cat('Please specify valid attributes.\n'); return() }
    
    # Perform query
    ensembl <- useMart("ensembl", dataset=as.character(dataset))
    attributes <- c(filter, attributes)
    res <- getBM(attributes=attributes, filters=filter, 
                    values=values, mart=ensembl)
    res <- trimws( as.matrix(res) )

    # Remove genes with no conversion
    idx <- which(res[, 2] == '')
    if ( length(idx) > 0 ) { res<- res[-idx, ] }

    rownames(res) <- 1:nrow(res)
    colnames(res) <- .getFilters()[colnames(res)]

    return(res)
}

.getFilters          <- function()
{
    filters     <- vector(mode="numeric", length=14)
    fnames      <- vector(mode="numeric", length=14)
    filters[1]  <- 'Ensembl Gene ID(s)'
    filters[2]  <- 'Ensemble transcript ID'
    filters[3]  <- 'Ensemble Protein ID'
    filters[4]  <- 'HGNC ID'
    filters[5]  <- 'HGNC symbol'
    filters[6]  <- 'HGNC transcript name'
    filters[7]  <- 'EntrezGene ID'
    filters[8]  <- 'Refseq mRNA ID(s)'
    filters[9]  <- 'Refseq protein ID(s)'
    filters[10] <- 'UniProt/Swissprot Accession(s)'
    filters[11] <- 'UniProt/Swissprot ID(s)'
    filters[12] <- 'UniGene ID(s)'
    filters[13] <- 'UniProt Genename ID(s)'
    filters[14] <- 'MGI symbol'
    fnames[1]   <- 'ensembl_gene_id'           
    fnames[2]   <- 'ensembl_transcript_id'
    fnames[3]   <- 'ensembl_peptide_id'
    fnames[4]   <- 'hgnc_id'
    fnames[5]   <- 'hgnc_symbol'
    fnames[6]   <- 'hgnc_transcript_name'
    fnames[7]   <- 'entrezgene'
    fnames[8]   <- 'refseq_mrna'
    fnames[9]   <- 'refseq_peptide'
    fnames[10]  <- 'uniprot_swissprot_accession'
    fnames[11]  <- 'uniprot_swissprot'
    fnames[12]  <- 'unigene'
    fnames[13]  <- 'uniprot_genename'
    fnames[14]  <- 'mgi_symbol'
    names(filters) <- fnames

    return(filters)
}
