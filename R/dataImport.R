
importExpressions      <- function( data, type, sep, org, mRNAnomenclature )
{
    #
    # Import gene or mirna expression data.
    # Data can be imported either by providing the filename of a txt file
    # containing them or by directly supplying a matrix.
    # 
    message('Importing ', type , ' data...', appendLF = FALSE)

    if (missing(type))         
    { 
        message('\nError. Please specify a type.')
        return()
    }
    if (missing(sep))               { sep  <- '\t' }

    if (missing(mRNAnomenclature))  { mRNAnomenclature  <- 'entrezgene' }
    
    if (!'miRNA' %in% type)  
    {
        mRNAExpr <- importGeneExpressions(org=org, sep=sep, file=data, 
                    mRNAnomenclature=mRNAnomenclature) 

        if (org == 'hsa')
        {
            cache$nomenclature <- c('from'=mRNAnomenclature, 
                                    'to'='hgnc_symbol')
        }
        if (org == 'mmu')
        {
            cache$nomenclature <- c('from'=mRNAnomenclature, 
                                    'to'='mgi_symbol')
        }

        message('done.')

        return(mRNAExpr)
    }
    if ('miRNA' %in% type) 
    {
        miRNAExpr <- importMiRNAExpressions(org=org, sep=sep, file=data)

        message('done.')

        return(miRNAExpr)
    }
}


importGeneExpressions  <- function( org, sep, file, mRNAnomenclature )
{
    # File can either be a filepath in input directory or a matrix
    if (class(file) == 'character')
    {
        # If a different filename than the default has been supplied, 
        # replace it at the filepath
        dir <- unlist(strsplit(cache$dirs$geneExpressions, '//'))
        dir <- dir[-length(dir)]
        dir <- paste(dir, collapse='//')
        cache$dirs$geneExpressions <- paste(dir, file, sep='//')

        if (!file.exists(cache$dirs$geneExpressions))
        { 
            message('No gene expression data availiable.')
            return(invisible()) 
        }

        # Get file extension
        ext <- tail(unlist(strsplit(file, '\\.')), 1)

        if (tolower(ext) == 'rdata')
        {
            e=new.env(); load(file=cache$dirs$geneExpressions, envir=e)

            if (length(ls(e)) == 0) 
            { 
                message('File ', file, ' is empty.\n') 
            }

            geneEx <- e[[ls(e)[1]]]

            if (class(geneEx) != 'matrix') 
            { 
                message('Please supply a matrix of gene expressions.')
                message('File ', file, ' does not contain a matrix.\n') 
            }
            if (is.null(rownames(geneEx)))
            {
                message('mRNA expression matrix must contain gene ids.')
            }
            if (ncol(geneEx) == 1)
            {
                message('mRNA expression matrix should contain more 
                    than one columns.')
            }

            return(geneEx) 
        }

        filepath <- cache$dirs$geneExpressions

        # Read expressions file
        data    <- as.matrix(read.table(filepath, sep = sep, as.is=TRUE))
        ids     <- data[,1]
        geneEx  <- matrix(as.numeric(data[,-1]), nrow=nrow(data))
        rownames(geneEx) <- ids
        
        # Remove NA's
        geneEx  <- na.omit(geneEx)
        attributes(geneEx)$'na.action' <- NULL

        # Replace txt file with processed RData file and save data
        file <- gsub('.txt', '.RData', cache$dirs$geneExpressions)
        cache$dirs$geneExpressions <- file 
    }

    if (class(file) == 'matrix')
    {
        geneEx  <- file
        dir     <- unlist(strsplit(cache$dirs$geneExpressions, '//'))
        dir     <- dir[-length(dir)]
        dir     <- paste(dir, collapse='//')
        file    <- paste(dir, 'geneExpressions.RData', sep='//')
        cache$dirs$geneExpressions <- file
    }

    if (mRNAnomenclature != 'entrezgene')
    {
        #
        # Convert nomenclature to entrez if necessary
        #
        
        # Check if supplied type is supported
        supportedNomenclature <- names(.getFilters())

        if (!(mRNAnomenclature %in% supportedNomenclature)) 
        {
            strErr <- sprintf('Error, %s nomenclature not supported\n', 
                                mRNAnomenclature)
            stop(strErr)
            return(NULL)
        }

        ids <- rownames(geneEx)
        res <- convertNomenclature(ids=ids, org=org, from=mRNAnomenclature, 
                                    to='entrezgene')

        # Remove conversion giving the same answer
        idx <- which((res[, 1] == res[, 2]))

        if (length(idx) > 0)
        {
            res <- res[-idx, ]    
        }

        # Create a hash to map distinct refseq ids to entrez ids.
        lib1        <- res[,2]
        names(lib1) <- res[,1]
        renamed     <- lib1[ids]
        
        # Remove gene expressions with ids with no entrez number
        idx         <- which(is.na(lib1[ids]))

        if (length(idx) > 0)
        {
            geneEx      <- geneEx[-idx, ]
            ids         <- ids[-idx]
            renamed     <- renamed[-idx]
        }
        rownames(geneEx) <- renamed

        # Sort expression matrix according to name
        geneEx <- geneEx[order(as.numeric(rownames(geneEx)) ),]
    }

    # Remove multiple entries
    geneEx <- geneEx[which(!duplicated(rownames(geneEx))), , drop=FALSE]
    
    save(geneEx, file=cache$dirs$geneExpressions)

    return(geneEx)
}


importMiRNAExpressions <- function(org, sep, file)
{
    # File can either be a filepath in input directory or a matrix
    if (class(file) == 'character')
    {
        # If a different filename than the default has been supplied, 
        # replace it at the filepath
        dir <- unlist(strsplit(cache$dirs$mirnaExpressions, '//'))
        dir <- dir[-length(dir)]
        dir <- paste(dir, collapse='//')
        cache$dirs$mirnaExpressions <- paste(dir, file, sep='//')

        if (!file.exists(cache$dirs$mirnaExpressions))
        { 
            message('No miRNA expression data availiable.')
            return(invisible()) 
        }

        # Get file extension
        ext <- tail(unlist(strsplit(file, '\\.')), 1)

        if (tolower(ext) == 'rdata')
        { 
            e=new.env(); load(file=cache$dirs$mirnaExpressions, envir=e)

            if (length(ls(e)) == 0) 
            { 
                message('File ', file, ' is empty.\n') 
            }

            miEx <- e[[ls(e)[1]]]

            if (class(miEx) != 'matrix') 
            { 
                message('Please supply a matrix of miRNA expressions.')
                message('File ', file, ' does not contain a matrix.\n') 
            }
            if (is.null(rownames(miEx)))
            {
                message('miRNA expression matrix must contain miRNA ids.')
            }
            if (ncol(miEx) == 1)
            {
                message('miRNA expression matrix should contain more 
                    than one columns.')            }

            return(miEx) 
        }

        filepath <- cache$dirs$mirnaExpressions

        # Read expressions file
        data    <- as.matrix(read.table(filepath, sep = sep, as.is=TRUE))

        ids     <- data[,1]
        miEx  <- matrix(as.numeric(data[,-1]), nrow=nrow(data))
        rownames(miEx) <- ids
        
        # Remove NA's
        miEx  <- na.omit(miEx)
        attributes(miEx)$'na.action' <- NULL

        # Replace txt file with processed RData file and save data
        cache$dirs$mirnaExpressions <- gsub('.txt', '.RData', 
                                    cache$dirs$mirnaExpressions) 

    }

    if (class(file) == 'matrix')
    {
        miEx    <- file
        dir     <- unlist(strsplit(cache$dirs$mirnaExpressions, '//'))
        dir     <- dir[-length(dir)]
        dir     <- paste(dir, collapse='//')
        cache$dirs$mirnaExpressions <- paste0(dir, '//miRNAExpressions.RData')
    }

    result <- tryCatch(
    {
        # Conform user miRNA annotations to the ones used by miRecords
        userMirs       <- rownames(miEx)
        finalMirs      <- convertMiRNANomenclature(org, userMirs)
        rownames(miEx) <- finalMirs[,2]
    }, warning = function(war) 
    {
    }, error = function(e) { }
    , finally = {}) 
    
    save(miEx, file=cache$dirs$mirnaExpressions)

    return(miEx)
}


importMiRNAFile        <- function(org)
{
    file <- paste(cache$dirs$miDir, org, cache$dirs$miFile, sep='//')

    e=new.env(); load(file=file, envir=e)
    df <- as.data.frame(e$targets, stringsAsFactors=FALSE)

    return(df)
}

saveMiRNAFile        <- function(org, miRNAinteractions)
{
    # Create organism subdirectory if it does not exist
    orgDir <- paste0(cache$dirs$miDir, '//' ,org) 
    if ( !file.exists(orgDir) )
    {
        dir.create(orgDir, showWarnings=FALSE, recursive=TRUE)
    }

    file <- paste0(orgDir, '//', cache$dirs$miFile)

    targets <- as.data.frame(miRNAinteractions, stringsAsFactors=FALSE)
    save(file=file, targets, compress='xz')

    return(targets)
}
