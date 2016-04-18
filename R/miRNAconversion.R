

convertMiRNANomenclature    <- function(org, miRNAs, update)
{
    if (missing(update)) { update <- FALSE }

    # Load or download 
    #  - Mapper from old miRNA annotations to new ones 
    #  - Availiable miRNA annotaions from miRecords

    annFile <- paste0(cache$dirs$anDir, '//miRNAmaps_', org, '.RData')
    if ( update || (!update && !file.exists(annFile)) )  
    {
        # Create mapper
        miRNA.conversions <- mirParser(org)
        # Download miRecords miRNA list
        miRNA.miRecords <- .getMiRecordsMiRNAs(org)
        # Save locally
        save(miRNA.conversions, miRNA.miRecords, file=annFile)        
    }
    if ( !update && file.exists(annFile) ) 
    {
        # If a new map has been downloaded, use it.
        e=new.env(); load(file=annFile, envir=e)
        miRNA.conversions <- e$miRNA.conversions
        miRNA.miRecords <- e$miRNA.miRecords
    }

    # Convert according to user miRNAs
    res <- mirRenamer(  org=org, map=miRNA.conversions, userMirs=miRNAs, 
                        miRecordsMirs=miRNA.miRecords)

    return ( res )
}


mirParser <- function( org )
{
    # Download data from mirbase
    mirbaseURL <- 'ftp://mirbase.org/pub/mirbase/CURRENT/database_files/'
    mirbaseURL <- paste0(mirbaseURL, 'mirna_mature.txt.gz')
    file <- tempfile()
    tryCatch( download.file(mirbaseURL, file, quiet=TRUE) , 
            error=function(e) ' ')
    data <- readLines(file)
    tryCatch( unlink(file), error=function(e) ' ')

    # Process data
    data <- strsplit(data, '\t')
    map <- do.call('rbind', lapply(data, function(x) { x[2:3]} ))
    colnames(map) <- c('Mature', 'Previous')
    # Find organism specific entries
    idx <- which( map[, 1] != gsub(org, '', map[, 1]) )
    if (length(idx) > 0 )
    {
        map <- map[idx, , drop=FALSE]
    }
    # Replace '' with NAs
    map[map[, 2] == '', 2] <- NA

    # Expand multiple old annotations delimited with ';'
    R   <- strsplit(map[,2], ';')
    len <- length(unlist(R))
    res <- matrix(, nrow=len, ncol=ncol(map))
    ctr <- 1
    for (i in 1:nrow(map))
    {
        lR <- length(R[[i]])
        res[ctr:(ctr+lR-1), ] <- as.matrix(expand.grid(map[i,1], R[[i]]))
        ctr <- ctr + lR
    }
    res <- res[-seq(ctr,nrow(res)),]
    colnames(res) <- c('Mature', 'Previous')

    # Remove entries where both old and new annotation are the same
    idx <- which(res[, 1] == res[, 2])
    if ( length(idx) > 0 ) { res <- res[-idx, , drop=FALSE] }
    # Remove duplicacted entries
    mapper <- unique(res)

    return( mapper )
}

mirRenamer                  <- function(map, userMirs, miRecordsMirs, org)
{
    #
    # Converts miRNA annotation to a miRecords-compatible one.
    # miRecords could have either a new or an old annotation.
    # Some other miRNAs might not be availiable from miRecords.
    # 

    # Cleanup data
    map      <- map[which(!is.na(map[,2])),]
    userMirs <- userMirs[which(!is.na(userMirs))]

    # Find which user miRNAs are incompatible with miRecords annotations
    oddMirs    <- unique(userMirs[which(!userMirs %in% miRecordsMirs)])

    # Find a compatible annotation using the old to new miRNA annotation map.

    # New annotations
    newOddMirs <- oddMirs[which(oddMirs %in% map[,1])]

    # Old annotations
    oldOddMirs <- oddMirs[which(oddMirs %in% map[,2])]

    # Not availiable
    naMirs     <- setdiff(oddMirs, c(newOddMirs, oldOddMirs))

    finalMirs  <- userMirs

    # Remove na annotations
    finalMirs[which(userMirs %in% naMirs)] <- NA

    # Replace old annotations with new (where necessary)
    M    <- map[which(map[,2] %in% oldOddMirs),,drop=FALSE]
    idx  <- order(M[,1])
    M    <- M[idx,,drop=FALSE]
    fidx <- which(userMirs %in% M[,2])
    finalMirs[fidx] <- M[,1]
    
    # Replace new annotations with old (where necessary)
    M    <- map[which(map[,1] %in% newOddMirs),,drop=FALSE]
    idx  <- order(M[,2])
    M    <- M[idx,,drop=FALSE]
    fidx <- which(userMirs %in% M[,1])
    finalMirs[fidx] <- M[,2]

    res <- data.frame('Initial'=userMirs, 'Final'=finalMirs, 
            stringsAsFactors=FALSE)

    return(res)
}

