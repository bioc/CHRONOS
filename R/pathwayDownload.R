

downloadKEGGPathwayList <- function(org)
{
    link <- 'http://rest.kegg.jp/list/pathway/'
    if (!missing(org))
    {
        link <- paste0(link, org, '/', org, '.txt')
        rp   <- paste0('path:', org)
    }else
    {
        rp   <- 'path:map'
    }

    r  <- unlist(strsplit(getURIAsynchronous(link), '\n'))

    r  <- mapply(strsplit, as.list(r), split = '\t')
    r1 <- sapply(r, function(x){gsub(rp,'',x[1])})
    r2 <- lapply(r, function(x){x[2]})
    # Remove organism name
    r2 <- sapply(r2, function(x) { unlist(strsplit(x, ' -'))[1] }) 

    s  <- sprintf('%s\t%s',r1,r2)
    # pathways        <- r1
    # names(pathways) <- r2
    pathways           <- data.frame(r1, r2, stringsAsFactors = FALSE)
    colnames(pathways) <- c('Id', 'Name')

    return(pathways)
}

downloadPathways  <- function(org, pathways)
{
    if (missing(pathways)) { pathways <- 'All'}
    paths <- pathways
    if ( length(pathways) == 1 && 
        pathways %in% c('All', 'Metabolic', 'Non-Metabolic'))
    {
        # Download availiable KEGG pathways list 
        paths        <- downloadKEGGPathwayList(org)
        metabolicThr <- 2000

        if ('Metabolic' %in% pathways)
        {   
            paths <- paths[which(as.integer(paths[,1]) < metabolicThr), 1]
        }
        if ('Non-Metabolic' %in% pathways)
        {
            paths <- paths[which(as.integer(paths[,1]) >= metabolicThr), 1]
        }
        if ('All' %in% pathways)
        {
            paths <- paths[, 1]
        }
        if (!'Metabolic' %in% pathways && !'Non-Metabolic' %in% pathways 
                && !'All' %in% pathways)
        {
            if (class(pathways) == 'character')
            {
                paths <- pathways
            }
            if (class(pathways) == 'numeric' || class(pathways) == 'integer')
            {
                paths <- as.vector(t(paths[pathways, 1, drop=FALSE]))
            }
        }
    }

    # Download files
    message('Downloading files...', appendLF = FALSE)

    paths <- .downloadPathways(org, paths)

    message('done.')

    return(paths)
}

.downloadPathways <- function(org, paths)
{
    xmlDir <- cache$dirs$xml
    xmlDir <- paste(xmlDir, org, sep='//')

    # Create the organism subdirectory if it does not exist
    if (!file.exists(xmlDir)) 
    { 
        dir.create(xmlDir, recursive=TRUE)
    }

    xmlList <- paste(org, paths, sep='')

    # If any of the requested files have already been downloaded, skip them.
    downloadedPaths <- paste(xmlDir, list.files(path=xmlDir), sep='//')
    downloadQueue   <- unname(xmlList)
    downloadedFiles <- gsub('.xml','', list.files(path=xmlDir))
    idx             <- which(downloadQueue %in% downloadedFiles)
    if (length(idx) > 0) 
    { 
        downloadQueue <- downloadQueue[-idx] 
    }
    if (length(downloadQueue) == 0) 
    { 
        return(paths)
    }

    # Current http links for retrieving KEGG xml file
    prefix <- 'http://www.kegg.jp/kegg-bin/download?entry=' 
    suffix <- '&format=kgml'
    # Create download links
    links  <- paste0(prefix, downloadQueue, suffix)
    # Create destination paths
    dest   <- paste0(xmlDir, '//', downloadQueue, '.xml')
    
    # Split linear queue to batches of jobs
    N      <- length(downloadQueue)
    k      <- 4
    jobs   <- split(1:N, cut(1:N, quantile(1:N, (0:k)/k), include.lowest=TRUE,
                    labels=FALSE))
    export <- c('getURIAsynchronous', 'foreach', '%do%')
    cores  <- ifelse( N > detectCores()*10, 'default', 1 )
    lens   <- .doSafeParallel(
                funcName=downloadPathway,
                export=export, combine='c',
                N=k, cores=cores,
                links, dest, jobs)

    # Do some error-reporting
    errIdx <- which( lens < 3 )
    files  <- list('valid'=c(), 'invalid'=c())
    files$valid <- downloadQueue
    if ( length(errIdx) > 0 )
    {
        files$valid   <- downloadQueue[-errIdx]
        files$invalid <- downloadQueue[errIdx]
    }

    return(paths)
}

downloadPathway   <- function(i, ...)
{
    # ------- Unpacking arguments -------
    args  <- list(...)
    links <- args[[1]]
    dest  <- args[[2]]
    jobs  <- args[[3]]
    # ----------------------------------
    idx <- jobs[[i]]

    raw <- getURIAsynchronous(links[idx], noprogress=TRUE, 
                            .opts = list(timeout = 200, verbose=FALSE))

    lens <- vector(mode='numeric', length=length(raw))
    j <- 0; foreach (j = 1:length(raw)) %do% 
    {
        write(raw[[j]], file = dest[idx][[j]])
        lens[j] <- nchar(raw[[j]])
    }
    return(lens)
}

