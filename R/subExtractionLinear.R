

extractLinearSubpathways  <- function(graphs, pathways, a, b, filter, export, 
                                        groupMode, verbose)
{
    if (missing(a))         { a         <- 3 }
    if (missing(b))         { b         <- 10 }
    if (missing(export))    { export    <- c('.RData') }
    if (missing(filter))    { filter    <- FALSE }
    if (missing(groupMode)) { groupMode <- 'expand'}

    tryCatch(
    {
        message('Extracting Linear Subpathways...', appendLF = FALSE)        
        linSubs <- .extractLinearSubpathways(graphs, pathways, a, b, 
                                                filter, export, groupMode)
        message('done.')

        return(linSubs)
    }, warning = function(war) 
    {
    }, error = function(e) 
    {
        closeAllConnections()
        message('\nError in linear subpathway extraction.')
    }
    , finally = { }
    )    
}


.extractLinearSubpathways <- function(graphs, pathways, a, b, filter, 
                                        export, groupMode)
{
    if (filter)
    {
        # Extract user genes from gene expressions to keep only relevant 
        # subpaths
        load(cache$dirs$geneExpressions, e <- new.env()) 
        geneExpr  <- e[[ls(e)[1]]]

        if (is.null(geneExpr))
        { 
            stop('Please import expression data or set filter=FALSE\n')
            return(NULL) 
        }
    }

    t         <- 10000
    limit     <- 10
    opts      <- list(filter=filter, a=a, b=b, t=t, limit=limit, 
                    export=export)
    cAdjMats  <- graphs$combined
    eAdjMats  <- graphs$expanded
    org       <- graphs$org
    userGenes <- NULL

    if (!missing(pathways))
    {
        # Keep only graphs for selected pathways
        pathways <- paste(org, pathways, sep='')
        cAdjMats <- cAdjMats[pathways]
        eAdjMats <- eAdjMats[pathways]
    }

    if (length(cAdjMats) == 0) { return(NULL) }

    if (filter)
    {
        # Export genes from gene expressions
        userGenes <- sort(as.numeric(unique(rownames(geneExpr))))
    }

    # Create temporary directories
    .createDirectories(org)
        
    # Print compact adjacency matrices to disk (java input).
    matrixToFile(cAdjMats, org, cache$dirs)

    # Extract pathway names
    pathNames <- paste0(names(cAdjMats), '.txt')

    # Extract compact linear subpaths
    export <- c('getLinearPairs', 'cache', 'write.table', 'tail')
    cores  <- ifelse( length(cAdjMats) > detectCores()*10, 'default', 1 )
    res    <- .doSafeParallel( 
                funcName=getLinearSubpath, export=export,
                combine='default', N=length(cAdjMats), cores=cores,
                cAdjMats, pathNames, org, opts$a, opts$b, opts$t, opts$limit)

    # Summarize results
    if (length(res) > 0) 
    { 
        linRes <- data.frame( matrix(unlist(res), nrow=length(res), 
                            byrow=TRUE), stringsAsFactors=FALSE)    
    }else{ linRes <- NULL }


    if (groupMode == 'expand')
    {
        # Expand compact linear subpaths to obtain final subpaths.    
        subDir   <- paste(cache$dirs$tmp$sub, org, sep='//')
        subFiles <- paste(subDir, '//',names(cAdjMats), '.txt', sep='')

        export <- c('getSplittingFactor', 'expansionEvaluator',
                    'fullCollapse', 'subpathExpansion', 'fk', 'cleanSubpath')
        eSubs  <- .doSafeParallel(funcName=expandSubpath, 
                        export=export,
                        combine='default',
                        N=length(cAdjMats),
                        cores=cores,
                        cAdjMats, subFiles, opts$t, opts$b)
        names(eSubs) <- names(cAdjMats)

        # Delete temporaty files and directories
        .cleanDirectories()

        # Conform subpaths to a single matrix
        subpaths  <- unlistToMatrix(eSubs)

        if (filter)
        {
            # Filter subpaths using user genes from gene expression data
            subpaths <- .filterLinearSubpaths(subpaths=subpaths, 
                                            userGenes=userGenes)
        }

        # Export subpaths 
        exportSubpaths(subpaths, dir=paste(cache$dirs$lnr, org, sep='//'), 
                        type=opts$export, verbose=FALSE)

        adjMats <- eAdjMats    
    }

    if (groupMode == 'collapse')
    {
        # Collapse multiple genes to a single entry
        subpaths <- .collapseSubpaths(org, cAdjMats, b, userGenes)
        adjMats  <- cAdjMats    
    }

    # Remove concecutive duplicated genes from subpathways
    subpaths   <- .uniquifySubpaths(subpaths)
    subsLinear <- list(subpaths=subpaths, adjMats=adjMats, 
                        groupMode=groupMode, org=org, subpathwayType='Linear')


    return(subsLinear)
}



##
## Linear Subpathway Extraction
## 


getLinearSubpath      <- function(i, ...)
{
    # ------- Unpacking arguments -------
    args      <- list(...)
    adjMat    <- args[[1]][[i]]
    pathName  <- args[[2]][i]
    org       <- args[[3]]
    a         <- args[[4]]
    b         <- args[[5]]
    threshold <- args[[6]]
    limit     <- args[[7]]
    # ----------------------------------

    firstRun <- TRUE
    while(TRUE)
    {
        if (!firstRun) 
        {
            adjMat[upper.tri(adjMat)] <- 0 
        }

        parDir <- paste(cache$dirs$tmp$par, org, sep='//')
        lPairs <- getLinearPairs(adjMat, pathName, parDir)
        res    <- c(gsub('.txt', '', pathName), -2, -2)

        if (!is.null(lPairs$valid))
        {
            pathToJar <- system.file('java//CHRONOS.jar', package='CHRONOS')

            javaPath <- cache$dirs$javapath
            cmd <- paste(javaPath, '-cp', pathToJar, '-Xmx640m LinearPaths')
            baseDir <- paste(gsub('//','/',cache$dirs$int),'/', sep='')
            cmd <- paste(cmd, lPairs$valid, org, baseDir, a, b, 
                                threshold, limit)
            res <- system(cmd, intern=TRUE, ignore.stderr=FALSE)

            if (length(res) > 0)
            {
                res[1] <- gsub('.txt', '', 
                                tail(unlist(strsplit(res[[1]], '/')), n=1))
            }else
            {
                res <- c(gsub('.txt', '', pathName), -1, -1)
            }
        }

        if (res[2] != '0') { break() }
        if (res[2] == '0' && firstRun)  { firstRun <- FALSE; next() }
        if (res[2] == '0' && !firstRun) { break() }
    }

    return (res)
}

getLinearPairs        <- function(adjMat, pathNames, pairDir)
{

    sourceNodes  <- vector()
    destinNodes  <- vector()
    valid = c(); invalid = c()

    # Source nodes do not have any ones in their respective columns
    # Destination nodes do not have ones in their respective rows.
    if(!is.matrix(adjMat)) 
    {
        return (list(valid=NULL, invalid=pathNames))
    }


    for (i in 1:ncol(adjMat))
    {
        # A source node has no incoming edges
        isSource = ( length(which(adjMat[,i] == 0)) == ncol(adjMat) )
        # A destination node has no outgoing edges
        isDestin = ( length(which(adjMat[i,] == 0)) == nrow(adjMat) )

        # If it both source and destination node it is a remote node
        if (isSource && isDestin) { next() }

        if (isSource) { sourceNodes <- c(sourceNodes, i) }
        if (isDestin) { destinNodes <- c(destinNodes, i) }
    }

    res <- list()
    ctr <- 1
    for ( s in sourceNodes )
    {
        source <- unname(s)
        for ( d in destinNodes )
        {
            destin     <- unname(d)
            if (source != destin) 
            {
                res[[ctr]] <- list(source, destin)
                ctr        <- ctr + 1
            }
        }
    }

    if (!is.null(unlist(res))) 
    { 
        results <-  t(matrix(unlist(res), nrow = 2)) 
    }else 
    { 
        results = '' 
    }


    # Print each set of valid pairs to file
    if ( is.matrix(results) ) 
    {
        valid <- c(valid, pathNames)
    }
    else
    {
        results <- as.matrix(vector(), ncol=2)
        invalid <- c(invalid, pathNames)
    }

    outfile <- paste(pairDir,  pathNames, sep = '//')
    write.table( results, file = outfile, sep = '\t', col.names=FALSE, 
                row.names=FALSE )
    files = list(valid=valid, invalid=invalid)

    return(files)
}

filterLinearSubpaths  <- function(subpaths, userGenes)
{

    return(.filterLinearSubpaths(subpaths, userGenes))
}

.filterLinearSubpaths <- function(subpaths, userGenes)
{
    if ( is.null(subpaths) )         { return(subpaths) }
    if ( is.null(userGenes) )        { return(subpaths) }
    if ( nrow(subpaths)==0 )         { return(subpaths) }


    # Conform subpaths to a single matrix
    if (is.list(subpaths)) 
    { 
        subpaths <- unlistToMatrix(subpaths) 
    }

    userGenes     <- c(0, userGenes)
    r             <- matrix(1, nrow=nrow(subpaths), ncol=ncol(subpaths))
    R             <- subpaths
    for (i in 1:ncol(subpaths))
    {
        idx <- which(subpaths[,i]  %in% userGenes)
        if (length(idx) > 0)     { r[idx, i] <- 0 }
        if (length(idx) == 0)     { r[, i]     <- 0 }
    }
    idx <- which(rowSums(r) > 0)

    if (length(idx) > 0) 
    {
        R           <- subpaths[-idx, ,drop=FALSE]
        rownames(R) <- rownames(subpaths)[-idx]
    }
            
    return(R)
}



## 
## Linear Subpathway Collapsing
## 

.collapseSubpaths    <- function(org, adjMats, b, userGenes)
{

    message('Collapsing Linear Subpaths...', appendLF = FALSE)

    pathNames <- names(adjMats)
    subDir    <- paste(cache$dirs$tmp$sub, org, sep='//')
    subFiles  <- paste(subDir, '//',pathNames, '.txt', sep='')
    subpaths  <- matrix(,nrow=0, ncol=b)
    
    for (i in 1:length(subFiles))
    {
        if (!file.exists(subFiles[i])) { next() }

        if (file.exists(subFiles[i])) 
        {
            cSubpaths <- scan(subFiles[i], what='', sep="\n", quiet=TRUE)
            if (length(cSubpaths) == 0)  { next() }
        }
        
        genes <- vector(mode='list', length=ncol(adjMats[[i]]))
        z     <- 1
        for (col in colnames(adjMats[[i]]))
        {
            genes[[z]] <- as.integer(unlist(strsplit(col, ' ')))
            z          <- z + 1
        }

        # Check each element (single gene or group of genes) 
        # and remove genes not in user list
        if (!missing(userGenes))
        {            
            for ( j in 1:length(genes) )
            {
                idx <- which(genes[[j]] %in% userGenes)
                if (length(idx) > 0)  { genes[[j]] <- genes[[j]][idx] }
                if (length(idx) == 0) { genes[[j]] <- NA }
            }
        }

        if (length(cSubpaths) > detectCores() * 10)
        {
            # Create and register parallel backend
            cl <- makeCluster(detectCores(logical=FALSE))
            on.exit(stopCluster(cl), add=TRUE)
            registerDoParallel(cl)
        }else
        {
            # Run sequentially for trivial inputs
            registerDoSEQ()
        }

        subs <- foreach(j=1:length(cSubpaths), .combine='rbind') %dopar%
        {
            # Create a subpath as a list of bins
            sub <- as.integer(unlist(strsplit(cSubpaths[[j]], ' ')))
            sub <- c(genes[sub], as.list(rep(0, b-length(sub))))

            # Reject subpath if it has a invalidated gene (not in user list)
            if(length(which(is.na(sub))) > 0) { return(NULL) }

            sub <- sapply(sub, function(x) { paste(x, collapse=' ') })
            sub <- matrix(sub, nrow=1)
            rownames(sub) <- names(adjMats)[i]

            return(sub)
        }

        subpaths <- rbind(subpaths, subs)
    }

    message('done.')

    return(subpaths)
}


##
## Linear Subpathway Expansion
##


expandSubpath      <- function( i, ... )
{
    # ------- Unpacking arguments -------
    args      <- list(...)
    adjMat    <- args[[1]][[i]]
    subFile   <- args[[2]][i]
    threshold <- args[[3]]
    b         <- args[[4]]
    # ----------------------------------

    # Read compact subpathway from local file of created subpaths 
    if (file.exists(subFile))
    {
        cSubpaths <- scan(subFile, what='', sep="\n", quiet=TRUE)
    }else { return(NULL) }
    
    j     <- 1
    genes <- vector(mode='list', length=ncol(adjMat))
    for (col in colnames(adjMat))
    {
        genes[[j]] <- as.integer(unlist(strsplit(col, ' ')))
        j <- j + 1
    }

    sf    <- getSplittingFactor(cSubpaths, genes, b, threshold)
    k     <- sf$k
    nsubs <- sf$nsubs
    genes <- sf$genes
    
    # Full collapsing
    if (k < 1)
    {
        res <- fullCollapse(cSubpaths, adjMat, b, threshold)
    }
    # Partial collapsing
    if (k > 0)
    {
        R   <- subpathExpansion(genes, cSubpaths, nsubs, threshold,
                                k=k, b=b)
        res <- R$eSubpaths
    }

    return(res)
}

getSplittingFactor <- function(subs, genes, b, threshold)
{
    if (length(subs) > 2000) 
    { 
        return( list(k=0, nsubs=length(subs), genes=genes) ) 
    }
    
    while(TRUE)
    {
        k  <- 1
        s0 <- 0
        while(TRUE)
        {
            j     <- 1
            nsubs <- vector(mode="numeric", length=length(subs))
            for (sub in subs)
            {
                nsubs[j<-j+1] <- expansionEvaluator(sub, genes, k=k, b=b, 
                                                    threshold)
            }
            s  <- sum(nsubs)
            if ( s <= threshold && k<2 )  { s0 <- s;   break() }
            if ( s <= threshold && k>1 )  { s0 <- s;   break() }
            if ( s > s0 && k > 1 )        { k  <- k-1; break() }
            k  <- k + 1
            s0 <- s
        }

        if (s0 > threshold)        
        { 
            z <- 1
            genesCut <- vector(mode='list', length=length(genes))
            for (x in genes)
            {
                s             <- 1:ceiling(length(x)/2)
                genesCut[[z]] <- x[s]     
                z             <- z + 1
            }
            genes <- genesCut
            k     <- 1 
            s0    <- 0
            next() 
        }else { break() }
    }


    return(list(k=k,nsubs=s0,genes=genes))
} 

expansionEvaluator <- function(compactSubpath, genes, k, b, threshold)
{
    if (is.null(compactSubpath)) { return(NULL) }
    
    # Create a subpath as a list of bins
    sub <- as.integer(unlist(strsplit(compactSubpath, ' ')))
    sub <- c(genes[sub], as.list(rep(0, b-length(sub))))
        
    # Spit the compact subpath to multiple components
    s   <- vector(mode='list', length=k)        
    for (j in 1:k) { s[[j]] <- lapply(sub, fk, k, j) }

    stats <- 0
    for (j in 1:k) 
    {
        gctr <- 1
        sub  <- s[[j]]
        for (z in 1:length(s[[j]]))
        { 
            if (length(sub[[z]]) > 1) { gctr <- gctr * length(sub[[z]]) } 
        }
        stats <- stats + gctr
        if (stats > threshold) { break() }
    } 

    return(stats=stats)
}

fullCollapse       <- function(compactSubpaths, adjMat, b, threshold)
{
    j     <- 1
    genes <- vector(mode='numeric', length=ncol(adjMat))
    for (cnames in colnames(adjMat))
    {
        x        <- unlist(strsplit(cnames, ' '))
        s        <- 1
        genes[j] <- as.integer(x[s])          
        j        <- j + 1
    }

    N <- length(compactSubpaths)

    if (N > threshold) 
    { 
        s               <- 1:threshold
        compactSubpaths <- compactSubpaths[s]   
        N               <- length(compactSubpaths)
    }
    eSubpaths <- matrix(0, nrow=N, ncol=b)
    j <- 1
    for (compactSubpath in compactSubpaths)
    {
        sub             <- as.integer(unlist(strsplit(compactSubpath, ' ')))
        sub             <- c(genes[sub], rep(0, b-length(sub)))
        eSubpaths[j,]   <- sub
        j               <- j + 1
    }

    return(eSubpaths)
}

subpathExpansion   <- function(genes, compactSubpaths, nsubs, threshold, 
                                k, b)
{
    if(length(compactSubpaths) == 0) { return(NULL) }
    
    rowCtr    <- 0
    eSubpaths <- matrix(0, nrow=nsubs, ncol=b)
    for (compactSubpath in compactSubpaths)
    {
        # Create a subpath as a list of bins
        sub <- as.integer(unlist(strsplit(compactSubpath, ' ')))
        sub <- c(genes[sub], as.list(rep(0, b-length(sub))))
        # Remove a possible second occurences of a gene in consecutive bins.
        sub <- cleanSubpath(sub)

        # Spit the compact subpath to multiple components
        s   <- vector(mode='list', length=k[1])
        for (i in 1:k) { s[[i]] <- lapply(sub, fk, k, i) }

        # Expand components
        for (i in 1:k) 
        {
            m <- t(t(expand.grid(s[[i]])))
            m <- unname(split(m, 1:nrow(m)))
            if ( !is(m, 'matrix') ) 
            { 
                m <- matrix(unlist(m), ncol=b, byrow=FALSE)
            }
            if ( is(m, 'matrix') )
            { 
                m <- matrix(unlist(m), ncol=b, byrow=TRUE) 
            }
            eSubpaths[(rowCtr+1):(rowCtr+nrow(m)),] <- m
            rowCtr <- rowCtr + nrow(m)
        }
    }
    # Remove exceeding rows (produced by emptying successive bins 
    # with same gene ids)
    eSubpaths <- eSubpaths[1:rowCtr,]
    if (!is.matrix(eSubpaths)) { eSubpaths <- matrix(eSubpaths, ncol=b) }
    

    return(list(eSubpaths=eSubpaths))
}

fk                 <- function(x,k,i)
{
    if (length(x) < k)  
    { 
        r <- x[1:length(x)] 
    }
    if (length(x) >= k) 
    { 
        r <- x[(floor((i-1)*length(x)/k)+1) : floor(i*length(x)/k)] 
    }
    return(r)
}

cleanSubpath       <- function(x)
{
    for (i in 1:(length(x)-1))
    {
        a <- x[[i]]
        b <- x[[i+1]]
        S <- intersect(a,b)
        for (s in S)
        {
            if (length(a) == 1 & length(b) > 1) 
            { 
                b <- b[-which(b == s)]; next() 
            }
            if (length(a) > 1 & length(b) == 1) 
            { 
                a <- a[-which(a == s)]; next() 
            }
            if (length(a) > 1 & length(b) > 1)  
            { 
                b <- b[-which(b == s)]; next() 
            }
        }
        x[[i]]   <- a
        x[[i+1]] <- b
    }

    return(x)
}


#
# Subpathway cleanup
#

.uniquifySubpaths   <- function(subpaths)
{

    if (is.null(subpaths))    { return(NULL) }
    if (nrow(subpaths) == 0)  { return(subpaths) }

    # Remove duplicates
    subpaths = subpaths[which(!duplicated(subpaths)),  , drop=FALSE]

    # Split subpaths according to their lengths
    len  <- getLengths(subpaths)
    didx <- c()

    for (i in 1:max(len))
    {
        idx0     <- which(len ==i)
        subs     <- subpaths[idx0, 1:i, drop=FALSE]

        if (nrow(subs) == 0) { next() }

        dupl     <- t(apply( subs , 1 , diff ))
        idx1     <- unique(which(dupl == 0, arr.ind = TRUE)[,'row'])
        didx     <- c(didx, idx0[idx1])
    }

    if (length(didx) > 0)
    {
        for (i in 1:length(didx))
        {
            idx     <- didx[i]
            ugenes  <- subpaths[idx, cumsum(rle(subpaths[idx,])$length)]
            subpaths[idx,]                  <- rep(0, ncol(subpaths))
            subpaths[idx, 1:length(ugenes)] <- ugenes
        }
    }

    # Remove too short subpaths
    lens     <- getLengths(subpaths)
    idx     <- which(lens < 3)
    if (length(idx) > 0)
    {
        subpaths <- subpaths[-idx, , drop=FALSE]
    }

    return(subpaths)
}


