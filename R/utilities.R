#
# Data consolidation
#


fillMatrixList <- function(L, maxLen)
{
    if (is.null(L))     { return(L) }

    # If input is not a list, return original data
    if (!is.list(L))    { return(L) }

    lens <- c()
    for (i in 1:length(L))
    {    
        if (!is.null(L[[i]]))
        {
            lens <- c(lens, getLengths(L[[i]]))
        }
    }

    if (is.null(lens)) { return(NULL) }


    # Set number of columns.
    if (missing(maxLen)) { maxLen <- max(lens, na.rm=TRUE) }


    # If no filling is necessary, return original data
    if (length(unique(lens)) == 1) { return(L) }

    # Enforce maximum length on each subpath of each pathway
    res <- L
    for (i in 1:length(L))
    {
        submat <- L[[i]]
        if (is.null(submat)) { next() }
        replmat <- matrix(0, nrow=nrow(submat), ncol=maxLen)
        for (j in 1:nrow(submat))
        {   
            replmat[j,] <- c(submat[j,], rep(0, maxLen - length(submat[j,])))
        }
        res[[i]] <- replmat
    }

    return(res)
}

unlistToMatrix <- function(L, mode='rbind')
{
    if (is.null(L))     { return(L) }

    # If input is already a matrix, return original data.
    if ( is(L,'matrix') ) { return(L) }

    # Find what type of data the list holds
    type <- NULL
    for (i in 1:length(L))
    {
        if(!is.null(L[[i]]))
        {
            if (is.vector(L[[i]])) type <- 'vector'
            if (is.matrix(L[[i]])) type <- 'matrix'
            break()
        }
    }

    if (is.null(type)) { return(NULL) }

    # List of vectors of variable size
    if (type == 'vector')
    {
        lens <- sapply(L, function(x) { length(x) })
        M    <- matrix(0, nrow=length(L), ncol=max(lens))
        for (i in 1:length(L)) 
        {
            M[i,1:lens[i]] <- L[[i]] 
        }
        rownames(M) <- names(L)
    }

    # A list of matrices with at least one common dimension size
    # rbind:same number of columns, cbind:same number of rows
    if (type == 'matrix')
    {
        M      <- do.call(mode, L)
        lnames <- names(L) 

        if (length(lnames) > 0)
        {
            lens   <- sapply(L, function(x) 
                                { 
                                    if (!is.null(x)) nrow(x) else 0 
                                } )
            rnames <- vector(mode='numeric', length=nrow(M))
            ctr   <- 0
            for (i in 1:length(lnames))
            {
                if (lens[i] > 0)
                {
                    idx         <- (ctr+1) : (ctr<-ctr+lens[i])
                    rnames[idx] <- rep(lnames[i], lens[i])               
                }
            }
            rownames(M) <- rnames
        }
    }
    return(M)
}

getLengths     <- function(M)
{
    # 
    # Count length of non zero elements in each row of a matrix.
    # Each row must consist of a prefix with non zero elements, 
    # and a suffix of consecutive zeros denote absence of data.
    # Ideal for matrices with a large number of rows.
    #
    esub <- cbind(M, rep(0, nrow(M)))
    esub <- rbind(esub, c(-1, rep(0, ncol(M))))
    df   <- which(diff(which(c(t(esub)) != 0)) - 1 > 0)
    len  <- c(df[1], diff(df))

    return(len)
}

splitWork      <- function(N, k)
{
    if (missing(k))
    {
        k <- detectCores(logical=FALSE)
    }
    if (N < k) 
    { 
        k <- N 
    }
    
    if (N > 1) 
    {
        jobs <- split(1:N, cut(1:N, quantile(1:N,(0:k)/k), 
                    include.lowest=TRUE, labels=FALSE))
    }
    if (N == 1)
    {
        jobs <- list(c(1)); names(jobs) <- '1'
    }

    return(list(jobs=jobs, k=k))
}


matrixToList   <- function(subpaths, block, k)
{
    if ( is(subpaths,'list') )
    {
        subpaths <- unlistToMatrix(subpaths)
    }
    if (missing(block)) { block <- 1000 }
    if (missing(k))     { k <- detectCores(logical=FALSE) }

    K   <- floor(nrow(subpaths)/block)
    xx  <- splitWork(nrow(subpaths), K)
    S   <- vector(mode='list', length=K)
    for (i in 1:K)
    { 
        S[[i]] <- subpaths[xx$jobs[[i]],]
    }

    kk  <- splitWork(length(S), k=k)
    res <- list(subpaths=S, jobs=kk$jobs, k=kk$k)

    return(res)
}


matrixToFile   <- function(adjMats, org, dirs )
{
    # Remove possible edge type information
    for (i in 1:length(adjMats))
    {
        adjMats[[i]][adjMats[[i]] > 1] <- 1
    }

    adjMatFiles <- paste( dirs$tmp$mat,'//', org, '//', names(adjMats) ,
                        '.txt', sep = '' )
    IdsFiles    <- paste( dirs$tmp$ids,'//', org, '//', names(adjMats) ,
                        'ids.txt', sep = '' )

    mapply( write.table, adjMats, file = adjMatFiles, sep = ",", 
            row=FALSE, col=FALSE )
    matIds <- lapply(adjMats, function(x) { colnames(x) } )
    mapply( write.table, matIds, file = IdsFiles, sep = ",", 
        row=FALSE, col=FALSE )

    return(NULL)
}



optimalLengths <- function(subpaths, mode)
{
    if (is.list(subpaths)) { subpaths <- unlistToMatrix(subpaths) }

    res  <- vector(mode='numeric', length=ncol(subpaths)) 
    lens <- getLengths(subpaths)
    for (len in 3:ncol(subpaths))
    {
        idx      <- which(lens <= len)
        res[len] <- length(unique(as.vector(subpaths[idx,])))
    }

    return(data.frame('optimalLengths'=res))
}


idxToMat       <- function(L, idx)
{
    rows <- ((idx-1) %% nrow(L)) + 1
    cols <- floor((idx-1)/nrow(L)) + 1

    return(list(rows=rows, cols=cols))
}


normalize      <- function(x)
{
    res <- (x-min(x)) / (max(x)-min(x))
    return(res)
}





#
# Data Export
#

exportSubpaths   <- function(subpaths, dir, type, verbose, mode)
{
    if (is.null(subpaths))  { return() }
    if (missing(verbose))   { verbose <- TRUE }
    if (missing(mode))      { mode <- 'clean'}
    if (missing(type))      { type <- c('.txt','.RData')}

    if (verbose)
    { 
        message('Exporting subpaths...', appendLF = FALSE) 
    }
    
    # Determine max serial number of output files already present in directory
    n1 <- 0; n2 <- 0
    f1 <- gsub('.RData', '', list.files(dir, pattern='.RData'))
    f2 <- gsub('.txt', '', list.files(dir, pattern='.txt'))
    if (length(f1) > 0)
    {
        n1 <- as.integer(unique(unlist(strsplit(f1, 'subpaths'))))
        n1 <- max(n1, na.rm=TRUE) 
    }
    if (length(f2) > 0)
    {
        n2 <- as.integer(unique(unlist(strsplit(f2, 'subpaths'))))
        n2 <- max(n2, na.rm=TRUE)
    }
    n  <- max(n1,n2)
        
    if ('.RData' %in% type)
    {   
        file <- paste0(dir, '//subpaths', sprintf("%03d", n+1),'.RData')
        file.create(file)
        save(subpaths, file=file)
    }
    if ('.txt' %in% type)
    {
        # Create output file
        file <- paste0(dir,'//subpaths', sprintf("%03d", n+1),'.txt')
        file.create(file)
        i      <- 1
        paths  <- names(subpaths)
        if (is.list(subpaths))
        {
            for (sub in subpaths)
            {
                # Append subpaths to output file
                if (is.null(sub))
                { 
                    i <- i + 1; next()
                }
                sub <- cbind(rep(paths[i], nrow(sub)), sub)
                write.table(sub, file=file, sep='\t', append=TRUE, 
                                row.names=FALSE, col.names=FALSE, quote=FALSE)
                i <- i + 1
            }
        }
        if (is.matrix(subpaths))
        {
            # Append subpaths to output file
            if (mode != 'justify')
            {
                write.table(subpaths, file=file, sep='\t', append=TRUE, 
                            row.names=TRUE, col.names=FALSE, quote=FALSE)   
            }
            if (mode == 'justify')
            {
                exportToFile(R=subpaths, expfile=file)  
            }
        }
    }

    if (verbose) { message('done.') }

    return(invisible())
}


createScoreFile  <- function(dir, type, fileType)
{
    if (missing(fileType)) { fileType <- '.txt'}

    # Determine max serial number of output files already present in directory
    n1 <- 0; n2 <- 0
    f1 <- gsub(paste('a', fileType, sep=''), '', 
                        list.files(dir, pattern=paste('a', fileType, sep='')))
    f2 <- gsub(paste('b', fileType, sep=''), '', 
                        list.files(dir, pattern=paste('b', fileType, sep='')))

    if (length(f1) > 0) 
    {
        N1 <- as.integer(unique(unlist(strsplit(f1, 'scores'))))
        n1 <- max(n1, na.rm=TRUE)
    }
    if (length(f2) > 0)
    {
        n2 <- as.integer(unique(unlist(strsplit(f2, 'scores'))))
        n2 <- max(n2, na.rm=TRUE) 
    }
    if (type == 'a')    { n <- n1 }
    if (type == 'b')    { n <- n2 }

    file <- paste0(dir,'//scores', sprintf("%03d", n+1), 
                    paste0(type, fileType))
    file.create(file)

    return(file)
}


exportToFile    <- function(R, expfile, append=FALSE)
{
    if (!file.exists(expfile)) {  file.create(expfile) }
    max.print <- getOption('max.print')
    width     <- getOption('width')
    options(width = 10000)
    if ( is(R,'matrix') )
    {
        options(max.print=nrow(R) * ncol(R))
        sink(file=expfile, append=append)
        # print.data.frame(noquote(as.data.frame(R)), row.names=T, q=F)
        print(noquote(R), row.names=TRUE, col.names=FALSE, q=FALSE)
        sink(NULL)
    }
    if ( is(R,'list') )
    {
        options(max.print=length(unlist(R)))
        sink(file=expfile); 
        print(R); sink();
    }
    if ( is(R,'numeric') )
    {
        x <- paste(names(R), '\t', R, collapse='\n')
        cat(x, file=expfile , sep = '\n')
    }
    if ( is(R,'integer') )
    {
        x <- paste(R, '\t', names(R), collapse='\n')
        cat(x, file=expfile , sep = '\n')
    }
    if ( is(R,'data.frame') )
    {
        options(max.print=nrow(R) * ncol(R))
        sink(file=expfile, append=append)
        print(noquote(R), row.names=FALSE, right=FALSE)
        sink()
        options(max.print=max.print)
    }
    if ( is(R,'character') )
    {
        options(max.print=length(R))
        sink(file=expfile)
        print(noquote(R), row.names=FALSE)
        sink()
        options(max.print=max.print)
    }
    options(width = width)
}



checkFile       <- function(dir, file, filepath='')
{
    if (filepath=='')
    {
        filepath <- paste(dir, file, sep='//')
    }
    if (!file.exists(filepath))
    {
        message('Wrong or missing filepath', filepath);
        return(FALSE)
    }
    else
    {
        return(TRUE)
    }
}



# 
# Basic fault tolerant parallel computation 
# 

.doSafeParallel <- function(funcName, export, combine, N, cores, ...)
{
    if ( missing(export) )    { export <- NULL }
    if ( cores == 'default' ) { cores <- detectCores(logical=FALSE) }


    if ( cores == 1 )
    {
        # Run sequentially for trivial inputs
        registerDoSEQ()
        exports <- NULL
    }else
    {
        # Create and register parallel backend
        cl <- makeCluster(cores)
        on.exit(stopCluster(cl), add=TRUE)
        registerDoParallel(cl)
        exports <- c(as.character(substitute(funcName)), export)
        env     <- parent.env(parent.frame())

        if ( !is.null(export) )
        {
            clusterExport(cl=cl, varlist=exports, envir=env)
        }
    }


    i <- 0
    while (TRUE)
    {
        tryCatch(
        {
            if ( combine == 'default')
            {
                res <- foreach( i=1:N, .export=exports) %do%
                {
                    funcName(i, ...)
                }
            }else
            {
                res <- foreach( i=1:N, .export=exports, 
                                .combine=combine) %do%
                {
                    funcName(i, ...)
                }
            }
            break;
        }, error = function(e) 
        {
            closeAllConnections()
            message('\n\t Master cannot establish communication with workers')
        })

        message('\t+ Reinitializing clusters...', appendLF = FALSE)
        on.exit(NULL)
        cl <- makeCluster(detectCores(logical=FALSE))
        on.exit(stopCluster(cl), add=TRUE)
        registerDoParallel(cl)
        if ( !is.null(export) )
        {
            clusterExport(cl=cl, varlist=exports, envir=env)
        }
        message('done.')
    }

    return(res)
}

#
# Obtain full path to an executable file
#
.Sys.which2 <- function(cmd) 
{
    stopifnot(length(cmd) == 1)
    if (.Platform$OS.type == "windows") 
    {
        suppressWarnings({
            pathname <- system(sprintf("where %s", cmd), intern=TRUE)[1]
        })
        if (!is.na(pathname)) return(setNames(pathname, cmd))
    }
    Sys.which(cmd)
}

