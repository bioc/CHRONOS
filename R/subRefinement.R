
##
## Score subpathways using an interacting scoring scheme (IS) for each time 
## point
##

scoreSubpathways   <- function(subpathways, filters, measures, parameters,
                                miRNAinteractions)
{
    tryCatch(
    {        
        res <- .scoreSubpathways(subpathways, filters, measures, parameters,
                                miRNAinteractions)
        return(res)
    }, warning = function(war) 
    {
    }, error = function(e) 
                { 
                    closeAllConnections()
                    message('\nParall error in refinement'); 
                }
    , finally = { }
    )
    
    return(res)
}

.scoreSubpathways  <- function(subpathways, filters, measures, parameters,
                                miRNAinteractions)
{
    if (missing(filters))        
    { 
        stop('Missing filters.')
    }
    if (missing(parameters))     
    { 
        parameters <- c('C'=1, 'K'=5, 'T'=0.3) 
    }
    if (missing(miRNAinteractions))        
    { 
        miRNAinteractions <- NULL
    }
    if ('measures' %in% names(filters) && missing(measures))
    {
        stop('Measures parameter missing.')
    }
    if (!'subScore' %in% names(filters))
    {
        filters <- c(filters, 'subScore'=0.0)
    }
    if ( 'mirScore' %in% names(filters) )
    {
        if ( is.null(miRNAinteractions) )
        {
            # Download miTargets file if it is not availiable
            org    <- subpathways$org
            miFile <- paste(cache$dirs$miDir, org, cache$dirs$miFile, sep='//')
                    
            if (!file.exists(miFile))
            {
                downloadMiRecords(org=org, update=FALSE, databases='All', pn=5)
            }
        }
        if ( !is.null(miRNAinteractions) )
        {
            org <- subpathways$org
            saveMiRNAFile(org=org, miRNAinteractions=miRNAinteractions)
        }
    }

    # Load imported gene expressions
    load(cache$dirs$geneExpressions, e <- new.env())
    geneExpr <- e[[ls(e)[1]]]

    if ( is.null(geneExpr) ) 
    { 
        stop('Error importing gene expression data.\n')
        return(NULL) 
    }

    if ( 'mirScore' %in% names(filters) )
    {
        if(!file.exists(cache$dirs$mirnaExpressions))
        {
            message('Error, no miRNA expressions data supplied.\n')
            return(NULL)
        }

        # Load imported miRNA expressions
        load(cache$dirs$mirnaExpressions, e <- new.env())

        mirnaExpr <- e[[ls(e)[1]]]
        if (is.null(mirnaExpr))
        {
            stop('Error importing miRNA expression data.\n')
            return(NULL)
        }         
    }

    res <- .evaluateSubpaths(subpathways, filters, measures, parameters,
                            geneExpr, mirnaExpr)

    return(res)
}

.evaluateSubpaths  <- function(subpathways, filters, measures, parameters, 
                                geExpr, miExpr)
{
    subpaths  <- subpathways$subpaths
    adjMats   <- subpathways$adjMats
    groupMode <- subpathways$groupMode
    type      <- subpathways$subpathwayType
    org       <- subpathways$org
    fci.C     <- parameters['C']
    fci.K     <- parameters['K']
    fci.T     <- parameters['T']
    res       <- list(org=org, filters=filters, type=type)

    if (is.null(subpaths)) { return(NULL) }

    if (type == 'Linear')
    {
        # Filter adjacency matrices with user genes
        adjMats <- filterMatrix( adjMats=adjMats, org=org, 
                                userGenes=rownames(geExpr) )

        # Build global gene and edge indexer
        lex <- matrixToLexicon( adjMats=adjMats, org=org, groupMode=groupMode )
    }

    if (type == 'Non-Linear')
    {

        lex <- subpathways$lex
    }
    

    # Filter gene expressions
    geExpr     <- .filterExpressions( geExpr, lex, groupMode )

    message('\t#', nrow(subpaths), ' subpathways initially')


    if ('p-value' %in% names(filters))         
    {
        # Find most significant subpaths using hypergeometric testing
        qValues     <- .getStatistics(subpaths, names(lex$genes), type=type, 
                                    pVal=filters['p-value'], pAdjust = '')
        # Keep subpaths with good q-value
        idx0         <- which(qValues < 0.05) 
        subpaths     <- subpaths[idx0, , drop=FALSE]
        qValues     <- qValues[idx0]
        res          <- c(res, list(pValues=qValues))

        message('\t#', nrow(subpaths), ' subpathways after p.value.')
    }


    if ('subScore' %in% names(filters))    
    {
        if( nrow(subpaths) == 0 ) { return(NULL) }
                
        message('Calculating subscores...', appendLF = FALSE)

        # Index subpaths
        edgeList     <- lex$edges
        iSubs        <- .indexSubpaths(subpaths, edgeList, type, org)

        # Score only subpath interactions    
        uEdges       <- as.integer(unique(as.vector(iSubs)))
        uEdges       <- uEdges[which(uEdges != 0)]
        edgeList     <- edgeList[uEdges,, drop=FALSE]

        # Calculate edge scores for all subpathway edges
        edgeScores  <- .getEdgeScores(expr=geExpr, edges=edgeList, 
                        fci.C=fci.C, fci.K=fci.K, fci.T=fci.T)

        # Calculate subscore for each subpathway
        subScores   <- .scoreSubpaths(edgeScores=edgeScores, iSubs=iSubs, 
                                    type=type)

        # Keep subpaths with specific subscore
        idx1        <- .filterSubScore(subpaths, subScores, 
                                        filters['subScore'])
        subpaths    <- subpaths[idx1, , drop=FALSE]
        subScores   <- subScores[idx1, , drop=FALSE]
        res$pValues <- res$pValues[idx1] # NULL if not specified.


        message('done.')

        message('\t#', nrow(subpaths), ' subpathways after subscore.')

        res <- c(res, list(subScores=subScores))
    }

    if ('mirScore' %in% names(filters))
    {
        if( nrow(subpaths) == 0 ) { return(NULL) }

        message('Calculating mirscores...', appendLF = FALSE)
        
        miLex       <- .createMirnaLexicon(org=org, genes=lex$genes)
        
        # Score subpathway specific interactions
        subGenes    <- getSubpathwayGenes(subpaths, type)
        idx         <- which(miLex$edges[,4] %in% subGenes)
        miLex$edges <- miLex$edges[idx, ]
        edgeList    <- miLex$edges
        
        edgeScores  <- .getMiScores(geExpr, miExpr, edgeList, lex, type, 
                            fci.C=fci.C, fci.K=fci.K, fci.T=fci.T)

        # Append edge scores as third column of mirna-gene edgelist
        miLex$edges <- cbind(miLex$edges, edgeScores)
        edgeList    <- miLex$edges[, -c(1,2)]

        summary     <- .filterMiScores(subpaths, miLex, filters['mirScore'], 
                            type=type, lex=lex) 
        res         <- c(res, summary)

        message('...end.')
    }


    if ('measures' %in% names(filters))
    {
        # Apply topological evaluator to data
        top <- .getMeasures(names(lex$genes), subpaths, measures, mode)
        res <- c(res, list(measures=top))
    }

    res <- c(res, list(subpaths=subpaths))

    return (res)
}

.filterExpressions <- function( expressions, lex, groupMode )
{

    # Keep only gene expressions of genes belonging to user pathway list
    idx  <- which(rownames(expressions) %in% names(lex$genes))
    expr <- expressions[idx, , drop=FALSE]

    if (groupMode == 'collapse')
    {    
        # Find nodes with multiple genes and organize them to a matrix
        gn  <- lapply(as.list(names(lex$genes)), function(x) 
                    { matrix(unlist(strsplit(x, ' ')), nrow=1) } )
        len <- sapply(gn, function(x) { ncol(x) } )

        d     <- 0: max(len)
        idxT  <- which(len <= d[length(d)])
        nexp  <- matrix(, nrow=length(idxT), ncol=ncol(expr))
        rownames(nexp) <- 1:nrow(nexp)
        for (j in 1:(length(d)-1))
        {
            idx <- intersect(which(len > d[j]), which(len <= d[j+1]))
            if (length(idx) == 0) 
            { 
                next() 
            }
            gni <- unlistToMatrix(fillMatrixList(gn[idx]))

            # Keep relevant gene expressions
            ugenes  <- unique(as.vector(gni))
            idx     <- which(rownames(expressions) %in% ugenes)
            expr    <- expressions[idx, , drop=FALSE]
            
            nei     <- matrix(, nrow=nrow(gni), ncol=ncol(expr))
            for (t in 1:ncol(expr))
            {
                lb        <- expr[, t]
                names(lb) <- rownames(expr)
                expr_t    <- matrix(lb[gni], nrow=nrow(gni), ncol=ncol(gni))
                nei[,t]   <- apply( expr_t, 1, function(x) 
                                                { median(x, na.rm=TRUE ) } )
            }
            rownames(nei) <- apply( gni, 1, function(x)
                                { 
                                    paste(x[which(x!=0)], collapse=' ') 
                                } )
            idx1                <- which(len[idxT] > d[j])
            idx2                <- which(len[idxT] <= d[j+1])
            idx                 <- intersect(idx1, idx2)
            nexp[idx,]          <- nei
            rownames(nexp)[idx] <- rownames(nei)
        }

        # Add complex expression values to original
        expr <- rbind(expr, nexp)         
    }

    return(expr)
}



## 
## Evaluate mRNA-miRNA interactions in subpathways
##

.indexSubpaths  <- function(subpaths, edgeList, type, org)
{
    if (type=='Linear')
    {
        # Replace every 2 consecutive genes with a corresponding interaction.
        # Create a named vector containing all pathway interactions.
        # Used to transform the entrez subpaths to indexed subpaths.
        # Index reference lexicon
        y   <- sprintf('%s-%s-%s', edgeList[,1], edgeList[,2], edgeList[,6]) 
        ids <- 1:length(y)
        names(ids) <- y

        # Split linear queue to batches of jobs
        k     <- floor(nrow(subpaths)/5)
        jobs  <- splitWork(N=nrow(subpaths), k)$jobs
        cores <- ifelse( nrow(subpaths) > detectCores()*1000, 'default', 1 )
        iSubs <- .doSafeParallel(
                    funcName=indexSubpath, export=NULL, 
                    combine='rbind', N=k, cores=cores,
                    subpaths, org, ids, jobs )
        return(iSubs)
    }

    if (type=='Non-Linear')
    {
        # Replace every 2 consecutive genes with a corresponding interaction.
        # Create a named vector containing all pathway interactions.
        # Used to transform the entrez subpaths to indexed subpaths.

        y   <- sprintf('%s-%s-%s', edgeList[,1], edgeList[,2], edgeList[,6]) 
        ids <- 1:length(y)
        names(ids) <- y
        rnames     <- gsub(org, '', rownames(subpaths))
        subpaths   <- apply(subpaths, 2, function(x, y) 
                                    {
                                        paste(x, y, sep='-') 
                                    }, rnames )
        rownames(subpaths) <- rnames        
        Ld                 <- matrix(ids[subpaths], ncol=ncol(subpaths), 
                                    byrow=FALSE)
        Ld[is.na(Ld)]      <- 0
        rownames(Ld)       <- rnames

        return(Ld)
    }
}


indexSubpath    <- function(i, ...)
{
    # ------- Unpacking arguments -------
    args     <- list(...)
    subpaths <- args[[1]]
    org      <- args[[2]]
    ids      <- args[[3]]
    jobs     <- args[[4]]
    # ----------------------------------
    idx      <- jobs[[i]]

    subs          <- subpaths[idx, ,drop=FALSE]
    x1            <- subs[, 1:(ncol(subs) - 1)]
    x2            <- subs[, 2:ncol(subs)]
    x3            <- gsub(org, '', rownames(subs))
    x             <- sprintf('%s-%s-%s', x1, x2, x3)
    x             <- x
    Ld            <- as.numeric(ids[x])
    Ld            <- matrix(Ld, ncol=ncol(subs)-1, byrow=FALSE)
    rownames(Ld)  <- rownames(subs)
    Ld[is.na(Ld)] <- 0

    return(Ld)
}

.scoreSubpaths  <- function(edgeScores, iSubs, type)
{
    if ( is.null(iSubs) || nrow(iSubs)==0 || is.null(edgeScores) ) 
    {
        return(NULL) 
    }

    # Recreate the original edgelist which corresponds to all interactions
    # and not just the ones in subpaths used to calculate edgescores,
    # since original subpath indexing was done using all interactions.
    # Otherwise, the subpaths would require reindexing.
    edgeScores  <- edgeScores[,-c(1:6)]
    maxIdx      <- max(as.integer(rownames(edgeScores)))
    E           <- matrix(, nrow=maxIdx, ncol=ncol(edgeScores))
    idx         <- as.integer(rownames(edgeScores))
    E[idx,]     <- as.matrix(edgeScores)
    # Zero as scores for each time point, for indexing of '0' to work
    edgeScores  <- rbind(E, rep(0, ncol(E))) 

    # Calculate subpath lengths
    len         <- getLengths(iSubs)
    t           <- ncol(edgeScores)
    
    # Replace '0' with a valid index which point to '0' again
    iSubs[iSubs == 0]  <- nrow(edgeScores)
    # Each column holds the score for each interaction from iSubs, 
    # columnwise (i.th time point).
    edgeScores         <- edgeScores[iSubs, ]

    subScores      <- matrix(0,nrow=nrow(iSubs), ncol=t)
    for (i in 1:t)
    {
        # iSubs turned into scores for i.th time point
        iSubsToScores     <- matrix(edgeScores[, i], ncol=ncol(iSubs))
        # Divide sum of scores with number of interactions to obtain subScore
        subScores[,i]     <- rowSums(iSubsToScores) / len
    }


    if (type == 'Non-Linear') { subScores <- normalize(subScores) }

    return(geneScores=subScores)    
}

.getEdgeScores  <- function(expr, edges, fci.C, fci.K, fci.T)
{
    if (is.null(edges) || is.null(expr))
    {
        return(NULL)
    }

    edgeList <- edges[,3:4]
    rel      <- edges[,5]
    pathway  <- edges[,6]

    # Calculate FCI score    
    fci <- foreach(i = 1:1, .export='scoreFCI',.combine='rbind') %do% 
    {
        scoreFCI(expr, edgeList, fci.C=fci.C, fci.K=fci.K, fci.T=fci.T)
    }
    t <- ncol(expr)

    # Calculate TVI score
    i   <- 0
    tvr <- foreach(i = 1:1, .export='scoreTVI', .combine='rbind') %do% 
    {
        scoreTVI(expr, edgeList)
    }

    # Calculate edge list's gene expressions
    x1 <- matrix(0,ncol=t, nrow=nrow(edgeList))
    x2 <- matrix(0,ncol=t, nrow=nrow(edgeList))
    for (j in 1:t)
    {
        x1[,j] <- expr[as.numeric(edgeList[,1]), j]
        x2[,j] <- expr[as.numeric(edgeList[,2]), j]    
    }

    # Calculate gene-gene interactions final score
    idx1       <- matrix(rel==1 & fci>0 & tvr>0, nrow=nrow(edgeList)) 
    idx2       <- matrix(rel==2 & fci<0 & tvr<0, nrow=nrow(edgeList))
    idx3       <- matrix(rel==3 & sign(fci*tvr)>0, nrow=nrow(edgeList))
    edgeScores <- (0.8*abs(fci) + 0.2*abs(tvr)) * (idx1 + idx2 + idx3)


    edgeScores           <- cbind(edges, edgeScores)
    rownames(edgeScores) <- rownames(edges)

    return(edgeScores)
}

.filterSubScore <- function(subpaths, scores, thr)
{
    if ( is.null(subpaths) || nrow(subpaths) == 0 || is.null(scores) )
    {
        return(NULL)
    }

    minNumberOfGoodScores <- 1

    # Calculate average of scores for all time points
    avgScorePerSub <- rowSums(scores)/ncol(scores)
    upaths         <- unique(rownames(subpaths))

    # Find highest ranking subpathway per pathway 
    idx0 <- c()
    for (i in 1:length(upaths))
    {
        # Find all subpaths of a single pathway
        lindx <- which(rownames(subpaths) == upaths[i])
        # Find the average scores for each of the subpaths
        lavg  <- avgScorePerSub[lindx]
        # Find max average value amomg subpaths
        mx    <- lindx[which.max(lavg)]

        if (avgScorePerSub[mx] > thr) 
        { 
            idx0 <- c(idx0, mx) 
        }
    }

    # Subpath with largest average score first
    idx0  <- idx0[order(avgScorePerSub[idx0], decreasing = TRUE)]
    
    # Keep subpaths with at least one good score
    rSums <- rowSums(matrix(scores > thr, nrow=nrow(scores)))
    
    # Remove top results of first phase
    iiidx <- which(rSums >= minNumberOfGoodScores)
    idx1  <- setdiff(iiidx, idx0)

    # Sort by decreasing average
    idx1  <- idx1[order(avgScorePerSub[idx1], decreasing = TRUE)]
    
    # Join previous two results categories
    idx   <- c(idx0, idx1)
    
    return(idx)
}

