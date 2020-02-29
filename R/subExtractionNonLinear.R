

extractNonLinearSubpathways  <- function(graphs, pathways, a, b, k, filter, 
                                        groupMode, export, verbose)
{
    if ( missing(verbose) ) { verbose <- TRUE }
    
    tryCatch(
    {        
        message('Extracting Non Linear Subpathways...', appendLF = FALSE)
        
        subsNonLinear <- .extractNonLinearSubpathways(graphs, pathways, 
                        a, b, k, filter, groupMode, export, verbose)

        message('done.')

        return(subsNonLinear)
    }, warning = function(war) 
    {
    }, error = function(e) 
                { 
                    closeAllConnections()
                    message('\nError in non linear extraction.'); 
                }
    , finally = { }
    )   
}


.extractNonLinearSubpathways <- function(graphs, pathways, a, b, k, filter, 
                                        groupMode, export, verbose)
{
    if (missing(verbose))   { verbose   <- TRUE }
    if (missing(a))         { a         <- 3  }
    if (missing(b))         { b         <- 10 }
    if (missing(k))         { k         <- 2  }
    if (missing(filter))    { filter    <- FALSE }
    if (missing(groupMode)) { groupMode <- c('expand') }
    if (missing(export))    { export    <- c('.RData') }

    opts     <- list(filter=filter, k=k, a=a, b=b, export=export)
    dirs     <- cache$dirs
    eAdjMats <- graphs$expanded
    cAdjMats <- graphs$combined
    org      <- graphs$org

    if (!missing(pathways))
    {
        # Keep a subset of pathways
        pathways <- paste(org, pathways, sep='')
        eAdjMats <- eAdjMats[pathways]
        cAdjMats <- cAdjMats[pathways]
    }

    .createDirectories(org)

    # Find user genes
    if (opts$filter)   
    {
        e=new.env(); load(cache$dirs$geneExpressions, e)
        geneEx     <- e[[ls(e)[1]]]
        genes     <- sort(as.numeric(unique(rownames(geneEx))))
    }
    if (!(opts$filter)) { genes <- NULL}


    # Filter adjacency matrices according to user genes
    if (groupMode == 'expand')   { adjMats <- eAdjMats }
    if (groupMode == 'collapse') { adjMats <- cAdjMats }

    adjMats  <- filterMatrix(adjMats=adjMats, org=org, userGenes=genes)

    # Build global gene and edge indexer
    lex      <- matrixToLexicon(adjMats=adjMats, org=org, 
                                    groupMode=groupMode)

    # Convert edgeType matrix to adjacency matrix
    adjMatsBin <- adjMats
    for (i in 1:length(adjMatsBin) ) 
    {
        if (is.matrix(adjMatsBin[[i]])) #if (!is.null(adjMatsBin[[i]]))
        {
            adjMatsBin[[i]][adjMatsBin[[i]] > 1] <- 1 
        }         
    }

    # Apply k-clique algorithm to find the k-subpaths for each pathway.
    export <- c('unlistToMatrix', 'as', 'graphAM', 'kCliques', 'ugraph')
    cores  <- ifelse( length(adjMatsBin) > detectCores()*10, 'default', 1 )
    subs   <- .doSafeParallel(funcName=.nlSubpaths,
                    export=export,
                    combine='default',
                    N=length(adjMatsBin),
                    cores=cores,
                    adjMatsBin, lex, k, names(adjMatsBin), a, b)

    if ( !is.null(subs) )        { names(subs) <- names(adjMatsBin) }
    if ( !is(subs,'list') ) { subs <- list(subs) }

    
    # Ensure that all subpaths have the same number of columns
    subs <- fillMatrixList(subs)

    # Join list of named matrices to a single named matrix
    subs <- unlistToMatrix(subs)

    # Delete temporaty files and directories
    .cleanDirectories()

    # Export subpaths
    exportSubpaths(subs, dir=paste(dirs$nlr, org, sep='//'), type=opts$export,
                    verbose=FALSE)

    subsNonLinear <- list(subpaths=subs, adjMats=adjMats, 
        groupMode=groupMode, org=org, lex=lex, subpathwayType='Non-Linear')

    return(subsNonLinear)
}


.nlSubpaths                   <- function(i, ...)
{
    # ------- Unpacking arguments -------
    args      <- list(...)
    adjMat    <- args[[1]][[i]]
    lex       <- args[[2]]
    k         <- args[[3]]
    pathname  <- args[[4]][i]
    a         <- args[[5]]
    b         <- args[[6]]
    # ----------------------------------

    org      <- lex$org
    edgeList <- lex$edges

    if (!is.matrix(adjMat)) { return() }

    tmp              <- paste(org,':',sep='')
    rownames(adjMat) <- gsub(tmp,'', rownames(adjMat))  #!!
    colnames(adjMat) <- gsub(tmp,'', colnames(adjMat))  #!!

    # Create a graphNEL graph from the adjacency matrix
    if (nrow(adjMat) * sum(adjMat) > 1e05)  
    { adjMat[lower.tri(adjMat, diag = FALSE)] <- 0 }
    G                <- as(graphAM(adjMat, 'directed'), 'graphNEL')
    ksubpath         <- kCliques(ugraph(G))[[paste(k, 'cliques', sep='-')]]
    if (is.null(ksubpath)) { return() }

    # Keep cliques of length a to b
    lgt             <- sapply(ksubpath, function(x) { length(x) })

    # Set max column count if right border is not supplied
    if(is.null(b)) { b <- max(lgt) }
    droppedSubsIdx  <- which(lgt<a  | lgt>b)
    acceptedSubsIdx <- which(lgt>=a & lgt<=b)

    # A subpath is stored in each line on Ld, and its maximum length is b.
    Ld <- matrix(0, nrow=length(ksubpath), ncol=b)

    for (j in acceptedSubsIdx)
    {
        Ld[j,] <- c(ksubpath[[j]], rep(0, b-lgt[j]))
    }
    if (length(droppedSubsIdx) > 0) 
    {  
        Ld <- matrix(Ld[-droppedSubsIdx,], ncol=b) 
    }
    if (nrow(Ld) == 0) {  return() }


    # Each subpath by now is a bucket of genes. Lets consider all the possible
    # interactions between those genes. Not all of them correspond to actual 
    # interactions. The ones that do are the final subpath. 
    R <- list()

    for (j in seq(nrow(Ld),1))
    {
        # Edge list indexes
        r1   <- as.numeric(is.element(edgeList[,1], Ld[j,]))
        r2   <- as.numeric(is.element(edgeList[,2], Ld[j,]))
        idx1 <- which((r1 + r2) == 2)
        idx2 <- which(edgeList[,6] == gsub(org, '', pathname))
        r    <- intersect(idx1, idx2)

        if ( length(r) > 0 ) 
        { 
            R <- c(R, list(r))
        }
    }
    names(R) <- rep(pathname, length(R))

    # Turn list of uneven subpaths to matrix (indexes of edgelist)
    iSubpaths <- unlistToMatrix(R)

    # Replace zeros with next availiable row index (max_index + 1). 
    iSubpaths[which(iSubpaths==0)] <- nrow(lex$edges) + 1

    # Create mapper from indexes to entrez numbers
    lib <- paste(lex$edges[,1], lex$edges[,2], sep='-')

    # Zero mapping (map (max_index + 1) to zero)
    lib <- c(lib, 0)

    # Turn indexes to entrez names edges
    subpaths <- matrix(lib[iSubpaths], nrow=nrow(iSubpaths))

    return(subpaths)
}


