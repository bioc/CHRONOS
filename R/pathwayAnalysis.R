##
## Genelist and edgelist maintainance
##


matrixToLexicon <- function( adjMats, org, groupMode )
{
    if ( missing(groupMode) )
    {
        groupMode <- 'collapse'
    }

    lex <- list('org'=org,'genes'=NULL,'edges'=NULL, 'adjMat'=NULL)
    if (length(adjMats) == 0) { return(lex) }

    #
    # Find unique genes in organism's pathways
    #
    genes <- c()
    for (i in 1:length(adjMats))
    {
        adjMat  <- adjMats[[i]]
        res     <- c() 

        for (m in 1:length(rownames(adjMat)))
        {
            s   <- rownames(adjMat)[m]
            s   <- gsub(paste(org,':',sep=''), '', s )
            res <- c(res, s)
        }
        genes <- c(genes, unique(unlist(res)) )
    }
    if (groupMode == 'expand') { genes <- sort(unique(as.numeric(genes))) }
    
    geneLex <- NULL
    if ( length(genes) > 0 )
    {
        geneLex         <- 1:length(genes)
        names(geneLex)  <- genes
    }

    #
    # Find unique edges in organism's pathways
    #
    edges   <- data.frame()
    pathway <- names(adjMats)
    for (i in 1:length(adjMats))
    {
        if (!is.matrix(adjMats[[i]])) { next() }
        if (nrow(adjMats[[i]]) == 0) { next() }

        X     <- t(adjMats[[i]])
        rn    <- gsub(paste(org,':',sep=''), '', rownames(X))
        cn    <- gsub(paste(org,':',sep=''), '', colnames(X))
        y     <- expand.grid(rn, cn, stringsAsFactors = FALSE)
        y     <- cbind(y, as.vector(X), gsub(org,'',pathway[i]))
        edges <- rbind(edges, y[which(X > 0),])
    }
    edges   <- unique(edges)
    edgeLex <- NULL
    if ( nrow(edges) > 0 )
    {
        # Indexed columns of edgelist, with entrez ids as names.
        edgeLex <- data.frame( 'E1'=edges[,2],'E2'=edges[,1],
                        'e1'=geneLex[ as.character(edges[,2]) ],
                        'e2'=geneLex[ as.character(edges[,1]) ], 
                        edgeType=edges[,3], pathway=edges[,4], 
                        stringsAsFactors = FALSE )
    }

    adjMat <- NULL
    if (nrow(edges) > 0)
    {
        # Keep only one interaction type
        idx     <- which(!duplicated(edgeLex[,1:2]))
        adjMat  <- ftM2adjM(as.matrix(edgeLex[idx, 1:2]), edgeLex[idx, 5], 
                            edgemode="directed")
        adjMat  <- adjMat[order(as.numeric(rownames(adjMat)) ), ]
        adjMat  <- adjMat[, order(as.numeric(colnames(adjMat)) )]
    }

    # Return values
    lex <- list('org'=org, 'genes'=geneLex, 'edges'=edgeLex, 'adjMat'=adjMat)
    
    return(lex)
}


filterMatrix <- function( adjMats, org, userGenes )
{

    if (is.null(userGenes))
    {
        return(adjMats)
    }

    filtAdjMats <- adjMats
    for (i in 1:length(adjMats))
    {
        if (is.null(adjMats[[i]]))  { next() }

        genes  <- rownames(adjMats[[i]])
        genes  <- lapply(genes, function(x) 
                    { matrix(as.numeric(unlist(strsplit(x, ' '))), nrow=1) })
        genes       <- unlistToMatrix(fillMatrixList(genes))
        idxGenes    <- matrix(genes %in% userGenes, nrow=nrow(genes)) * 1
        filtGenes   <- genes * idxGenes
        filtGenes[filtGenes == 0] <- ''
        rnames      <- apply(filtGenes, 1, function(x) 
                    { paste(x[which(x != '')], collapse=' ') })

        rownames(filtAdjMats[[i]])  <- rnames
        colnames(filtAdjMats[[i]])  <- rnames
        idx                         <- which(rnames != '')
        filtAdjMats[[i]]            <- filtAdjMats[[i]][idx, idx, drop=FALSE]

        # Remove rows and columns with duplicated names        
        idx              <- which(!duplicated(rownames(filtAdjMats[[i]])))
        filtAdjMats[[i]] <- filtAdjMats[[i]][idx, idx, drop=FALSE]
    }

    return(filtAdjMats)
}


##
## MiRNAs
##
.createMirnaLexicon <- function(org, genes)
{
    # Given a list of genes, and a file which maps miRNA's with their gene 
    # targets. Construct a vector containing the genes, followed by the 
    # miRNA's that target them. Their indexes are used to construct a two 
    # column matrix, containing the mirna's on the first column, and their 
    # respective gene targets on the second (index-index). 
    
    # Import mirgene map
    targets     <- importMiRNAFile(org)
    
    # Keep user miRNAs 
    e=new.env(); load(file=cache$dirs$mirnaExpressions, envir=e)
    miEx        <- e[[ls(e)[1]]]
    targets     <- targets[which((targets[,1] %in% rownames(miEx))),]
    
    # Find mirnas that target the given genes
    targets     <- targets[which((targets[,2] %in% names(genes))),]
    userMirs    <- unique(targets[,1])

    # Mirna names are accompanied by their indexes 
    mirs        <- (length(genes) + 1) : (length(genes) + length(userMirs))
    names(mirs) <- userMirs

    # Create edge list between mirnas and gene targets 
    edgeList           <- unique(targets[, 1:2])
    rownames(edgeList) <- 1:nrow(edgeList)

    # Create indexer
    interact     <- c(names(genes), userMirs)
    lib          <- 1:length(interact) 
    names(lib)   <- as.character(interact)

    # Edgelist
    idx           <- matrix(0, nrow=nrow(edgeList), ncol=2)
    # Index mirna's column
    idx[,1]       <- as.numeric(lib[edgeList[,1]])
    # Index gene targets column
    idx[,2]       <- as.numeric(lib[as.character(edgeList[,2])])
    colnames(idx) <- c(' ', ' ')
    edgeList      <- cbind(idx, edgeList)

    return(list('nodes'=mirs, 'edges'=edgeList))
}

    

##
## Others
##
getSubpathwayGenes <- function(subpaths, type)
{
    if ( !is(subpaths,'matrix') )
    {
        subpaths <- matrix(subpaths, nrow=1)
    }

    if (type == 'Linear')     
    {
        genes <- unique(as.vector(subpaths))
        genes <- genes[which(genes!=0)]
    }
    if (type == 'Non-Linear') 
    {
        usubs <- vector(mode='list', length=nrow(subpaths))
        for(i in 1:nrow(subpaths))
        {
            usub       <- unique(unlist(strsplit(subpaths[i,], split='-')))
            usubs[[i]] <- matrix(as.numeric(usub[which(usub!='0')]), nrow=1)
        }

        # Join list of named matrices to a single named matrix
        usubs <- unlistToMatrix(fillMatrixList(usubs))
        genes <- unique(as.vector(usubs))
    }


    return(genes)
}
