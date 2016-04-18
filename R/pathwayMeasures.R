

pathwayMeasures       <- function(graphs)
{
    message('Topological evaluation of genes in pathways...',appendLF = FALSE)
    
    rnd     <- 10^4
    adjMats <- graphs$expanded
    org     <- graphs$org
    f       <- vector(mode='list', length=3)
    lex     <- matrixToLexicon( adjMats=adjMats, org=org, groupMode='expand')

    # Topological evaluation of genes in pathways 1 - Bridgeness
    res1 <- .bridgeness(adjMats, org)
    # Topological evaluation of genes in pathways 2 - Betweenness Centrality
    res2 <- .betweennessCentrality(lex)
    # Topological evaluation of genes in pathways 3 - Degree
    res3 <- .degree(lex)

    res           <- matrix(,ncol=3, nrow=length(lex$genes))
    colnames(res) <- c('pathness', 'BC', 'DEG')
    rownames(res) <- names(lex$genes)
    res[,1]       <- round(normalize(res1) * rnd)/rnd
    res[,2]       <- round(normalize(res2) * rnd)/rnd
    res[,3]       <- round(normalize(res3) * rnd)/rnd
    
    message('done.')

    return(res)
}



.bridgeness            <- function(adjMats, org)
{

    .pathwayGenes          <- function(adjMats, org, lex)
    {
        # Create a matrix which contains the unique genes of a pathway in each 
        # column. Columns contain the organism's pathways. 

        N    <- length(adjMats)
        ids  <- vector(mode='list', length=N)
        lens <- vector(mode='numeric', length=N)
        i   <- 1
        for (adjMat in adjMats)
        {
            if(!is.matrix(adjMat)) { i <- i+1; next() }

            locIds <- c() 
            for (m in 1:length(rownames(adjMat)))
            {
                s      <- rownames(adjMat)[m]
                locIds <- c(locIds, gsub(paste(org,':',sep=''),'',s ))     
            }
            # Filter genes keeping only the ones that appear in the subpaths.
            if (!missing(lex)) 
            { 
                locIds <- intersect(locIds, names(lex$genes)) 
            }
            ids[[i]] <- as.numeric(locIds)
            lens[i]  <- length(locIds)
            i        <- i + 1
        }

        # Turn list of uneven subpaths to matrix
        M <- matrix(0, nrow=max(lens), ncol=N)
        for (i in 1:N)
        {
            if (lens[i] > 0) { M[1:lens[i],i] <- ids[[i]] }
        }

        return(M)
    }


    # Calculates the bridgeness of all provided genes based on 
    # their intersections betweeen pathways.
    
    cat('Calculating bridgeness of pathway genes...')
    
    modules <- .pathwayGenes(adjMats, org)
    m       <- ncol(modules)
    nodes   <- sort(unique(as.vector(modules)))
    if (length(nodes) < 1) { return(NULL) }

    if (nodes[1] == 0) { nodes <- nodes[2:length(nodes)] }

    b      <- vector(mode='numeric', length=length(nodes))

    for (i in 1:(m-1))
    {
        for (j in (i+1):m)
        {
            common <- intersect(modules[,i], modules[,j])
            idx    <- nodes %in% common
            b[idx] <- b[idx] + 1/length(common)
        }
    }
    res        <- b
    names(res) <- nodes

    cat('done\n')

    return(res)
}


.betweennessCentrality <- function(lexicon)
{
    # Calculates the betweeness centrality for all provided genes based on 
    # their corresponding iteractions.

    cat('Calculating betweeness centrality of pathway genes...')

    genes <- lexicon$genes
    edges <- lexicon$edges

    # Create the adjacency matrix for the whole organism
    N      <- length(genes)
    adjMat <- matrix(0, nrow=N, ncol=N)
    for (i in 1:nrow(edges))
    {
        adjMat[edges[i,3], edges[i,4]] <- 1
    }
    rownames(adjMat) <- names(genes)
    colnames(adjMat) <- names(genes)

    # Create a graph based on the adjacency matrix
    g   <- ugraph(as(graphAM(adjMat, 'directed'),'graphNEL'))

    # Calculate the betweeness centrality for each gene
    r   <- brandes.betweenness.centrality(g)
    bc  <- t(r$betweenness.centrality.vertices)
    res <- unname(as.vector(bc))
    names(res) <- names(genes)

    cat('done\n')

    return(res)
}


.degree                <- function(lexicon)
{
    cat('Calculating degree of pathway genes...')

    genes <- lexicon$genes
    edges <- lexicon$edges

    # Create the adjacency matrix for the whole organism
    N      <- length(genes)
    adjMat <- matrix(0, nrow=N, ncol=N)
    for (i in 1:nrow(edges))
    {
        adjMat[edges[i,3], edges[i,4]] <- 1
    }

    res <- vector(mode='numeric', length=length(genes))
    for (i in 1:nrow(adjMat))
    {
        res[i] <- sum(adjMat[i,]) + sum(adjMat[,i])
    }

    cat('done\n')

    return(res)
}

