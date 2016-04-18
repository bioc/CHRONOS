#
# Fold Change Interactivity (FCI) Score
#


scoreFCI <- function(expr, edgeList, fci.C, fci.K, fci.T)
{
    # Captures pairs of nodes with similar high positive, negative, 
    # or strongly opposed fold change values .

    if (missing(fci.C)) { fci.C <- 1 }
    if (missing(fci.K)) { fci.K <- 5 }
    if (missing(fci.T)) { fci.T <- 0.3 }


    fciScore <- array(0, dim=c(nrow(edgeList), ncol(expr)))
    for (e in 1:nrow(edgeList))
    {
        for (t in 1:ncol(expr))
        {
            i <- expr[edgeList[e,1], t] # Indexing by integer
            j <- expr[edgeList[e,2], t]

            s <- 1
            if (i*j < 0)
            {
                i <- abs(i)
                j <- abs(j)
                s <- -1;
            }

            s1 <- 1 + fci.C*(exp(-fci.K*(i-fci.T)) + exp(-fci.K*(j-fci.T)))
            s2 <- 1 + fci.C*(exp(-fci.K*(-i-fci.T)) + exp(-fci.K*(-j-fci.T)))
            fscore <- s * abs(s1^(-1) - s2^(-1))

            fciScore[e,t] <- fscore
        }
    }

    return(fciScore)
}


#
# Time Varying Interactivity (TVI) Score
#
scoreTVI  <- function(expressions, edgeList)
{
    #
    # Analyzes the interactivity dynamics among nodes by 
    # considering the dynamics of their expression changes.
    #
    wRange     <- 1
    weights    <- seq(-wRange,wRange)
    w          <- length(weights)
    expressNrm <- expressions

    nE         <- nrow(edgeList)
    prob_W     <- array(rep(0,w*1*nE), dim=c(w,1,nE))


    for (iq in 1:nrow(edgeList))
    {
        a1              <- expressNrm[edgeList[iq,1], 1]
        a2              <- expressNrm[edgeList[iq,2], 1]
        a               <- a1 * a2 
        prob_W[, 1, iq] <- exp(weights*a) / sum(exp(weights*a))
    }

    tv <- .markov(expressNrm, edgeList, weights, prob_W)  

    return(tv)
}

.markov   <- function(X, edges_mat, W, prob_w)
{
    # This function computes the Forward Backward iterates for the 
    # specified edge, E=(i,j), given observed data, X, and estimate for 
    # transition probability
    nE        <- nrow(edges_mat)
    nW        <- length(W)
    var_T     <- ncol(X)
    # memory allocation
    var_O     <- array(rep(0,nW*var_T*nE), dim=c(nW,var_T,nE))
    var_F     <- array(rep(0,nW*var_T*nE), dim=c(nW,var_T,nE))
    # pre-compute
    W             <- t(as.matrix(W))
    for (i in 1:var_T)
    {
        Xmult         <- as.matrix(X[edges_mat[,1],i] * X[edges_mat[,2], i])
        Xmult         <- as.matrix(apply(Xmult, 1, function(x) sum(x) ))
        Wmult         <- exp(t(W) %*% t(Xmult))  ## -t?
        xx             <- t(as.matrix(apply(Wmult, 2, function(x) sum(x) )))
        Wmult          <- Wmult/kronecker(matrix(1,nW,1), xx)
        var_O[,i,]     <- matrix(as.vector(Wmult), nrow=nW, ncol=nE) 
    }

    # initialization
    var_F[,1,]  <- prob_w

    for (t in 2:var_T)
    {
        for (e in 1:nE)
        {
            var_F[, t, e] <- var_O[, t, e] * var_F[, t-1, e]
            var_F[, t, e] <- var_F[, t, e] / sum( var_F[, t, e] )
        }
    }

    pW     <- var_F
    cW     <- matrix(0, nE, var_T)
    for (i in 1:var_T)
    {
        r      <- apply(pW[,i,] , 2, function(x) which.max(x))
        cW[,i] <- W[r]
    }

    return(cW)
}



#
# Evaluate miRNA-mRNA interactions
#
.getMiScores    <- function(geneExpr, mirExpr, edges, lex, mode, fci.C, fci.K,
                            fci.T)
{
    # Create joint expression matrix
    expressions     <- rbind(geneExpr, mirExpr) 
    edgeList        <- edges[,1:2]

    k <- 1
    if (nrow(edgeList) > 3000) 
    { 
        k <- floor(nrow(edgeList)/3000) 
    }
    S  <- splitWork(N=nrow(edgeList), k)

    i  <- 0
    export <- c('scoreTVI', 'scoreFCI')
    IS <- foreach (i = 1:S$k, .export=export, .combine='rbind') %do%
    {
        # Calculate TVI score
        tvr <- scoreTVI(expressions, edgeList[S$jobs[[i]], , drop=FALSE])

        # Calculate FCI score
        fci <- scoreFCI(expressions, edgeList[S$jobs[[i]], , drop=FALSE], 
                            fci.C=fci.C, fci.K=fci.K, fci.T=fci.T) 

        # Calculate mir-gene interactions final score
        # len     <- nrow(edgeList[S$jobs[[i]], , drop=FALSE])
        idx     <- matrix(sign(fci*tvr)>0, nrow=length(S$jobs[[i]]))
        score   <- (0.8*abs(fci) + 0.2*abs(tvr))*idx

        return(score)
    }

    return(edgeScore=IS)
}

.filterMiScores <- function(subpaths, mlex, mirnaScoreThres, type, lex)
{
    
    scoreThr     <- mirnaScoreThres

    if (is.null(subpaths) || nrow(subpaths) < 1) { return(subpaths) }

    # Keep interaction with at least one good score in any time point
    E           <- mlex$edges
    miScore     <- E[, 5:ncol(E)]
    t           <- ncol(miScore)
    idx         <- which(miScore > scoreThr)
    M           <- matrix(0, nrow=nrow(miScore), ncol=ncol(miScore)) 
    M[idx]      <- 1
    midx        <- which(rowSums(M) > 0)
    names(midx) <- midx
    E           <- E[midx, 3:ncol(E)]

    # Find miRNA-gene interactions relevant to subpathway genes.
    genes       <- getSubpathwayGenes(subpaths, type)
    E           <- E[which(E[,2] %in% genes), ]
    Ns          <- getLengths(subpaths)

    # Find unique genes in each subpathway    
    subGenes <- vector(mode='list', length=nrow(subpaths))
    for (i in 1:nrow(subpaths))
    {
        subGenes[[i]] <- getSubpathwayGenes(subpaths[i,, drop=FALSE], type)
    }

    # Find miRNAs targeting each gene
    uGenes  <- unique(E[,2])
    R       <- vector(mode='list', length=length(uGenes))
    for (i in 1:length(R))
    {
        idx     <- which(E[,2] %in% uGenes[i])
        R[[i]]  <- idx
    }
    names(R) <- uGenes


    # Find gene targets for each miRNA
    uniqueMiRNAs    <- unique(E[,1])
    miRNATargets    <- vector(mode='list', length=length(uniqueMiRNAs))
    for (i in 1:length(uniqueMiRNAs))
    {
        miRNATargets[[i]] <- E[which(E[,1] == uniqueMiRNAs[i]), 2]
    }
    names(miRNATargets) <- uniqueMiRNAs

    k <- 1
    if (nrow(subpaths) > 5) 
    {  
        k <- floor(nrow(subpaths)/5) 
    }
    jobs   <- splitWork(N=nrow(subpaths), k)$jobs
    subIds <- 1:nrow(subpaths)

    # Extract compact linear subpaths
    export <- c('getSubpathwayGenes', 'phyper', 'na.omit')
    cores  <- ifelse( nrow(subpaths) > detectCores()*1000, 'default', 1 )

    r <- .doSafeParallel( 
            funcName=miRNAsubScoring, export=export,
            combine='c', N=k, cores=cores,
            subpaths, E, R, genes, Ns, scoreThr, type, subIds, miRNATargets,
            jobs)


    # Unlist first level of results
    subMirs            <- do.call('c', r[seq(1, length(r), 3)])
    mirScores          <- do.call('rbind', r[seq(2, length(r), 3)])
    edgeList           <- do.call('rbind', r[seq(3, length(r), 3)])
    subMirs            <- matrix(subMirs, ncol=1)
    rownames(subMirs)  <- 1:nrow(subMirs)
    colnames(subMirs)  <- 'miRNAsOverSubpathway'
    rownames(edgeList) <- 1:nrow(edgeList)

    # Filter multiple edges
    idx <- which(duplicated(edgeList[, 1:2]))
    if ( length(idx) > 0 )
    {
        edgeList <- edgeList[-idx, ]
    }
    # Sort edgeList according to miRNA
    idx <- order(as.numeric(gsub("\\D", "", edgeList[, 1])))
    edgeList <- edgeList[idx,  ]
    rownames(edgeList) <- 1:nrow(edgeList)

    # Remove unused rows
    mirScores <- mirScores[rowSums(is.na(mirScores)) != ncol(mirScores), ]
    colnames(mirScores) <- c('miRNA', 'Entrez', paste0('T', 1:t))

    res <- list(miRNAsOverSubpathway=subMirs, mirScores=mirScores, 
                edgeList=edgeList)

    return(res)
}

miRNAsubScoring <- function(i, ...)
{

    # ------- Unpacking arguments -------
    args         <- list(...)
    subpaths     <- args[[1]]
    E            <- args[[2]]
    R            <- args[[3]]
    genes        <- args[[4]]
    Ns           <- args[[5]]
    scoreThr     <- args[[6]]
    type         <- args[[7]]
    subIds       <- args[[8]]
    miRNATargets <- args[[9]]
    jobs         <- args[[10]]
    # ----------------------------------
    idx      <- jobs[[i]]
    subpaths <- subpaths[idx, ,drop=FALSE]
    Ns       <- Ns[idx]
    subIds   <- subIds[idx]
    # ----------------------------------


    t           <- ncol(E) - 2 
    subMirs     <- rep('', nrow(subpaths))
    mx          <- nrow(subpaths)*length(unique(E[,1]))
    mirScores   <- matrix(, ncol=(2+t), nrow=mx)

    finalInteractions <- c()
    cRow              <- 1
    for (si in 1:nrow(subpaths)) 
    {
        # Collapse the miRNA-gene scoring matrix to miRNAs targeting the sub.
        locGenes <- getSubpathwayGenes(subpaths[si,], type)
        map <- E[unique(unname(unlist(R[as.character(locGenes)]))), ]

        # Find the miRNAs that target any of the subpathway
        miOverSub   <- unique(map[,1])

        if (length(miOverSub) == 0) { next() }

        # Find statistical significance of each miRNA targeting the sub.
        sigMiOverSub <- miOverSub
        for (j in 1:length(miOverSub))
        {
            # Genes in subpath targeted by the miRNA.
            A       <- length(unique(map[, 2]))
            # User genes
            B       <- length(genes)    
            # Number of gene targets of miRNA in miRecords.
            C       <- length(unlist(miRNATargets[miOverSub[j]]))
            # Subpath length
            D       <- Ns[si]

            if ( (1 - phyper(A, C, B-C, D)) > 0.05) { sigMiOverSub[j] <- NA }
        }

        sigMiOverSub                         <- na.omit(sigMiOverSub)
        attributes(sigMiOverSub)$'na.action' <- NULL

        # Calculate miRNA-subpathway scores
        subMsg              <- NULL
        mirScore            <- matrix(, nrow=length(sigMiOverSub), ncol=(2+t))
        colnames(mirScore)  <- c('miRNA', 'sub', paste('T', 1:t, sep=''))

        for (j in 1:length(sigMiOverSub))
        {
            mir <- sigMiOverSub[j]
            e   <- map[which(map[,1] == mir), , drop=FALSE]

            if (nrow(e) > 0)
            {                
                fScores             <- e[,-c(1,2), drop=FALSE]
                colnames(fScores)   <- paste('T', 1:ncol(fScores), sep='')
            }

            r   <- colMeans(e[-c(1,2)],)
            idx <- unname(which(r >= scoreThr))

            if (length(idx) > 0)
            {
                tP     <- paste(paste('t', idx, sep=''), collapse=',')
                msg    <- sprintf('%s(%s)',mir, tP)
                subMsg <- c(subMsg, msg)

                mirScore[j, 1]         <- sigMiOverSub[j]
                mirScore[j, 2]         <- subIds[si]
                mirScore[j, 3:(2 + t)] <- round(r*10^4)/10^4
            }

            # Save edgelist indexes that pass the bar
            if ( nrow(e) > 0 )
            {
                for ( z in 1:nrow(e) )
                {
                    idx <- unname(which(e[z, , drop=FALSE] >= scoreThr))
                    if (length(idx) > 0)
                    {
                        rname <- rownames(e[, z, drop=FALSE])
                        finalInteractions <- c(finalInteractions, rname)
                    }
                }
            }
        }
        if (nrow(mirScore) > 0)
        {
            idx                 <- cRow : (cRow + nrow(mirScore) - 1)
            mirScores[idx, ]    <- mirScore
            cRow                <- cRow + nrow(mirScore)
        }

        subMirs[si] <- paste(subMsg,  collapse=' ')
    }

    # Return edgeList containing only strong miRNA-mRNA interactions
    edgeList <- E
    idx <- which(rownames(E) %in% unique(finalInteractions))
    if ( length(idx) > 0 )
    {
        edgeList <- E[idx, , drop=FALSE]
    }

    return(list(subMirs=subMirs, mirScores=mirScores, edgeList=edgeList))
}


