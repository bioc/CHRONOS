#
# Evaluate statistical significance
#


.getStatistics     <- function(subpaths, genes, type, pVal, pAdjust)
{
    # The statistical significance of the subscore is evaluated via
    # cumulative hypergeometric distribution so as to examine the enrichment
    # of a subpathway in any user-defined set of interesting genes 
    # (e.g. list of differentially expressed genes (DEGs)). 
    #

    # d: Interesting genes involved in the subpathway
    # D: User defined set of interesting genes in pathway list
    # G: User defined genes in pathway list 
    # l: All genes involved in the subpathway

    message('Evaluating statistical significance of subscore...', 
        appendLF = FALSE)
    if (missing(pVal))    { pVal <- 0.05 }
    if (missing(pAdjust)) { pAdjust <- '' }

    # Read interesting genes from directory
    inFile <- cache$dirs$interestingGenes
    if (!file.exists(inFile)) 
    { 
        stop('No interesting genes provided\n') 
    }
    
    interGenes        <- as.integer(readLines(inFile))
    subGenes          <- getSubpathwayGenes(subpaths, type)
    interGenesInSubs  <- interGenes[interGenes %in% subGenes] 
    interGenesInPaths <- interGenes[interGenes %in% genes]

    if ('Linear' %in% type)
    {
        # Genes of each subpathway considered as interesting
        x     <- rowSums(matrix(subpaths %in% interGenesInSubs, 
                        nrow=nrow(subpaths)))
        # Interesting genes in all pathways
        D     <- length(interGenesInPaths)
        # All user genes (in pathways)
        G     <- length(genes)
        # Number of genes in each subpathway
        l     <- getLengths(subpaths)
        # P-value calculation
        PV <- 1 - phyper(x, D, G-D, l)
    }

    if ('Non-Linear' %in% type)
    {
        PV <- rep(0, nrow(subpaths))
        for(i in 1:nrow(subpaths))
        {
            subGenes     <- getSubpathwayGenes(subpaths[i,], type)
            # Genes of each subpathway considered as interesting
            nH             <- length(which(subGenes %in% interGenesInSubs))
            # Interesting genes in all pathways
            nD             <- length(interGenesInPaths)
            # All user genes
            nG             <- length(genes)
            # Number of genes in each subpathway
            nL             <- length(subGenes)
            # P-value calculation
            PV[i]         <- 1 - phyper(nH, nD, nG-nD, nL)
        }
    }

    # The false positive discovery rate of P results is reduced 
    # by providing False Discovery Rate (FDR) corrected P-values (q-values)
    if (pAdjust == 'fdr') { QV <- p.adjust(PV, method='fdr') }
    if (pAdjust != 'fdr') { QV <- PV }

    message('done.')

    return(qValues=QV)
}



#
# Evaluate topological relations
#
.getMeasures       <- function(genes, subpaths, measures, mode)
{
    message('Evaluating subpaths based on gene topology...',appendLF = FALSE)
    if (is.null(subpaths) || nrow(subpaths)==0)
    { 
        message('done.'); return(NULL) 
    }

    T1 <- measures[,1]
    T2 <- measures[,2]
    T3 <- measures[,3]

    # Calculate the bridgeness and the betweenness centrality for all 
    # organism genes. Returns a list of two vectors. One named vector with 
    # the bridgeness values, and one with the betweenness centrality values. 
    # The order of these values is by increasing entrez gene id. 
    # subpathsGenes matrix contains the unique genes of each pathway.

    if (mode == 'Linear')
    {
        # Linear subpaths consist of consecutive genes. To index them, 
        # we assign each gene to real number and use this index to access 
        # the topological values vectors, which contain values in ascending 
        # order of entrez gene id.

        subpathsGenes <- subpaths
    }
    
    if (mode == 'Non-Linear')
    {
        # Non Linear subpaths are a set of interactions. Since ee are 
        # interested in the genes that interact and not the interactions 
        # themselves, we keep the set of unique genes for each subpath and 
        # score as shown previously.

        subpathsGenes <- matrix(0, nrow=nrow(subpaths), ncol=ncol(subpaths))
        for (i in 1:nrow(subpaths))
        {
            uGenes <- unique(unlist(strsplit(subpaths[i,], split='-')))
            subpathsGenes[i,1:length(uGenes)] <- uGenes
        }
    }


    # Evaluate subpaths with respect to the genes they consist of, each 
    # of which has its own bridgeness and betweenness centrality score.
    y          <- as.character(genes)
    y          <- c(y, 0)
    ids        <- 1:(length(y)+1)
    names(ids) <- y
    iSubs      <- matrix(ids[as.character(subpathsGenes)], 
                        ncol=ncol(subpathsGenes), byrow=FALSE)

    # Create indexing structures
    lib1 <- c(T1, 0)
    lib2 <- c(T2, 0)
    lib3 <- c(T3, 0)

    # Subpath length per line
    L  <- getLengths(iSubs)

    # Convert the matrix of gene indexes to a matrix of bridgeness values
    S1 <- matrix(lib1[iSubs], nrow=nrow(iSubs), ncol=ncol(iSubs))
    r1 <- rowSums(S1) / L
    # Scale bridgeness values to [0,1]
    r1 <- r1 / max(r1)

    # Convert the matrix of gene indexes to a matrix of betweenness values
    S2 <- matrix(lib2[iSubs], nrow=nrow(iSubs), ncol=ncol(iSubs))
    r2 <- rowSums(S2) / L
    # Scale betweenness values to [0,1]
    r2 <- r2 / max(r2)  

    # Convert the matrix of gene indexes to a matrix of betweenness values
    S3 <- matrix(lib3[iSubs], nrow=nrow(iSubs), ncol=ncol(iSubs))
    r3 <- rowSums(S3) / L

    res           <- matrix(,ncol=3, nrow=length(r1))
    res[,1]       <- r1
    res[,2]       <- r2
    res[,3]       <- r3
    colnames(res) <- c('subPathness', 'subBC', 'subDEG')
    rownames(res) <- 1:nrow(res)

    message('done.')

    return(res)
}

