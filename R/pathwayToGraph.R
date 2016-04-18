

createPathwayGraphs <- function(org, pathways, edgeTypes, doubleEdges, choice,
                                groupMode)
{
    if ( missing(edgeTypes) )
    {
        edgeTypes   <- getEdgeTypes()
    }
    if ( missing(doubleEdges) )
    {
        doubleEdges <- getDoubleEdges()
    }
    if ( missing(choice) )
    {
        choice <- 'reactions'
    }
    if ( missing(groupMode) )
    {
        groupMode <- 'expand'    
    }
    if (missing(org))
    {
        message('Please specify a valid organism id.')
        return(NULL)
    } 

    # Xml directory for organism
    xmlDir    <- paste(cache$dirs$xml, org, sep='//')

    # Choose valid pathways
    if (missing(pathways))  
    { 
        paths <- list.files(xmlDir) 
    }

    if (!missing(pathways)) 
    { 
        paths <- paste(org, pathways, '.xml', sep='') 
    }

    message('Creating adjacency matrices...', appendLF = FALSE)

    # Create compact adjacency matrices for given pathways.
    types  <- getPathwayType(paste(xmlDir, paths, sep='//'))
    N <- length(paths)
    cores  <- ifelse( N > detectCores()*10, 'default', 1 )
    export <- c(
        'nonMetabolicPathwayToGraph', 'expandMetabolicGraph', 
        'removeCompoundsNonMetabolicGraph', 
        'metabolicPathwayToGraph', 'expandNonMetabolicGraph',
        'removeCompoundsMetabolicGraph',
        'xmlTreeParse', 'xmlRoot', 'xmlGetAttr', 'xmlName', 'xmlChildren')

    cAdjMats <- .doSafeParallel(
                    funcName=pathwayToGraph,
                    export=export, 
                    combine='default',
                    N=length(paths),
                    cores=cores,
                    xmlDir, paths, types, FALSE, edgeTypes, 
                    doubleEdges, choice, groupMode)

    names(cAdjMats) <- gsub('.xml', '', paths)

    # Create expanded adjacency matrices with edge type information 
    eAdjMats <- .doSafeParallel(funcName=pathwayToGraph,
                                export=export, 
                                combine='default',
                                N=length(paths),
                                cores=cores,
                                xmlDir, paths, types, TRUE, edgeTypes, 
                                doubleEdges, choice, groupMode)

    names(eAdjMats) <- gsub('.xml', '', paths)

    message('done.')

    return(list(combined=cAdjMats, expanded=eAdjMats, org=org))
}


pathwayToGraph <- function (i, ...)
{

    # ------- Unpacking arguments -------
    args        <- list(...)
    dir         <- args[[1]]
    file        <- args[[2]][i]
    type        <- args[[3]][i]
    expand      <- args[[4]]
    edgeTypes   <- args[[5]]
    doubleEdges <- args[[6]]
    choice      <- args[[7]]
    groupMode   <- args[[8]]
    # ----------------------------------

    adjMat <- NULL
    if ( type == 'Metabolic'  &&  choice == 'reactions' )
    {
        path <- paste(dir, file, sep='//')
        gr   <- metabolicPathwayToGraph(path)
        if (length(gr) > 0)
        {
            if (!expand)
            {
                adjMat <- removeCompoundsMetabolicGraph(gr)
            }
            else
            {
                adjMat <- removeCompoundsMetabolicGraph(gr)
                adjMat <- expandMetabolicGraph(adjMat)
            }
        } 
    }
    if ( ( type == 'Metabolic'  &&  choice == 'relations' )  
        ||  type == 'Non-Metabolic' )
    {
        path <- paste(dir, file, sep='//')
        gr   <- nonMetabolicPathwayToGraph(path, doubleEdges, groupMode)

        if (!expand)
        {
            adjMat <- removeCompoundsNonMetabolicGraph(gr, unique=FALSE,
                                                    edgeTypes=edgeTypes)
        }
        else
        {
            adjMat <- removeCompoundsNonMetabolicGraph(gr, unique=FALSE,
                                                    edgeTypes=edgeTypes)
            adjMat <- expandNonMetabolicGraph(adjMat)
        }
    }

    return(adjMat)
}

#
# Find Graph Type
#
getPathwayType                    <- function(filepath, file)
{
    types <- vector(mode='numeric', length=length(filepath))
    for (i in 1:length(filepath))
    {
        num <- tail(unlist(strsplit(filepath[i], '//')), 1)
        num <- as.numeric(gsub("[^0-9]", "", num))
        if (num <= 2000)    { type = 'Metabolic'    }
        if (num > 2000)     { type = 'Non-Metabolic'}
        types[i] <- type
    }

    return(types)
}
#
# Graph from Metabolic Pathways
#

metabolicPathwayToGraph           <- function(filepath)
{
    xmlDoc <- tryCatch(xmlTreeParse(filepath,error=NULL),
            error=function(e) "error")
    if(!('XMLDocument' %in% class(xmlDoc) )) { return() }

    root          <- xmlRoot(xmlDoc)
    org           <- xmlGetAttr(root,"org")
    nodes         <- as.vector(mapply(xmlName, root))
    entryNodes    <- root[which(nodes == 'entry')]
    reactionNodes <- root[which(nodes == 'reaction')]

    # Create vertices from entry nodes
    idx           <- sapply(1:length(entryNodes), function(x) x)
    eid           <- sapply(entryNodes, xmlGetAttr, "id")        
    enames        <- gsub(paste(org,':',sep=''),'', 
                        sapply(entryNodes, xmlGetAttr, "name")) 
    etype         <- sapply(entryNodes, xmlGetAttr, "type")
    vertices      <- data.frame(ID=eid,id=eid,names=enames, type=etype, 
                    stringsAsFactors=FALSE)

    # Create edges from reaction nodes
    edges <- list(); e <- 1
    for (r in reactionNodes)
    {
        rids        <- sapply(xmlChildren(r), xmlGetAttr, "id")
        subsId      <- which(names(rids) == "substrate")
        prodId      <- which(names(rids) == "product")
        substrate   <- list()
        product     <- list()  
        for (s in 1:length(subsId)) 
        { substrate[[s]] <- list(id=rids[subsId[s]]) }
        for (p in 1:length(prodId)) 
        { product[[p]]   <- list(id=rids[prodId[p]]) }

        enzyme      <- xmlGetAttr(r,"id")
        substrateId <- unname(sapply(substrate, '[[', "id"))
        productId   <- unname(sapply(product,   '[[', "id"))

        #  Non Reversible 
        for (i in 1:length(substrateId))     
        { edges[[e]] <- c(substrateId[i], enzyme); e <- e + 1 }
        for (i in 1:length(productId))       
        { edges[[e]] <- c(enzyme, productId[i]);   e <- e + 1 }
        
        #  Double the edges if reversible  
        if(xmlGetAttr(r,"type") == "reversible")
        {             
            for (i in 1:length(productId))   
            { edges[[e]] <- c(productId[i], enzyme);   e <- e + 1 }
            for (i in 1:length(substrateId)) 
            { edges[[e]] <- c(enzyme, substrateId[i]); e <- e + 1 }
        }
    }  
    if (length(edges) > 0 )
    {
        edges <- unique(data.frame(do.call(rbind, edges), 
                                    stringsAsFactors=FALSE))
        names(edges) <- c('e1','e2')  
    }
    res <- list(name=xmlGetAttr(root,"name"), edges=edges,vertices=vertices)


    return(res)
}

removeCompoundsMetabolicGraph     <- function(path)
{
    nodeType <- 'gene'   
    if(path$name != gsub('ec','',path$name)) { nodeType<-"enzyme" }
    enzymes  <- which(path$vertices$type == nodeType)
    vid      <- path$vertices$id
    entry1 <- c(); entry2 <- c()

    if ( length(path$edges) > 0 )
    {
        # For each enzyme node find the neighbors. 
        # The neighbors that are of type 'compound' are to be removed, and 
        # newconnections have to be established with the appropriate enzymes.
        for(j in 1:length(enzymes) )   
        {       
            L <- c(); R <- c(); 
            ec  <- vid[enzymes[j]]
            # For each neighbor of the enzyme
            for (r1 in path$edges[path$edges$e1 == 
                                path$vertices[,'id'][enzymes[j]],]$e2)  
            {
                # if neighbor is a compound, find its neighbors
                for (r2 in path$edges[path$edges$e1 == 
                                path$vertices[,'id'][which(vid == r1)],]$e2)
                {
                    # for enzyme neighbors of the compound others than 
                    # the initial enzyme
                    nid <- vid[which(path$vertices$id == r2)]
                    if (!(nid %in% R) && nid != ec) 
                    { 
                        L <- c(L,ec)
                        R <- c(R,nid) 
                    }
                }
            }
            entry1 <- c(entry1,L)
            entry2 <- c(entry2,R)
        }
    } 


    # Construct the adjacency matrix
    xid    <- path$vertices$id[enzymes]
    names  <- path$vertices$names[enzymes]       
    adjMat <- matrix(rep(0,length(enzymes)^2), nrow=length(enzymes))
    colnames(adjMat) <- names
    rownames(adjMat) <- names
    for(j in 1:length(entry1) )  
    {
        adjMat[which(xid == entry1[j]), which(xid == entry2[j])] <- 1
    }
    
    return(adjMat) 
}

expandMetabolicGraph              <- function(adjMat)
{
    names  <- c()
    Vcount <- nrow(adjMat)
    Ecount <- sum(adjMat)
    if (Vcount == 0) { return(NULL); }

    n1 <- c(); n2 <- c()
    for ( j in 1:nrow(adjMat) )
    {
        for ( z in 1:ncol(adjMat) )
        {
            if (adjMat[j,z] == 1)
            {
                n1  <- c(n1, colnames(adjMat)[j] )
                n2  <- c(n2, colnames(adjMat)[z] )
            }     
        }
    }

    # Expand gene names
    for(j in 1:Vcount)
    {
        for(vk in unlist(strsplit(colnames(adjMat)[j]," ")) )
        {   
            names <- c(names, vk)
        }
    }
    names <- unique(names)

    L <- c(); R <- c()
    if (Ecount > 0)
    {
        for(j in 1:Ecount)
        {
            for(lm in unlist(strsplit(n1[j]," ")))
            {
                for(rn in unlist(strsplit(n2[j]," ")))
                {
                    L <- c(L, lm)
                    R <- c(R, rn)
                }
            }
        }
    }

    # Construct the adjacency matrix
    N      <- length(names)
    adjMat <- matrix(rep(0,N^2), nrow=N)
    for(j in 1:length(L) )  
    {
        t <- 3  # Unknown  
        adjMat[which(names == L[j]), which(names == R[j])] <- t
    }
    colnames(adjMat) <- names 
    rownames(adjMat) <- names

    # Remove self cycles
    for(j in 1:nrow(adjMat) )  { adjMat[j,j] <- 0 }

    return(adjMat)
}

#
# Graph from Mon Metabolic Pathways
#

nonMetabolicPathwayToGraph <- function(filepath, doubleEdges, groupMode)
{

    xmlDoc         <- tryCatch(xmlTreeParse(filepath,error=NULL),
                                error=function(e) "error")
    if(!('XMLDocument' %in% class(xmlDoc) )) { return(NULL) }

    # Create node families
    root           <- xmlRoot(xmlDoc)
    org            <- xmlGetAttr(root,"org")
    nodes          <- as.vector(mapply(xmlName, root))
    entryNodes     <- root[which(nodes == 'entry')]
    type           <- sapply(entryNodes, xmlGetAttr, "type")
    groupNodes     <- entryNodes[which(type == 'group')]
    orthoNodes     <- entryNodes[which(type != 'group')]
    relationNodes  <- root[which(nodes == 'relation')]
    
    # Create vertices from entry nodes
    vertices       <- list()
    vertices$id    <- as.integer(unname(sapply(orthoNodes, xmlGetAttr, "id")))
    vertices$names <- gsub(paste(org,':',sep=''), '', 
                        unname(sapply(orthoNodes, xmlGetAttr, "name")))
    vertices$type  <- unname(sapply(orthoNodes, xmlGetAttr, "type"))
    group          <- list()

    if (groupMode == 'expand')
    {
        # Create groups
        for(gn in groupNodes) 
        { 
            group <- c( group, 
                            list(sapply(xmlChildren(gn), xmlGetAttr, "id")) )
        }
        names(group) <- sapply(groupNodes, xmlGetAttr, "id")
    }

    if (groupMode == 'collapse')
    {   
        # Replace group nodes with all genes of their components
        for(gn in groupNodes) 
        { 
            entries <- unlist(lapply(xmlChildren(gn), xmlGetAttr, "id"))
            genes   <- vertices$names[which(vertices$id %in% entries)]
            genes   <- unique(unlist(lapply(as.list(genes), 
                        function(x) { strsplit(x, ' ') } )))
            genes   <- paste(genes, collapse=' ')

            vertices$id     <- c(vertices$id,    xmlGetAttr(gn, "id"))
            vertices$names  <- c(vertices$names, genes)
            vertices$type   <- c(vertices$type,  "gene")
        }
    }

    # Create relations
    rentry1  <- sapply(relationNodes, xmlGetAttr, "entry1")
    rentry2  <- sapply(relationNodes, xmlGetAttr, "entry2")
    rtype    <- sapply(relationNodes, xmlGetAttr, "type")
    rents    <- lapply(relationNodes, function(x) { xmlChildren(x) })
    rsubtype <- list()

    for(ent in rents)
    {
        subtype <- mapply(list, name='unknow', value='unknow', SIMPLIFY=FALSE)
        if(length(ent) > 0)
        {
            X       <- unname(sapply(ent, xmlGetAttr, "name"))
            Y       <- unname(sapply(ent, xmlGetAttr, "value"))
            subtype <- mapply(list, name=X, value=Y, SIMPLIFY=FALSE)
        }       
        rsubtype <- append(rsubtype, list(subtype))
    }

    bla <- function(e1, e2, name='unknow', type='noType', duplicate = FALSE)
    {
        L1 <- list(e1=e1, e2=e2, name=name, type=type)
        L2 <- list(e1=e2, e2=e1, name=name, type=type)

        if (!duplicate) { return( list( L1 ))  }
        if (duplicate)  { return( list( L1, L2 )) }
    }

    # Create edges
    E <- list(); i <- 0      
    for(subtypes in rsubtype)
    {
        i       <- i + 1
        inGroup <- FALSE
        # type    <- paste(names(subtypes), collapse='_')        
        if(length(group) != 0) 
        {
            # If a node is of type group, then create all actual edges between
            # its components and the other node. 
            GroupL <- groupL <- unlist(group[names(group) %in%  rentry1[i]])
            GroupR <- groupR <- unlist(group[names(group) %in%  rentry2[i]]) 
            if(length(groupL) == 0 && length(groupR)  > 0) 
            { GroupL <- rentry1[i] }
            if(length(groupL) >0   && length(groupR) == 0) 
            { GroupR <- rentry2[i] }

            # Create an single edge by default; if the edge is of type 
            # compound, create the actual edges connecting the compounds.
            # For an ambiguous edge type, create an edge in the other 
            # direction. 
            for(rn in subtypes)
            {
                dupl <- rn$name %in% doubleEdges
                
                for(gl in GroupL)
                {
                    for(gr in GroupR)
                    {                    
                        if (rn$name == "compound")                     
                        {    
                            E <- c(E, bla(gl,       rn$value, rn$name, 
                                    dup = dupl, type=rn$name))
                            E <- c(E, bla(rn$value, gr,       rn$name, 
                                    dup = dupl, type=rn$name))
                        }
                        if (rn$name != "compound")
                        {
                            E <- c(E, bla(gl, gr, rn$name, 
                                    dup = dupl, type=rn$name))
                        }
                    }
                }

                if (length(groupL) >0 || length(groupR) > 0) 
                { 
                    inGroup = TRUE 
                }
            }
        }
        if (inGroup) { next() }
        

        # No group nodes
        for(rn in subtypes)
        {
            dupl <- rn$name %in% doubleEdges            
            if (rn$name == "compound")
            {
                E <- c(E, bla(rentry1[i],rn$value,   dup=dupl, type=rn$name))
                E <- c(E, bla(rn$value,  rentry2[i], dup=dupl, type=rn$name))
            }                        
            if (rn$name != "compound")
            {
                E <- c(E, bla(rentry1[i], rentry2[i], rn$name, 
                        dup = dupl, type=rn$name))  
            }
        }
    }
    
    e1    <- sapply(E, function(x){ x$e1 })
    e2    <- sapply(E, function(x){ x$e2   })
    ename <- sapply(E, function(x){ x$name })
    etype <- sapply(E, function(x){ x$type })
    edges <- unique(data.frame(e1=e1, e2=e2, ename=ename, etype=etype))
    if (nrow(edges) > 0)
    {
        edges <- list(e1=edges[,1], e2=edges[,2], type=edges[,3])
    }    
    res <- list(vertices=vertices, edges=edges, name=xmlGetAttr(root,"name"))

    return(res)
}

removeCompoundsNonMetabolicGraph <- function(path, unique, edgeTypes)
{
    if (is.null(path)) return(NULL)
    vid      <- as.numeric(path$vertices$id)
    etype    <- path$vertices$type
    nodeType <- 'gene'   
    if(path$name != gsub('ko','',path$name)) { nodeType <- "ortholog" }

    source    <- c(); destin <- c(); types <- c()          
    genesIndx <- which(path$vertices$type == nodeType)
    for(gi in genesIndx)
    {               
        # For each gene, find its neighboring vertices. 
        # KEGG entries are not always numbered successively. 
        L <- c(); R <- c(); TT <- c()
        neighbors <- path$edges$e2[path$edges$e1 == vid[gi]]

        for(neighbor in neighbors)
        { 
            nbrId <- which( vid == neighbor )

            # The neighbor is a gene which has yet to be added to 
            # the edgelist.
            if( etype[nbrId] == nodeType &&  !(vid[nbrId] %in% R) )
            {
                # Connect the first gene with the neighbor        
                L    <- c( L, vid[gi]    )
                R    <- c( R, vid[nbrId] )

                # Edges with more than one types have been stored seperately.
                # Join all their types in a single entry the first time the 
                # edge has been found and make sure it is not added multiple 
                # times (2nd if clause).
                idx1 <- which( path$edges$e1 == vid[gi] )
                idx2 <- which( path$edges$e2 == vid[nbrId] )
                idx  <- intersect(idx1, idx2)
                TT   <- c( TT, paste((path$edges$type[idx]), collapse='_') )
            }

            # The neighbor is a compound. 
            if( etype[nbrId] == "compound" )
            {                
                # The compound has to be removed, while the source gene is 
                # connected with the compound's neighboring gene.
                cpdNeighbors <- path$edges$e2[ 
                                        which(path$edges$e1 == vid[nbrId]) ]

                for( cpdNeighbor in cpdNeighbors )
                {
                    idx  <- which( vid == cpdNeighbor )
                    # The compound's neighbor is a gene which has yet to be 
                    # added to the edgelist.
                    if( etype[idx] == nodeType && !(vid[idx] %in% R) )
                    {
                        # Connect source gene with the neighbor 
                        L <- c( L, vid[gi]   )
                        R <- c( R, vid[idx]  )
                        # The edgetype is 'compound'
                        TT <- c( TT, "compound")
                    }
                }
            }
        }     
        source <- c( source, L )
        destin <- c( destin, R )
        types  <- c( types, TT  )
    }


    # Construct the adjacency matrix
    if (unique)
    {
        # Recalculate indexes according to unique gene names
        names            <- unique(path$vertices$names[genesIndx])
        adjMat           <- matrix(rep(0,length(names)^2), nrow=length(names))
        colnames(adjMat) <- names
        rownames(adjMat) <- names

        if (is.null(source) || is.null(destin))
        {
            return(adjMat)
        }

        for (i in 1:length(source))
        {   
            idx1 <- which(path$vertices$id == source[i])
            idx2 <- which(path$vertices$id == destin[i])
            source[i] <- names[ names == path$vertices$names[idx1] ]
            destin[i] <- names[ names == path$vertices$names[idx2] ]
        }

        for(j in 1:length(source) )  
        {
            idx <- which(edgeTypes[,1] == types[j])
            if (length(idx) == 0)
            {
                # Set new interaction types to apathetic
                t <- 3
            }
            if (length(idx) > 0)
            {
                t <- edgeTypes[idx, 2]
            }
            
            adjMat[which(names == source[j]), which(names == destin[j])] <- t
        }
    }
    else
    {
        # Construct the adjacency matrix
        gids                <- path$vertices$id[genesIndx]
        names               <- unname(path$vertices$names[genesIndx])
        adjMat              <- matrix(rep(0,length(genesIndx)^2), 
                                        nrow=length(genesIndx))
        colnames(adjMat)    <- names 
        rownames(adjMat)    <- names

        if (is.null(source) || is.null(destin))
        {
            return(adjMat)
        }

        for(j in 1:length(source) )  
        {
            
            idx <- which(edgeTypes[,1] == types[j])
            if (length(idx) == 0)
            {
                t <- 3; # New interaction type
            }
            if (length(idx) > 0)
            {
                t <- edgeTypes[idx, 2]
            }

            adjMat[which(gids == source[j]), which(gids == destin[j])] <- t
        }
    }

    # Remove self loops
    for (j in 1:nrow(adjMat)) { adjMat[j,j] <- 0 }


    return(adjMat) 
}

expandNonMetabolicGraph          <- function(adjMat)
{
    if (!is.matrix(adjMat)) return(NULL)

    names  <- c()
    Vcount <- nrow(adjMat)
    Ecount <- sum(adjMat)
    if (Vcount == 0) { return(NULL); }

    n1 <- c(); n2 <- c(); t <- c()
    for ( j in 1:nrow(adjMat) )
    {
        for ( z in 1:ncol(adjMat) )
        {
            if (adjMat[j,z] > 0)
            {
                n1   <- c(n1, colnames(adjMat)[j] )
                n2   <- c(n2, colnames(adjMat)[z] )
                t    <- c(t, adjMat[j,z])
            }     
        }
    }

    # Expand gene names
    names <- unique(unlist(lapply(colnames(adjMat), 
                function(x) {strsplit(x,' ')} )))

    L <- c(); R <- c(); TT <- c()
    if (Ecount > 0)
    {
        for(j in 1:Ecount)
        {
            l1 <- unlist(strsplit(n1[j],' '))
            l2 <- unlist(strsplit(n2[j],' '))
            for(lm in l1) 
            { 
                for(rn in l2) 
                { 
                    L   <- c(L, lm)
                    R   <- c(R, rn) 
                    TT  <- c(TT, t[j])
                } 
            }
        }
    }
    L   <- if(!is.null(L)) { L <- na.omit((L)) } 
    R   <- if(!is.null(L)) { R <- na.omit((R)) } 
    TT  <- if(!is.null(TT)) { TT <- na.omit((TT)) } 


    # Construct the adjacency matrix
    N      <- length(names)
    adjMat <- matrix(rep(0,N^2), nrow=N)
    for(j in 1:length(L) )  
    {
        adjMat[which(names == L[j]), which(names == R[j])] <- TT[j]
    }
    colnames(adjMat) <- names
    rownames(adjMat) <- names

    # Remove self cycles
    for(j in 1:nrow(adjMat) )  { adjMat[j,j] <- 0 }

    return(adjMat)
}



getEdgeTypes   <- function(type)             
{
    # noEdge     0
    # activation 1
    # inhibition 2
    # apathetic  3

    data = c(
        'unknow', 3,
        'activation', 1,
        'inhibition', 2,
        'binding/association', 3,
        'expression', 1,
        'repression', 2,
        'phosphorylation', 3, #1
        'dephosphorylation', 3,
        'ubiquitination', 3,
        'dissociation', 3,
        'indirect effect', 3, #1
        'state change', 3,
        'compound', 3,  #1
        'hidden compound', 3,
        'missing interaction', 3,
        'activation_phosphorylation', 1,
        'activation_dephosphorylation', 1,
        'activation_ubiquitination', 1,
        'activation_indirect effect', 1,
        'activation_binding/association', 1,
        'activation_inhibition', 3,
        'activation_methylation', 1,
        'activation_phosphorylation_binding/association', 1,
        'activation_phosphorylation_indirect effect', 1,
        'inhibition_phosphorylation', 2,
        'inhibition_dephosphorylation', 2,
        'inhibition_ubiquitination', 2,
        'inhibition_indirect effect', 2,
        'inhibition_binding/association', 2,
        'inhibition_expression', 2,
        'inhibition_methylation', 2,
        'compound_expression', 1,
        'compound_activation', 1,
        'compound_inhibition', 2,
        'compound_activation_indirect effect', 1,
        'compound_activation_phosphorylation', 1,
        'phosphorylation_indirect effect', 3,
        'phosphorylation_binding/association', 3,
        'phosphorylation_dissociation', 3,
        'dephosphorylation_indirect effect', 3,
        'binding/association_missing interaction', 3,
        'binding/association_indirect effect', 3,
        'expression_indirect effect', 1,
        'repression_indirect effect', 2,
        'ubiquitination_inhibition', 2,
        'dissociation_missing interaction', 3,
        'indirect effect_phosphorylation', 3,
        'inhibition_dephosphorylation_indirect effect', 3,  #new
        'inhibition_glycosylation',3,                       #new
        'glycosylation',3                                   #new
        )


    E        <- as.numeric(data[seq(2,length(data),2)])
    names(E) <- data[seq(1,length(data),2)]


    if (!missing(type))
    {
        return(E[type])
    }
    if (missing(type))
    {
        df <- data.frame(names(E),E)
        rownames(df) <- 1:nrow(df)
        colnames(df) <- c('Edge', 'Type')
        
        return(df)
    }

    return( E[type] )
}


getDoubleEdges <- function()
{
    doubleEdges <- c('compound', 
                    'binding/association', 
                    'dissociation', 
                    'dissociation_missing interaction', 
                    'state change', 
                    'binding/association_missing interaction', 
                    'missing interaction', 'hidden compound')

    return(doubleEdges)
}

