

visualizeResults <- function(summary, export, expand, colors, from, to)
{
    message('Summarising scores...', appendLF = FALSE)
    # 
    # Merging data to presentable form.
    # 
    subpaths      <- summary$subpaths
    pValues       <- summary$pValues
    subScore      <- summary$subScore
    sb            <- summary$measures[,1] # subPathness
    sp            <- summary$measures[,2] # subBC
    sd            <- summary$measures[,3] # subDEG 
    org           <- summary$org
    filters       <- summary$filters
    type          <- summary$type
    miRNAoverSubs <- as.vector(summary$miRNAsOverSubpathway)


    if (is.null(subpaths) || nrow(subpaths) == 0) 
    { 
        message('done.'); return(invisible()) 
    }
    if (missing(export)) 
    { 
        export <- c('.txt', '.xlsx') 
    }
    if (missing(expand)) 
    { 
        expand <- TRUE 
    }
    if (missing(colors)) 
    { 
        colors <- c('#FFFFFF', '#DCFFDC', '#99FF99', '#00EB00', '#00E100') 
    }
    if (missing(from)) { from <- 'entrezgene' } 
    if (missing(to)) { to <- 'entrezgene' } 

    # Round data
    pValues     <- format(ceiling(pValues*10^2)/10^2, scientific=FALSE)
    subScore    <- ceiling(subScore*10^3)/10^3
    sb          <- ceiling(sb*10^2)/10^1
    sp          <- ceiling(sp*10^2)/10^1
    sd          <- ceiling(sd*10^2)/10^1 

    # Creating directories    
    .createDirectories(org=org)


    # Map subpath ids to pathway names 
    pathways    <- downloadKEGGPathwayList(org)
    subpathIds  <- gsub(org, '', rownames(subpaths))
    y           <- pathways[,2]
    names(y)    <- pathways[,1]
    subNames    <- y[subpathIds]


    if (type == 'Linear')     
    {
        d   <- paste(cache$dirs$slnr, org, sep='//')
        f1  <- createScoreFile(dir=d, type='a', fileType='.txt')
        f2  <- createScoreFile(dir=d, type='a', fileType='.xlsx')
    }
    if (type == 'Non-Linear') 
    { 
        d   <- paste(cache$dirs$snlr, org, sep='//')
        f1  <- createScoreFile(dir=d, type='a', fileType='.txt')
        f2  <- createScoreFile(dir=d, type='a', fileType='.xlsx') 
    }

    pathnames <- list()

    # Subpath members (by removing zero placeholders)
    if (type=='Linear')
    {
        # Convert nomenclature
        ids <- unique(as.vector(subpaths))

        psubs2 <- NULL
        if ( from != to ) 
        {
            y <- convertNomenclature(ids=ids, org=org, from=from, to=to)

            if (nrow(y) > 1)
            {
                lib         <- y[,2]
                names(lib)  <- y[,1]
                
                z <- matrix(,nrow=nrow(subpaths), ncol=ncol(subpaths))
                for (i in 1:nrow(subpaths))
                {
                    z[i,] <- lib[as.character(subpaths[i,])]
                }    
                psubs2 <- apply(z, 1, paste, collapse=' ')
                psubs2 <- gsub(' NA',' ',psubs2)  # Remove zeros
            }else 
            {
                psubs2 <- apply(subpaths, 1, paste, collapse=' ')
                psubs2 <- gsub(' 0',' ',psubs2)  # Remove zeros          
            }
        }

        psubs1 <- apply(subpaths, 1, paste, collapse=' ')
        psubs1 <- gsub(' 0',' ',psubs1)  # Remove zeros
    }

    if (type=='Non-Linear')
    {
        usubs <- vector(mode='list', length=nrow(subpaths))
        for(i in 1:nrow(subpaths))
        {
            usub       <- unique(unlist(strsplit(subpaths[i,], split='-')))
            usubs[[i]] <- matrix(as.numeric(usub[which(usub!='0')]), nrow=1)

        }
        subs <- unlistToMatrix(fillMatrixList(usubs))
        ids  <- unique(as.vector(subs))

        psubs2 <- NULL
        if ( from != to ) 
        {
            y <- convertNomenclature(ids=ids, org=org, from=from, to=to)
            if (nrow(y) > 1)
            {
                lib         <- y[,2]
                names(lib)  <- y[,1]
                
                z <- matrix(,nrow=nrow(subs), ncol=ncol(subs))
                for (i in 1:nrow(subs))
                {
                    z[i,] <- lib[as.character(subs[i,])]
                }    
                psubs2 <- apply(z, 1, paste, collapse=' ')
                psubs2 <- gsub(' NA',' ',psubs2)  # Remove zeros
            }else 
            {
                psubs2 <- apply(subs, 1, paste, collapse=' ')
                psubs2 <- gsub(' 0',' ',psubs2)  # Remove zeros          
            }
        }

        psubs1 <- apply(subs, 1, paste, collapse=' ')
        psubs1 <- gsub(' 0',' ',psubs1)  # Remove zeros
    }

    # Create the links
    links        <- createKEGGLinks(subpaths, type=type, openInBrowser=FALSE)
    names(links) <- paste('S', 1:nrow(subpaths), sep='')

    # Create the main data frame
    titles <- paste('(', paste(org, names(subNames), ') ', sep=''), subNames, 
                    sep='')
    df     <- data.frame('Subpath'=links, 'Pathway'=titles, 
                    stringsAsFactors=FALSE)
    class(df$Subpath) <- c('hyperlink')

    # isMetabolic
    rnames        <- as.numeric(gsub(org, '', rownames(subpaths)))
    names(rnames) <- rownames(subpaths)
    isMetabolic   <- ifelse(rnames < 2000, '+', '-')
    df            <- cbind(df, data.frame('isMetabolic'=isMetabolic, 
                                            stringsAsFactors=FALSE)) 

    if ('subScore' %in% names(filters))  
    {
        # Append the gene scores for all time points
        dft           <- data.frame(subScore, stringsAsFactors=FALSE)
        colnames(dft) <- paste('T', 1:ncol(subScore), sep='')
        df            <- cbind(df, dft)
    }

    if ('p-value' %in% names(filters))           
    {
        df <- cbind(df, data.frame('q-value'=pValues, stringsAsFactors=FALSE))
    }


    if ('measures' %in% names(filters))           
    {
        df <- cbind(df, data.frame('subPathness'=sb, stringsAsFactors=FALSE))
        df <- cbind(df, data.frame('subBC'=sp, stringsAsFactors=FALSE)) 
        df <- cbind(df, data.frame('subDEG'=sd, stringsAsFactors=FALSE))
    }

    from.name      <- .getFilters()[from]
    to.name        <- .getFilters()[to]
    dft1           <- data.frame(psubs1, stringsAsFactors=FALSE)
    colnames(dft1) <- from.name
    if ( !is.null(psubs2) )
    {
        dft2           <- data.frame(psubs2, stringsAsFactors=FALSE)
        colnames(dft2) <- to.name
        df             <- cbind(df, dft2, dft1)
        rr <- 2
    }
    if ( is.null(psubs2) )
    {
        df <- cbind(df, dft1)
        rr <- 1
    }


    # Create a header mapper
    header <- c('baseStart'=1, 'baseStop'=3)
    header <- c(header, 'scoreStart'=(3+1), 'scoreStop'=(3+ncol(subScore)))
    right  <- (3 + ncol(subScore))
    if ('p-value' %in% names(filters))
    {
        header     <- c(header, 'p-value'=(right + 1))
        right     <- (right + 1)
    }
    if ('measures' %in% names(filters))
    {
        header <- c(header, 'measuresStart'=(right + 1), 
                    'measuresStop'=(right + 3))
        right  <- (right + 3)
    }    
    header     <- c(header, 'memberStart'=(right + 1), 
                    'memberStop'=(right + rr))
    right     <- (right + rr)
    if ('mirScore' %in% names(filters))
    {
        header <- c(header, 'mirStart'=(right + 1), 'mirStop'=(right + 1))
        right  <- (right + 1)
    }

    if ('mirScore' %in% names(filters))        
    {
        df <- cbind(df, data.frame('miRNA'=miRNAoverSubs, 
                                stringsAsFactors=FALSE))
    }
    if ('.xlsx' %in% export  &&  expand)
    { 
        .createExcelOutputMulti(df=df, links=links, file=f2, filters=filters,
                                header=header, colors=colors)
    }
    if ('.xlsx' %in% export  &&  !expand)
    { 
        .createExcelOutputSingle(df=df, links=links, file=f2, filters=filters,
                                header=header, colors=colors)
    }
    if ('.txt' %in% export)
    {
        df[,1] <- names(links)
        exportToFile(df, f1, append=FALSE)        
    }


    message('done.')

    return(df)
}

.createExcelOutputSingle <- function(df, links, file, filters, header, colors)
{
    # Specify columns according to specified filters
    ncols   <- 2 # Members
    hypCols <- 0; mesCols <- 0; mirCols <- 0
    if ('p-value' %in% names(filters))       
    { 
        ncols <- ncols + 1
        hypCols <- 1
    }
    if ('measures' %in% names(filters))      
    { 
        ncols <- ncols + 3 
        mesCols <- 3
    }
    if ('mirScore' %in% names(filters))      
    { 
        ncols <- ncols + 1
        mirCols <- 1
    }

    # Determime columns for each group
    subScoreRange <- 4 : (ncol(df) - ncols)

    if ('measures' %in% names(filters))      
    {
        measuresRange <- (ncol(df) - ncols + 1 + hypCols) : 
                            (ncol(df) - ncols + hypCols + mesCols) 
    }

    a             <- (ncol(df) - ncols + hypCols + mesCols)
    membersRange  <- (a + 1) : (a + 2)


    # Bind a row on top of df with column groupings 
    # (Subscores, Measures, miRNA Mediated Subpathway Members)
    cnames <- colnames(df)
    cnames <- paste(' ', cnames, ' ', sep='')
    df     <- rbind(cnames,df)
    cnames <- c('', '', '')
    cnames <- c(cnames, rep('Subscores', length(subScoreRange)), '')
    
    if ('measures' %in% names(filters))       
    {
        cnames <- c(cnames, rep('Measures',3))    
    }

    cnames <- c(cnames, rep('miRNA Mediated Subpathway Members',2), '')
    df     <- rbind(cnames, df)

    wb    <- createWorkbook()
    addWorksheet(wb, 'Results')
    cl    <- lapply(df, function(x) tolower(class(x)))
    wb$writeData(sheet='Results', df=df, startRow=1, startCol=1, 
                colClasses = cl, hlinkNames = c('','SubId',names(links)), 
                colNames=FALSE, keepNA=FALSE) 


    # Customize output
    setColWidths(wb, sheet='Results', cols=1:ncol(df), widths='auto')
    modifyBaseFont(wb, fontSize=10, fontColour='#000000', fontName="Segoe UI")
    freezePane(wb, sheet='Results', firstActiveRow=3, firstActiveCol=2)


    # Main Body
    style <- createStyle(fontSize=10, fontName='Segoe UI', halign="left", 
                        valign = "center", fgFill='#ffffff', 
                        borderColour='#000000')
    addStyle(wb, sheet='Results', style=style, rows=2:(nrow(df)+1),
                        cols=1:ncol(df), gridExpand=TRUE)

    # Gene scores and topologics
    style <- createStyle(fontSize=10, fontName='Segoe UI', halign="center", 
                        valign = "center", fgFill='#ffffff', 
                        borderColour='#000000')
    addStyle(wb, sheet='Results', style=style, rows=2:(nrow(df)+1), 
                cols=3:ncol(df), gridExpand=TRUE)

    # Gene members and mirnas
    style <- createStyle(fontSize=10, fontName='Segoe UI', halign="left", 
                        valign = "center", fgFill='#ffffff', 
                        borderColour='#000000')

    if (!'mirScore' %in% names(filters)) 
    { 
        tcols <- ncol(df) 
    } 
    if ('mirScore' %in% names(filters))   
    { 
        tcols <- (ncol(df)-1):ncol(df) 
    }
    addStyle(wb, sheet='Results', style=style, rows=2:(nrow(df)+1), 
        cols=tcols, gridExpand=TRUE)


    # Header
    mergeCells(wb, sheet='Results', rows=1, cols=subScoreRange)
    if ('measures' %in% names(filters))          
    {
        mergeCells(wb, sheet='Results', rows=1, cols=measuresRange)
    }
    mergeCells(wb, sheet='Results', rows=1, cols=membersRange)


    style <- createStyle(fontSize=11, fontName='Segoe UI', 
            fontColour='#000000', halign="center", valign = "center", 
            fgFill='#ffffff', borderColour = '#000000')
    addStyle(wb, sheet='Results', style=style, rows=1, cols=1:ncol(df),
            gridExpand=TRUE)


    # Header Column
    style <- createStyle(fontSize=11, fontName='Segoe UI', halign="left", 
                        valign = "center", fgFill='#ffffff', 
                        borderColour='#000000')
    addStyle(wb, sheet='Results', style=style, rows=1:(nrow(df)+1), cols=1:1, 
                        gridExpand=TRUE)

    setRowHeights(wb, sheet='Results', rows = 1, heights = 35)
    setRowHeights(wb, sheet='Results', rows = 2:(nrow(df)+1), heights = 25)

    # Split gene score matrix to part based on the score
    S      <- df[, 4:(ncol(df) - ncols), drop=FALSE]
    rn     <- 1:12
    breaks <- unname(sapply(split(rn,cut(rn,quantile(rn,(0:k)/k), 
                    include.lowest=TRUE)), function(x) { x[2] }))
    breaks <- c(0, breaks/10)
    k      <- length(colors) + 1

    for (i in 2:k)
    {
        idx  <- which((S >= breaks[i-1]) & (S < breaks[i]))
        rows <- ((idx-1) %% nrow(S)) + 1
        cols <- floor((idx-1)/nrow(S)) + 1
        # Set cell colors of gene scores according to their scores
        style <- createStyle(fontSize=10, fontName='Segoe UI', 
            halign="center", valign = "center", fgFill=colors[i-1])
        addStyle(wb, sheet='Results', style=style, rows=0+rows, cols=3+cols,
            gridExpand=FALSE)
    }

    saveWorkbook(wb, file, overwrite=TRUE)  
}

.createExcelOutputMulti  <- function(df, links, file, filters, header, colors)
{

    # Find number of necessary miRNA columns 
    mirCols <- 0
    if ('mirScore' %in% names(filters))      
    {
        mirCols <- 1
        for (i in 1:nrow(df))
        {
            mirs = unlist(strsplit(df[i,header['mirStart']], ' '))
            if (length(mirs) > mirCols) { mirCols <- length(mirs) }
        }
    }

    # Find number of necessary members columns 
    genCols <- 0
    for (i in 1:nrow(df))
    {
        gens = unlist(strsplit(df[i,header['memberStart']], ' '))
        if (length(gens) > genCols) { genCols <- length(gens) }
    }
    genCols <- (header['memberStop'] - header['memberStart'] + 1)*genCols

    x           <- matrix('', nrow=nrow(df), ncol=mirCols + genCols)
    colnames(x) <- rep(' ', mirCols + genCols)

    # Fill first b columns of expansion with individual genes
    for (i in 1:nrow(x))
    {
        if (header['memberStop'] - header['memberStart'] > 0)
        {
            # Symbol
            geneSymbols <- unlist(strsplit(df[i,header['memberStart']], ' '))
            x[i, 1:length(geneSymbols)]  <- geneSymbols

            # Entrez
            geneEntrez <- unlist(strsplit(df[i,header['memberStop']], ' '))
            idx1 <- length(geneSymbols) + 1
            idx2 <- length(geneSymbols) + length(geneEntrez)
            x[i, idx1: idx2]  <- as.integer(geneEntrez)
        }
        if (header['memberStop'] - header['memberStart'] == 0)
        {
            geneEntrez <- unlist(strsplit(df[i,header['memberStart']], ' '))
            x[i, 1:length(geneEntrez)]  <- geneEntrez
        }

        if ('mirScore' %in% names(filters))      
        {    
            mirs  <- unlist(strsplit(df[i,header['mirStart']], ' '))
            if (length(mirs) > 0)
            {
                x[i, (genCols + 1) : (genCols + length(mirs))] <- mirs    
            }
        }
    }
    # Remove empty columns
    idx <- which(is.na(x))
    if ( length(idx) > 0 ) { x[idx] <- ''  }
    idx <- table(which(x == '', arr.ind=TRUE)[, 'col'])
    idx <- as.numeric(names(idx)[which(unname(idx) == nrow(x))])

    if ( length(idx) > 0 ) 
    { 
        x <- x[, -idx, drop=FALSE]
        genCols <- genCols - length(idx)
    }

    from     <- colnames(df)[header['memberStart']]
    to       <- colnames(df)[header['memberStop']]
    qvalCols <- 0
    mesCols  <- 0
    if ('p-value'  %in% names(filters))   { qvalCols<- 1 }
    if ('measures'  %in% names(filters))  { mesCols <- 3 }
    if ('mirScore' %in% names(filters))
    {
        idx <- c(header['memberStart'], header['memberStop'], 
                header['mirStart'])
    }
    if (!'mirScore' %in% names(filters))  
    {
        idx <- c(header['memberStart'], header['memberStop'])
    }
    df <- cbind(df[,-idx], as.data.frame(x, stringsAsFactors=FALSE))


    # Specify columns according to specified filters
    subScoreRange <- 4 : (ncol(df) - (mirCols + genCols + qvalCols + mesCols))
    a             <- tail(subScoreRange,1) + 1

    if ('p-value' %in% names(filters))   
    {
        qvaluesRange <- a : a
        a            <- tail(qvaluesRange,1) + 1  
    }
    if ('measures' %in% names(filters))      
    {
        measuresRange <- a : (a + 2)
        a             <- tail(measuresRange,1) + 1
    }

    membersRange <- a : (a + genCols - 1) 

    if ('mirScore' %in% names(filters))      
    { 
        a        <- tail(membersRange,1) + 1
        mirRange <- a : (a + mirCols - 1)
    }

    # Bind a row on top of df with column groupings 
    # (Subscores, Measures, miRNA Mediated Subpathway Members)
    df     <- rbind(paste(' ', colnames(df), ' ', sep=''), df)
    cnames <- c(rep('',3), rep('Subscores', length(subScoreRange)))

    if ('p-value' %in% names(filters))    
    { 
        cnames <- c(cnames, rep('', 1)) 
    }
    if ('measures' %in% names(filters))    
    { 
        cnames <- c(cnames, rep('Measures', 3)) 
    }

    if (from != to)
    {
        title <- paste0('miRNA Mediated Subpathway Members  (', from, ', ', 
                    to, ')')
    }
    if (from == to)
    {
        title <- paste0('miRNA Mediated Subpathway Members  (', from, ')')
    }
    cnames <- c(cnames, rep(title, genCols), rep('miRNAs', mirCols))
    df     <- rbind(cnames, df)
    wb     <- createWorkbook()
    addWorksheet(wb, 'Results')
    cl     <- lapply(df, function(x) tolower(class(x)))
    wb$writeData(sheet='Results', df=df, startRow=1, startCol=1, 
                colClasses = cl, hlinkNames = c('','SubId',names(links)), 
                colNames=FALSE, keepNA=FALSE) 


    # Customize output
    cols <- c(1:3, subScoreRange)
    if ('p-value' %in% names(filters))   { cols <- c(cols, qvaluesRange) }
    if ('measures' %in% names(filters))  { cols <- c(cols, measuresRange) }
    if ('mirScore' %in% names(filters))  { cols <- c(cols, mirRange) }


    setColWidths(wb, sheet='Results', cols=cols, widths='auto')
    modifyBaseFont(wb, fontSize=10, fontColour='#000000', fontName="Segoe UI")
    freezePane(wb, sheet='Results', firstActiveRow=3, firstActiveCol=2)


    # Main Body
    style <- createStyle(fontSize=10, fontName='Segoe UI', halign="left", 
                        valign = "center", fgFill='#ffffff', 
                        borderColour='#000000')
    addStyle(wb, sheet='Results', style=style, rows=2:(nrow(df)+1), 
                cols=1:ncol(df), gridExpand=TRUE)

    # Gene scores and topologics
    style <- createStyle(fontSize=10, fontName='Segoe UI', halign="center", 
                        valign = "center", fgFill='#ffffff', 
                        borderColour='#000000')
    addStyle(wb, sheet='Results', style=style, rows=2:(nrow(df)+1), 
                        cols=3:ncol(df), gridExpand=TRUE)

    # Gene members and mirnas
    style <- createStyle(fontSize=10, fontName='Segoe UI', halign="left", 
                        valign = "center", fgFill='#ffffff', 
                        borderColour='#000000')

    addStyle(wb, sheet='Results', style=style, rows=3:nrow(df), 
                cols=membersRange, gridExpand=TRUE)

    # Header
    mergeCells(wb, sheet='Results', rows=1, cols=subScoreRange)
    if ('measures' %in% names(filters))          
    {
        mergeCells(wb, sheet='Results', rows=1, cols=measuresRange)
    }

    mergeCells(wb, sheet='Results', rows=1, cols=membersRange)

    if ('mirScore' %in% names(filters))      
    { 
        mergeCells(wb, sheet='Results', rows=1, cols=mirRange)
    }

    style <- createStyle(fontSize=11, fontName='Segoe UI', 
            fontColour='#000000', halign="center", valign = "center", 
            fgFill='#ffffff', borderColour = '#000000')
    addStyle(wb, sheet='Results', style=style, rows=1, cols=1:ncol(df), 
            gridExpand=TRUE)

    # Header Column
    style <- createStyle(fontSize=11, fontName='Segoe UI', halign="left", 
                        valign = "center", fgFill='#ffffff', 
                        borderColour='#000000')
    addStyle(wb, sheet='Results', style=style, rows=1:(nrow(df)+1), cols=1:1, 
                        gridExpand=TRUE)

    setRowHeights(wb, sheet='Results', rows = 1, heights = 35)
    setRowHeights(wb, sheet='Results', rows = 2:(nrow(df)+1), heights = 25)


    # Split gene score matrix to part based on the score
    S          <- df[, subScoreRange, drop=FALSE]
    k         <- length(colors) + 1
    rn         <- 1:(2*k)
    breaks     <- unname(sapply(split(rn,cut(rn,quantile(rn,(0:k)/k), 
                include.lowest=TRUE)), function(x) { x[2] }))
    breaks     <- c(0, breaks/((k-1)*2))
    
    for (i in 2:k)
    {
        idx  <- which((S >= breaks[i-1]) & (S < breaks[i]))
        # Index to row and column
        rows <- ((idx-1) %% nrow(S)) + 1
        cols <- floor((idx-1)/nrow(S)) + 1
        # Set cell colors of gene scores according to their scores
        style <- createStyle(fontSize=10, fontName='Segoe UI', 
                    halign="center", valign = "center", fgFill=colors[i-1])
        addStyle(wb, sheet='Results', style=style, rows=rows, cols=3+cols, 
                    gridExpand=FALSE)
    }

    saveWorkbook(wb, file, overwrite=TRUE)  
}


subpathwayKEGGmap <- function(subpathways, type, openInBrowser)
{
    
    createKEGGLinks(subpathways, type, openInBrowser)
}

createKEGGLinks   <- function(subpaths, type, openInBrowser)
{
    links <- vector(mode='numeric', length=nrow(subpaths))
    for(i in 1:nrow(subpaths))
    {
        if (type == 'Linear')
        {
            links[i] <- createKEGGLink(subpaths[i,,drop=FALSE], 
                        openInBrowser=openInBrowser)
        }
        if (type == 'Non-Linear')
        {
            sub           <- unlist(strsplit(subpaths[i,], split='-'))
            sub           <- matrix(unique(sub), nrow=1)
            rownames(sub) <- rownames(subpaths[i,,drop=FALSE])
            links[i]      <- createKEGGLink(sub, openInBrowser)
        }
    }

    return(links)
}

createKEGGLink    <- function(subpath, openInBrowser)
{
    #
    # Create KEGG map link
    #

    if (missing(openInBrowser)) { openInBrowser <- TRUE}

    pathway <- rownames(subpath)
    prefix  <- paste0('http://www.kegg.jp/kegg-bin/show_pathway?', pathway)
    suffix  <- c()
    for (i in 1:ncol(subpath))
    {
        sf     <- paste(as.numeric(subpath[1,i]), '%09yellow,red')
        suffix <- paste(suffix, sf, sep='/')
    }
    link <- paste(prefix, suffix, sep='')

    if (openInBrowser) { browseURL(link) }
    
    return(link)
}



subpathwayMiRNAs  <- function(summary, subIdx, timePoints)
{
    subpath  <- summary$subpaths[subIdx, , drop=FALSE]
    edgeList <- summary$edgeList
    tt       <- ncol(edgeList) - 2
    thres    <- unname(summary$filters['mirScore'])
    org      <- summary$org
    type     <- summary$type

    if (type == 'Linear')     { baseDir <- cache$dirs$vlnr }
    if (type == 'Non-Linear') { baseDir <- cache$dirs$vnlr }

    filename <- paste0(baseDir, '//subpathway', subIdx, '.pdf')

    if ( missing(timePoints) )
    {
        timePoints <- 1:tt
    }
    if ( !missing(timePoints) )
    {        
        if ( length(which(!timePoints %in% 1:tt)) > 0 )
        {
            message('Invalid time point.')
            return(NULL)
        }
    }

    # Final miRNA-mRNA interactions (miRNA x mRNA)
    sub.edgeList <- edgeList[edgeList[, 2] %in% subpath, ]
    sub.edgeList[, 1] <- gsub(paste0(org,'-'), '', sub.edgeList[, 1])

    mats <- list()
    for ( t in timePoints)
    {
        sub.t.edgeList <- sub.edgeList[, c(1:2, 2+t)]
        idx            <- which(sub.t.edgeList[, 3,  drop=TRUE] > thres)
        sub.t.edgeList <- sub.edgeList[idx, , drop=FALSE]
        sub.miRNAs  <- unique(sub.t.edgeList[, 1])
        sub.mRNAs   <- getSubpathwayGenes(subpath, 'Linear')
        nrows       <- length(sub.miRNAs)
        ncols       <- length(sub.mRNAs)
        sub.edgeMat <- matrix(0, nrow=nrows, ncol=ncols)
        rownames(sub.edgeMat) <- sub.miRNAs
        colnames(sub.edgeMat) <- sub.mRNAs

        for ( i in 1:nrow(sub.t.edgeList) )
        {
            idx1 <- which(rownames(sub.edgeMat) == sub.t.edgeList[i, 1])
            idx2 <- which(colnames(sub.edgeMat) == sub.t.edgeList[i, 2])
            sub.edgeMat[idx1, idx2] <- sub.edgeMat[idx1, idx2] + 1
        }

        # sub.edgeMat <- sub.edgeMat[which(rowSums(sub.edgeMat) > 0), , 
        #                             drop=FALSE]
        # sub.edgeMat <- sub.edgeMat[, which(colSums(sub.edgeMat) > 0), 
        #                             drop=FALSE]

        # Add a null miRNA that targets all genes by default so that mRNAs
        # with no interactions do not have to be removed.
        mat           <- sub.edgeMat * 100
        mat           <- rbind(rep(1, ncol(mat)), mat)
        rownames(mat) <- c(' ', rownames(sub.edgeMat) )
        mats         <- c(mats, list(mat))
    }

    doCirclize(mats, filename)

    return( mats )
}

doCirclize        <- function(mats, filename)
{
    # doCirclize
    tt <- length(mats)
    pdf(filename, width=3*tt, height=3, onefile=TRUE)
    par(mfrow=c(1, tt))

    for ( t in 1:tt )
    {
        mat <- mats[[t]]
        order <- c(rownames(mat), rev(colnames(mat)))
        
        # miRNA colors
        colors        <- c('#ffffff', rainbow(nrow(mat)-1))
        names(colors) <- rownames(mat)

        bgap <- 50

        gd   <- c(rep(2, nrow(mat)-1), bgap, rep(2, ncol(mat)-1), bgap)
        gcol <- NULL
        gcol[rownames(mat)]   <- colors[rownames(mat)]
        gcol[colnames(mat)]   <- 'lightgray'

        # Example
        degree <- 163

        circos.par(start.degree = degree, gap.degree = gd)
        chordDiagram(mat,order=order, grid.col = gcol, transparency = 0.5, 
                    annotationTrack="grid", annotationTrackHeight=0.01, 
                    preAllocateTracks=1)

        # since default text facing in `chordDiagram` is fixed, we need to 
        # manually add text in track 1
        for(si in get.all.sector.index()) 
        {
            xlim = get.cell.meta.data("xlim", sector.index = si, 
                                    track.index = 1)
            ylim = get.cell.meta.data("ylim", sector.index = si, 
                                    track.index = 1)
            circos.text(mean(xlim), ylim[1], si, facing = "clockwise", 
                        adj = c(0, 1.5),
            niceFacing = TRUE, cex = 0.7, col = "black", sector.index = si, 
                        track.index = 1)
        }

        circos.clear() 
    }
    
    dev.off()
}
