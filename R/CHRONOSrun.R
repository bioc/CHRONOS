
CHRONOSrun <- function(mRNAexp, mRNAlabel, miRNAexp, pathType, subType, 
                        measures, thresholds, org, export, verbose,
                        miRNAinteractions)
{
    if (missing(mRNAexp))    { mRNAexp  <- 'mRNAexpressions.txt' }
    if (missing(mRNAlabel))  { mRNAlabel <- 'entrezgene' }
    if (missing(miRNAexp))   { miRNAexp <- 'miRNAexpressions.txt' }
    if (missing(pathType))   { pathType <- 'All' }
    if (missing(subType))    { subType  <- 'All' }
    if (missing(measures))   { measures <- 'FALSE' }
    if (missing(org))        { org      <- 'hsa' }
    if (missing(export))     { export   <- '.txt' }
    if (missing(verbose))    { verbose  <- FALSE }
    if (missing(thresholds)) 
    {
        thresholds <- c('p-value'=0.04, 'subScore'=0.4, 'mirScore'=0.4)
    }
    if (missing(miRNAinteractions)) { miRNAinteractions <- NULL }

    if (!verbose) { sink( tempfile() ) }

    importExpressions(data=mRNAexp, type=mRNAlabel, sep='\t', org=org)
    importExpressions(data=miRNAexp, org=org, type='miRNA', sep='\t')

    # Download Insulin Signaling Pathway
    paths   <- downloadPathways(org=org, pathways=pathType)
    
    # Create pathway graph
    graphs  <- createPathwayGraphs(org=org, pathways=paths)

    if ( 'Linear' %in% subType || 'All' %in% subType )
    {
        # Extract linear subpathways
        linSubs <- extractLinearSubpathways(graphs=graphs, filter=TRUE)

        # Filter linear subpathways
        linSubsScored <- scoreSubpathways(subpathways=linSubs, 
                            filters=thresholds,
                            miRNAinteractions=miRNAinteractions)
        if ( '.xlsx' %in% export )
        {
            to <- if (org == 'hsa') 'hgnc_symbol' else mRNAlabel
            to <- if (org == 'mmu') 'mgi_symbol' else to
            visualizeResults( summary=linSubsScored, export=export,
                            from=mRNAlabel, to=to ) 
        }
    }
    if ( 'Non-Linear' %in% subType || 'All' %in% subType )
    {
        # Extract non linear subpathways
        nliSubs <- extractNonLinearSubpathways(graphs=graphs, filter=TRUE)

        # Filter non linear subpathways
        nliSubsScored <- scoreSubpathways(subpathways=nliSubs, 
                            filters=thresholds,
                            miRNAinteractions=miRNAinteractions)
        if ( '.xlsx' %in% export )
        {
            to <- if (org == 'hsa') 'hgnc_symbol' else mRNAlabel
            to <- if (org == 'mmu') 'mgi_symbol' else to
            visualizeResults( summary=nliSubsScored, export=export,
                            from=mRNAlabel, to=to )
        }
        
    }
    out <- list('graphs'=graphs)
    if ( 'Linear' %in% subType )
    {
        out <- list('graphs'=graphs, 'linSubs'=linSubs, 
                'linSubsScored'=linSubsScored) 
    }
    if ( 'Non-Linear' %in% subType )
    {
        out <- list('graphs'=graphs, 'nliSubs'=nliSubs, 
                'nliSubsScored'=nliSubsScored)
    }
    if ( 'All' %in% subType )
    {
        out <- list('graphs'=graphs, 'linSubs'=linSubs, 'nliSubs'=nliSubs, 
                'linSubsScored'=linSubsScored, 'nliSubsScored'=nliSubsScored)
    }
    if (!verbose) { sink() }

    return(out)
}


