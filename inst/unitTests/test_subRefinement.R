test_subRefinement <- function() 
{
    message('Testing test_subRefinement...', appendLF = FALSE)

    # Load extracted subpathways from toy data
    load(system.file('extdata', 'Examples//data.RData', package='CHRONOS'))

    sink( tempfile() ) 

    # Import mRNA expressions
    mRNAexpr  <- importExpressions(data=mRNAexpr, type='entrezgene', org='hsa')
    # Import miRNA expressions
    miRNAexpr <- importExpressions(data=miRNAexpr, type='miRNA', org='hsa')

    # Score extracted linear subpathways
    filters          <- c('subScore'=0.4, 'mirScore'=0.4)
    linSubsScoredRun <- scoreSubpathways(subpathways=linSubs, filters=filters,
                        miRNAinteractions=miRNAinteractions)
    all.equal(linSubsScoredRun, linSubsScored)    

    sink()

    message('done.')
}
