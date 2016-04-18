test_subExtractionLinear <- function() 
{
    message('Testing test_subExtractionLinear...', appendLF = FALSE)

    # Load pathway graphs from toy data
    load(system.file('extdata', 'Examples//data.RData', package='CHRONOS'))

    sink( tempfile() ) 

    # Import mRNA expressions
    mRNAexpr <- importExpressions(data=mRNAexpr, type='entrezgene', org='hsa')
    # Extract linear subpathways
    linSubsRun <- extractLinearSubpathways(graphs=graphs, filter=TRUE)
    
    sink()

    all.equal(linSubsRun, linSubs)

    message('done.')
}
