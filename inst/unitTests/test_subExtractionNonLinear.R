test_subExtractionNonLinear <- function() 
{
    message('Testing test_subExtractionNonLinear...', appendLF = FALSE)

    # Load pathway graphs from toy data
    load(system.file('extdata', 'Examples//data.RData', package='CHRONOS'))

	sink( tempfile() ) 

    # Import mRNA expressions
    mRNAexpr <- importExpressions(data=mRNAexpr, type='entrezgene', org='hsa')

    # Extract linear subpathways
    nliSubsRun  <- extractNonLinearSubpathways(graphs=graphs, filter=TRUE)

    sink()

    all.equal(nliSubsRun, nliSubs)

    message('done.')
}
