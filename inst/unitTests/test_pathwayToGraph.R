test_pathwayToGraph <- function() 
{
	message('Testing test_pathwayToGraph...', appendLF = FALSE)

    sink( tempfile() ) 

    load(system.file('extdata', 'Examples//data.RData', package='CHRONOS'))

    # Create pathway graph
    pathways  <- c('04915', '04917', '04930', '05031')
    graphsRun <- createPathwayGraphs(org='hsa', pathways=pathways)

    sink()

    all.equal(graphsRun, graphs)

    message('done.')
}


