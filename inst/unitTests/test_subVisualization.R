test_subVisualization <- function() 
{
    message('Testing test_subVisualization...', appendLF = FALSE)

    load(system.file('extdata', 'Examples//data.RData', package='CHRONOS'))
    
    sink( tempfile() ) 

    importExpressions(data=mRNAexpr, type='entrezgene', org='hsa')
    
    # Export the final subpathways to a xlsx file 
    linSubsVisualRun <- visualizeResults(summary=linSubsScored, export='.xlsx',
                        from='entrezgene', to='hgnc_symbol')
    nliSubsVisualRun <- visualizeResults(summary=nliSubsScored, export='.xlsx',
                        from='entrezgene', to='hgnc_symbol')

    # Visualize a subpathway
    linSubsmiRNAsRun <- subpathwayMiRNAs(summary=linSubsScored, subIdx=2, 
                                        timePoints=3)

    # Opening selected subpathways in default browser
    linLinksRun <- subpathwayKEGGmap(subpathways=linSubs$subpaths[1:2, , drop=FALSE], 
                                    type='Linear', openInBrowser=FALSE)

    sink()

    all.equal(linSubsVisualRun[, -9], linSubsVisual[, -9])
    all.equal(linSubsmiRNAsRun, linSubsmiRNAs)
    all.equal(linLinksRun, linLinks)

    message('done.')
}
