test_CHRONOSrun <- function() 
{
	message('Testing test_CHRONOSrun...', appendLF = FALSE)
    
    sink( tempfile() ) 

    load(system.file('extdata', 'Examples//data.RData', package='CHRONOS'))


    out <- CHRONOSrun( 
                mRNAexp=mRNAexpr,
    			mRNAlabel='entrezgene',
    			miRNAexp=miRNAexpr,
    			pathType=c('04915', '04917', '04930', '05031'),
				org='hsa',
				subType='All',
				thresholds=c('subScore'=0.4, 'mirScore'=0.4),
                miRNAinteractions=miRNAinteractions)
    
    sink()

    all.equal(CHRONOSdemo, out)

    message('done.')
}

