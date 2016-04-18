cache <- new.env()

.onLoad            <- function(libname, pkgname) 
{
    #
    #
}

.onAttach          <- function(libname, pkgname)
{
    # Create ext directory internal structure.
    .buildDirectories()
}

.onUnload          <- function(libname) 
{

    # 
}

.buildDirectories  <- function()
{
    # The function creates data directories for user usage.
    #
    # extdata
    #    -- Downloads
    #       -- KEGG
    #          -- hsa
    #             + hsa00010.xml, hsa00020.xml ...
    #          -- mmu
    #                   + mmu00010.xml, mmu00020.xml ...
    #       -- miRecords
    #          -- hsa
    #          -- mmu
    #    -- Input
    #       -- Default
    #          - mirnaTargets.txt
    #       -- User
    #               + geneExpressions.txt
    #               + mirnaExpressions.txt
    #    -- Output
    #       -- Subpaths
    #          -- Linear
    #             -- hsa
    #                + subpaths1.RData
    #                + subpaths1.txt
    #                + subpaths2.txt
    #             -- mmu
    #                + subpaths1.RData
    #                + subpaths1.txt
    #          -- Non Linear
    #             -- hsa
    #                + subpaths1.txt
    #             -- mmu
    #                + subpaths1.RData
    #       -- Scores
    #          -- Linear
    #             -- hsa
    #                + score1a.xlsx
    #                + score1b.xlsx
    #                + score2a.xlsx
    #                + score2b.xlsx
    #             -- mmu
    #                + score1a.xlsx
    #                + score1b.xlsx
    #          -- Non Linear
    #             -- hsa
    #                + score1a.xlsx
    #                + score1b.xlsx
    #             -- mmu
    #                + score1a.xlsx
    #                + score1b.xlsx
    #                + score2a.xlsx
    #                + score2b.xlsx
    dirs  <- list()

    path <- switch(.Platform$OS.type, unix = path.expand("~"),
                    windows= file.path(gsub("\\\\", "/",
                    Sys.getenv("USERPROFILE")), "AppData"))
    opt <- getOption("CHRONOS_CACHE", file.path(path, '.CHRONOS'))
    chronosDir <- Sys.getenv("CHRONOS_CACHE", opt)

    # Directory for user files
    inDir <- paste(chronosDir, '//extdata//User',sep='')
    dirs  <- c(dirs, list(usr=inDir))
    msg0 <- sapply(list(inDir), function(x) 
        { dir.create(x, showWarnings=FALSE, recursive=TRUE)} )    

    if ( is.null(dirs$geneExpressions) )
    {
        gf1   <- paste(inDir, 'mRNAExpressions.txt', sep='//')
        dirs  <- c(dirs, list(geneExpressions=gf1) )
    }
    if ( is.null(dirs$mirnaExpressions) )
    {
        mf1   <- paste(inDir, 'miRNAExpressions.txt', sep='//')
        dirs  <- c(dirs, list(mirnaExpressions=mf1) )
    }
    if ( is.null(dirs$interestingGenes) )
    {
        igf   <- paste(inDir, 'interestingGenes.txt', sep='//')
        dirs  <- c(dirs, list(interestingGenes=igf) )
    }
    if ( is.null(dirs$miRecordsFile) )
    {
        mrc   <- paste(inDir, 'miRecords.xlsx', sep='//')
        dirs  <- c(dirs, list(miRecordsFile=mrc) )
    }
    if ( is.null(dirs$tarBaseFile) )
    {
        trb   <- paste(inDir, 'TarBase.xlsx', sep='//')
        dirs  <- c(dirs, list(tarBaseFile=trb) )
    }

    # Directory for downloaded data
    pDir    <- paste(chronosDir, '//extdata', sep='')
    xmlDir  <- paste(pDir, '//Downloads//KEGG', sep = '')
    miDir   <- paste(pDir, '//Downloads//miRecords', sep = '')
    anDir   <- paste(pDir, '//Downloads//mirbase', sep = '')
    msg1    <- sapply(list(xmlDir, miDir, anDir), function(x) 
        { dir.create(x, showWarnings=FALSE, recursive=TRUE)} )
    dirs    <- c(dirs, list(xml=xmlDir, miDir=miDir, anDir=anDir,
                            miFile='miRNATargets.RData'))

    # Directory for temporary files
    intDir  <- paste(chronosDir, '//extdata',
                '//Internal', sep='')
    matDir  <- paste(intDir, '//mat', sep='') #
    idsDir  <- paste(intDir, '//ids', sep='') #
    parDir  <- paste(intDir, '//par', sep='') #
    subDir  <- paste(intDir, '//out', sep='') #
    jDir    <- paste(chronosDir, '//java', sep='')
    tmpDirs <- list(mat=matDir,ids=idsDir,par=parDir, sub=subDir)
    dirs    <- c(dirs, list(int=intDir, tmp=tmpDirs,java=jDir))

    #
    # Directories for output
    # 
    outDir  <- paste0(chronosDir, '//extdata','//Output')
    # Subpaths
    nlDir   <- paste(outDir, '//Subpaths//Non Linear', sep = '')
    lnDir   <- paste(outDir, '//Subpaths//Linear', sep='')
    dirs    <- c(dirs, list(lnr=lnDir,nlr=nlDir))
    msg2    <- sapply(list(nlDir, lnDir), function(x) 
        { dir.create(x, showWarnings=FALSE, recursive=TRUE)} )

    # Scores
    snlDir  <- paste(outDir, '//Scores//Non Linear', sep = '')
    slnDir  <- paste(outDir, '//Scores//Linear', sep='')
    dirs    <- c(dirs, list(slnr=slnDir, snlr=snlDir))
    msg2    <- sapply(list(snlDir, slnDir), function(x) 
        { dir.create(x, showWarnings=FALSE, recursive=TRUE)} )

    # Visualizations
    vnlDir  <- paste(outDir, '//Visualizations//Non Linear', sep = '')
    vlnDir  <- paste(outDir, '//Visualizations//Linear', sep='')
    dirs    <- c(dirs, list(vlnr=vlnDir, vnlr=vnlDir))
    msg3    <- sapply(list(vnlDir, vlnDir), function(x) 
        { dir.create(x, showWarnings=FALSE, recursive=TRUE)} )

    # Find full path of java executable
    dirs$javapath <- shQuote(.Sys.which2('java'))

    cache$dirs <- dirs

    # Copy demo files from package directory to user directort
    file.copy(  from=system.file('extdata', package='CHRONOS'), 
                to=file.path(chronosDir), 
                overwrite = TRUE, recursive = TRUE, copy.mode = TRUE)

    return(dirs)
}

.createDirectories <- function(org)
{
    #
    # Create subdirectories for specific organism inside each parent folder
    #
    d1 <- paste(cache$dirs$tmp$mat, org, sep='//') 
    d2 <- paste(cache$dirs$tmp$ids, org, sep='//') 
    d3 <- paste(cache$dirs$tmp$par, org, sep='//') 
    d4 <- paste(cache$dirs$tmp$sub, org, sep='//') 
    d5 <- paste(cache$dirs$lnr, org, sep='//')
    d6 <- paste(cache$dirs$nlr, org, sep='//')
    d7 <- paste(cache$dirs$slnr, org, sep='//')
    d8 <- paste(cache$dirs$snlr, org, sep='//')

    mapply(dir.create, list(d1, d2, d3, d4, d5, d6, d7, d8), 
            showWarnings=FALSE, rec=TRUE)
}

.cleanDirectories  <- function()
{
    #
    unlink(cache$dirs$int, recursive=TRUE)
}

