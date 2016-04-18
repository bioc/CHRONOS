

downloadMiRecords     <- function(org, pn, update, databases)
{
    #
    # Downloads mirna data from miRecords and summarizes all organism 
    # interactions between mirnas and gene targets. 
    #

    message('Downloading miRecords data...', appendLF = FALSE)
    if (missing(update))    { update    <- FALSE }
    if (missing(databases)) { databases <- 'All' }
    if (missing(pn))        { pn        <- 5 }


    # Create the organism download directory and file
    dir  <- cache$dirs$miDir
    file <- cache$dirs$miFile
    dir  <- paste(dir, org, sep='//')
    file <- paste(dir, file, sep='//')
    if (!file.exists(dir))  { dir.create(dir, recursive=TRUE) }
    if (!file.exists(file)) { file.create(file) }

    # Create full records and store them locally 
    targets <- tryCatch( 
    { 
        if (!update)
        {
            # Download from remote server
            df <- .getProcessedData(pn=pn, file=file, databases=databases,
                org=org)
        }
        if (update)
        {
            miRecordsFile <- cache$dirs$miRecordsFile
            tarBaseFile   <- cache$dirs$tarBaseFile

            df1 <- .getValMiRecordsTargets(miRecordsFile=miRecordsFile, 
                                            org=org)

            df2 <- .getValTarbaseTargets(tarBaseFile=tarBaseFile, org=org)

            df3 <- .getMiRecordsData(org=org, pn=pn, k=24)
            
            # Filter predicted targets according to selected databases.
            df3 <- .filterRecords(targets=df3, databases=databases)

            df  <- rbind(df1, df2, df3)
            df  <- df[!duplicated(df[,1:2]), ]
        }
        targets <- df
    }, 
        error = function(e) 
        {
            df <- .getProcessedData(pn=pn, file=file, mode='error', 
                    databases=databases, org=org)
        }
    ) 

    save(file=file, targets, compress='xz')

    message('done.')

    return(targets)
}


#
# Download miRNA predicted targets from miRecords
#
.getMiRecordsData     <- function(org, pn, k)
{
    # Create and register parallel backend
    cl <- makeCluster(detectCores(logical=FALSE))
    on.exit(stopCluster(cl), add=TRUE)
    registerDoParallel(cl)

    # Get the list of all availiable mirnas for organism
    lines  <- .getMiRecordsMiRNAs(org)

    # Current http links for retrieving KEGG xml file
    mirror <- 'http://c1.accurascience.com/miRecords'
    prefix <- paste(mirror, '/download_data.php?p=1', sep='') 
    suffix <- '&targetgene_type=refseq_acc'
    suffix <- paste0(suffix, '&targetgene_info=&search_int=Search&pn=')
    
    # Create download links 
    species <- paste('&species=', org, sep='') 
    mirnas  <- paste('&mirna_acc=', lines, sep='')
    links   <- paste(prefix, species, mirnas, suffix, pn, sep='')    

    # Download data 
    S <- splitWork(N=length(lines), k=k)
    i <- 0
    R <- foreach (i = 1:S$k, .combine='rbind') %do%
    {
        rawText <- getURL(links[S$jobs[[i]]], 
                            .opts = list(timeout = 20, verbose=FALSE))
        Sys.sleep(1)
        # Data cleanup
        Rt <- matrix(nrow=0, ncol=14)
        for (r in rawText)
        {
            r1 <- as.list(unlist(strsplit(r, '\n'))[-c(1,2)])
            r2 <- sapply(r1, function(x)   { strsplit(x, '\t') })
            rt <- t(sapply(r2, function(x) { x[-c(2,5)] }))

            if(ncol(rt) == ncol(Rt)) { Rt <- rbind(Rt, rt) }
        }

        return(Rt)
    }
    if ( nrow(R) > 0 ) { rownames(R) <- 1:nrow(R) } 
    
    db <- c('diana', 'microinspector', 'miranda', 'mirtarget2', 'mitarget', 
            'nbmirtar', 'pictar', 'pita', 'rna22', 'rnahybrid', 'targetscan')
    colnames(R) <- c('mirna', 'Refseq', 'GeneName', db)
    R         <- as.data.frame(R, stringsAsFactors=FALSE)
    
    if ( nrow(R) == 0 )
    {
        colnames(R)[2] <- 'Entrez'; return(R)
    }

    #
    # Annotate genes from NCBI Refseq to Entrez Gene Accession
    #
    targets <- R[, 2]
    
    # Annotate targets to entrez ids.
    annMap <- convertNomenclature(ids=targets, org=org, from='refseq_mrna', 
                                    to='entrezgene')

    # Annotation answer from ensembl server contains unique matches, 
    # whereas a gene can be the target of multiple miRNAs. 
    # Create a hash to map distinct refseq ids to entrez ids.
    lib1        <- annMap[, 2]
    names(lib1) <- annMap[, 1]
    hitIdx      <- !is.na(lib1[targets])
    R           <- R[hitIdx, ]
    targets     <- targets[hitIdx]

    # Map all original refseq targets to entrez
    R[,2]          <- lib1[targets]
    colnames(R)[2] <- 'Entrez'

    if ( nrow(R) == 0 ) { return(R) }

    # Remove duplicates
    idx <- which(duplicated(R[, 2:3]) )
    if (length(idx) > 0) 
    { 
        R       <- R[-idx, , drop=FALSE]
        targets <- targets[-idx]
    }
        
    # Remove interactions between miRNAs and genes with no entrez annotation
    idx <- which(R[,2] == targets)
    if (length(idx) > 0) { R <- R[-idx, , drop=FALSE] }

    if ( nrow(R) > 0 )
    {
        R           <- R[!duplicated(R[,1:2]), ]
        rownames(R) <- 1:nrow(R)
    }

    return(R)
}


.getMiRecordsMiRNAs  <- function(org)
{
    #
    # Downloads HTML code and parses it for all availiable mirnas of organism
    #
    link    <- 'http://c1.accurascience.com/miRecords/prediction_query.php'
    rawText <- getURL(link)
    
    # Extract mirnas for each organism 
    R   <- strsplit(rawText, 'species.value==')
    R   <- lapply(R, function(x) { unlist(strsplit(x, 'else')) })
    R   <- sapply(R, function(x) { x[-seq(1,length(x)-1,2)] })
    R   <- sapply(R, function(x) { unlist(strsplit(x, '[()]')) } )
    orgs <- sapply(R, function(x) { gsub('\'','',x[1]) } )
    R   <- sapply(R, function(x) { x[seq(1,length(x)-1,2)] })
    R   <- sapply(R, function(x) { gsub('\'','',x) })
    R   <- sapply(R, function(x) { unlist(strsplit(x, ',')) } )
    R   <- sapply(R, function(x) { x[seq(3,length(x)-1,2)] })
    R   <- R[-length(R)]

    names(R) <- sapply(strsplit(orgs[-length(orgs)],' '), function(x) 
                { 
                    tolower(paste(substr(x[1],1,1), substr(x[2],1,2), sep=''))
                })


    # If an organism has been specified, return only corresponding mirnas.
    if (!org %in% names(R)) 
    { 
        message(org, 'mirnas not availiable in miRecords.')
        return(NULL) 
    }
    if (!missing(org))
    { 
        R <- R[[org]] 
    }

    return(R)
}


#
# Download miRNA validated targets from miRecords
#
.getValTarbaseTargets   <- function(tarBaseFile, org)
{
    data    <- read.xlsx(xlsxFile=tarBaseFile)
    data    <- data[, c('Organism', 'miRNA', 'Gene', 'Ensembl')]
    orgName <- ''

    if (org == 'hsa')   { orgName <- 'Human' }
    if (org == 'mmu')   { orgName <- 'Mouse' }
    if (org == 'rno')   { orgName <- 'Rat' }
    if (orgName == '')  { return(NULL) }

    # Keep only organism specific interactions
    idx       <- which(data[,1] == orgName)
    if (length(idx) > 0) { data <- data[idx, -1, drop=FALSE] }

    df        <- cbind(data[,1], data[,3], data[,2])
    summary   <- matrix(0, nrow=nrow(df), ncol=11)
    R         <- cbind(df, summary)
    db        <- c('diana', 'microinspector', 'miranda', 'mirtarget2', 
                        'mitarget', 'nbmirtar', 'pictar', 'pita', 'rna22', 
                        'rnahybrid', 'targetscan')
    colnames(R)  <- c('mirna',  'Ensembl', 'GeneName', db)
    rownames(R)  <- 1:nrow(R)
    R            <- as.data.frame(R, stringsAsFactors=FALSE)
    R[,'diana']  <- 1

    # Convert miRNA names to the ones currently used by miRecords
    mirOld    <- paste(org, '-', unique(R[,1]), sep='')
    mirConv   <- convertMiRNANomenclature(org=org, miRNAs=mirOld)

    # Find miRNAs with no conversion
    idx       <- which(is.na(mirConv[,2]))
    if (length(idx) > 0)
    {
        mirDeleted  <- mirConv[idx, 1]
        mirDeleted  <- gsub(paste(org, '-', sep=''), '', mirDeleted)
        # Remove these miRNAs from results
        idx   <- which(R[,1] %in% mirDeleted)
        if (length(idx) > 0)
        { 
            R <- R[-idx, ,drop=FALSE] 
        }
        if (nrow(R) == 0)     { return(R) }
    }

    # Remove interactions with no ensembl targets
    idx   <- which(R[,2] == 'n_a')
    if (length(idx) > 0)
    {
        R <- R[-idx, ,drop=FALSE]
    }
    if (nrow(R) == 0) { return(R) }

    R[, 1]      <- paste(org, '-', R[,1], sep='') 
    rownames(R)  <- 1:nrow(R)

    #
    # Annotate genes from GeneName to Entrez Gene Accession
    #

    targets <- R[,2]

    # Annotate targets to entrez ids.
    annMap <- convertNomenclature(ids=targets, org=org, 
                        from='ensembl_gene_id', to='entrezgene')
    if (class(annMap) != 'matrix') 
    { 
        return(NULL) 
    }

    # Annotation answer from ensembl server contains unique matches, 
    # whereas a gene can be the target of multiple miRNAs. 
    # Create a hash to map distinct refseq ids to entrez ids.
    lib1            <- annMap[,2]
    names(lib1)     <- annMap[,1]
    # Map all original refseq targets to entrez
    R[, 2]          <- lib1[targets]
    colnames(R)[2]  <- 'Entrez'
    hitIdx          <- !is.na(R[, 2])
    R               <- R[hitIdx, ,drop=FALSE]
    targets         <- targets[hitIdx]

    # Remove duplicates
    idx       <- duplicated(R[, 2:3]) 
    if (length(idx) > 0) 
    { 
        R       <- R[-idx, , drop=FALSE]
        targets <- targets[-idx]
    }
    
    # Remove interactions between miRNAs and genes with no entrez annotation
    idx         <- which(R[,2] == targets)
    if (length(idx) > 0) { R <- R[-idx, , drop=FALSE] }
    R           <- R[!duplicated(R[,1:2]), ]
    rownames(R) <- 1:nrow(R)

    # Keep only relevant information
    R           <- R[, 1:3, drop=FALSE]

    return(R)
}


#
# Download miRNA validated targets from TarBase
#
.getValMiRecordsTargets <- function(miRecordsFile, org)
{
    data     <- read.xlsx(xlsxFile=miRecordsFile)
    selCols  <- c('Target.gene_species_scientific', 'miRNA_mature_ID', 
                    'Target.gene_Refseq_acc', 'Target.gene_name')
    data        <- data[, selCols]

    orgNames    <- unlist(sapply(unique(data[,1]), function(x) 
    {
        data    <- unlist(strsplit(x, ' ') )
        res  <- NULL
        if (length(data) == 2)
        {
            res <- tolower(paste(substr(data[1],1,1), substr(data[2],1,2), 
                            sep=''))
        }
    }))

    orgName <- names(orgNames)[which(orgNames == org)]
    idx <- which(data[, 1] == orgName)
    if (length(idx) == 0) { return(NULL) }
    df  <- data[idx, -1, drop=FALSE]

    summary <- matrix(-1, nrow=nrow(df), ncol=11)
    R       <- cbind(df, summary)
    db      <- c('diana', 'microinspector', 'miranda', 'mirtarget2', 
                'mitarget', 'nbmirtar', 'pictar', 'pita', 'rna22', 
                'rnahybrid', 'targetscan')
    colnames(R)  <- c('mirna',  'Ensembl', 'GeneName', db)
    rownames(R)  <- 1:nrow(R)
    R            <- as.data.frame(R, stringsAsFactors=FALSE)

    # Convert miRNA names to the ones currently used by miRecords
    mirOld    <- unique(R[,1])
    mirConv   <- convertMiRNANomenclature(org=org, miRNAs=mirOld)

    # Find miRNAs with no conversion
    idx       <- which(is.na(mirConv[,2]))

    if (length(idx) > 0)
    {
        mirDeleted  <- mirConv[idx, 1]
        # Remove these miRNAs from results
        idx   <- which(R[,1] %in% mirDeleted)
        if (length(idx) > 0)    { R <- R[-idx, ,drop=FALSE] }
        if (nrow(R) == 0)     { return(R) }
    }

    #
    # Annotate genes from GeneName to Entrez Gene Accession
    #

    targets   <- R[,2]

    # Annotate targets to entrez ids.
    annMap      <- convertNomenclature(ids=targets, org=org, 
                        from='refseq_mrna', to='entrezgene')
    if (class(annMap) != 'matrix') 
    { 
        return(NULL) 
    }

    # Annotation answer from ensembl server contains unique matches, 
    # whereas a gene can be the target of multiple miRNAs. 
    # Create a hash to map distinct refseq ids to entrez ids.
    lib1        <- annMap[,2]
    names(lib1) <- annMap[,1]
    # Map all original refseq targets to entrez
    R[, 2]      <- lib1[targets]
    colnames(R)[2] <- 'Entrez'
    hitIdx  <- !is.na(R[, 2])
    R       <- R[hitIdx, ,drop=FALSE]
    targets <- targets[hitIdx]

    # Remove duplicates
    idx       <- duplicated(R[, 2:3]) 
    if (length(idx) > 0) 
    { 
        R     <- R[-idx, , drop=FALSE]
        targets <- targets[-idx]
    }
    
    # Remove interactions between miRNAs and genes with no entrez annotation
    idx       <- which(R[,2] == targets)
    if (length(idx) > 0) { R <- R[-idx, , drop=FALSE] }

    R             <- R[!duplicated(R[,1:2]), ]
    rownames(R)  <- 1:nrow(R)
    
    # Keep only relevant information
    R             <- R[, 1:3, drop=FALSE]
    
    return(R)
}



#
# Filter downloaded records keeping a subset of targets for a set of 
# selected databases
#
.filterRecords     <- function(targets, databases)
{
    if (missing(databases) || databases=='All')
    {
        targets <- targets[, 1:3, drop=FALSE]
        return(targets)
    }
    
    # Databases indexes belong in 1:11
    if (length(which((!unique(databases) %in% 1:11))))
    {
        message('Please supply databases indexes from 1 to 11.')
        return(NULL)
    }

    # Targets are outputs of at least 'pn' databases. 
    # Databases argument can further limit the results.
    targets <- targets[, c(1:3, 3 + databases)]
    
    # Remove targets with no matches
    data     <- matrix(as.integer(as.matrix(targets[, -(1:3)])), 
                        nrow=nrow(targets)) 
    idx  <- which(rowSums( data ) > 0)
    targets <- targets[idx,  1:3, drop=FALSE]
    targets <- as.data.frame(targets, stringsAsFactors=FALSE)
    

    return(targets)
}



#
# Download backup data from server
#
.getProcessedData   <- function(pn, file, mode, databases, org) 
{
    if (missing(mode)) { mode <- '' }
    error1 <- sprintf('\nUnable to connect to miRecords, 
                        downloading backup data.')
    error2 <- sprintf('\nError while downloading data.\n')
    error3 <- sprintf('\nUnable to download backup miRNA targets file.')

    # Downloading from server due to failure to connect to miRecords
    if (mode=='error') { message(error1) }

    # Try downloading from server
    url <- .getMiRNAtargetsUrl(pn=pn, org=org)

    con <- file(file, open = "wb")
    tryCatch(writeBin(getBinaryURL(url, ssl.verifypeer = FALSE), con),
        error = function(e) { message(error3) })
    close(con)
    e=new.env(); tryCatch(load(file=file, envir=e), error=function(e) {   })
    df <- e$targets

    if (!is.null(df))
    {
        # Filter targets according to selected databases.
        df <- .filterRecords(targets=df, databases=databases)
    }

    return(df)
}

.getMiRNAtargetsUrl <- function(pn, org) 
{
    orgs <- c('cel', 'cfa', 'dme', 'dre', 'gga', 'hsa', 'mmu', 'oar', 'rno')

    miKeys  <- vector(mode='list', length=6)
    if (org =='cel')
    {
        miKeys[[1]]  <- c('lk3ws6k5269ihn1', 'miRNAtargets_cel_1.RData')
        miKeys[[2]]  <- c('8tqmr5mh1xqenvd', 'miRNAtargets_cel_2.RData')
        miKeys[[3]]  <- c('z52jmlrwkenwcjf', 'miRNAtargets_cel_3.RData')
        miKeys[[4]]  <- c('8yxktk8cwunnpgx', 'miRNAtargets_cel_4.RData')
        miKeys[[5]]  <- c('pa1u4j81ghcxui1', 'miRNAtargets_cel_5.RData')
        miKeys[[6]]  <- c('f8vczynapddlcrq', 'miRNAtargets_cel_6.RData')
        miKeys[[7]]  <- c('3ns58ooy2p7zefm', 'miRNAtargets_cel_7.RData')
        miKeys[[8]]  <- c('qj46i938a8g7osk', 'miRNAtargets_cel_8.RData')
        miKeys[[9]]  <- c('59ivtx2l7bxfu2u', 'miRNAtargets_cel_9.RData')
        miKeys[[10]] <- c('c75wqk7dllcaysx', 'miRNAtargets_cel_10.RData')
        miKeys[[11]] <- c('q6bvxfqrz6bzzqp', 'miRNAtargets_cel_11.RData')
    }
    if (org =='cfa')
    {
        miKeys[[1]]  <- c('eevs9uh2dz8t4ll', 'miRNAtargets_cfa_1.RData')
        miKeys[[2]]  <- c('rt5ye1bzqap457t', 'miRNAtargets_cfa_2.RData')
        miKeys[[3]]  <- c('ooqbjpfcsavelhj', 'miRNAtargets_cfa_3.RData')
        miKeys[[4]]  <- c('gvi2fi0880lstb6', 'miRNAtargets_cfa_4.RData')
        miKeys[[5]]  <- c('xkewhltcouk24uo', 'miRNAtargets_cfa_5.RData')
        miKeys[[6]]  <- c('nhyz1438ei1flzl', 'miRNAtargets_cfa_6.RData')
        miKeys[[7]]  <- c('zvnz77xyv9tm8mi', 'miRNAtargets_cfa_7.RData')
        miKeys[[8]]  <- c('yxjxwzu7ln30ie3', 'miRNAtargets_cfa_8.RData')
        miKeys[[9]]  <- c('0qogawcoeon87rd', 'miRNAtargets_cfa_9.RData')
        miKeys[[10]] <- c('m7lyjft51xhlppy', 'miRNAtargets_cfa_10.RData')
        miKeys[[11]] <- c('rnaz7476wyq58tf', 'miRNAtargets_cfa_11.RData')
    }
    if (org =='dme')
    {
        miKeys[[1]]  <- c('bv98y2ou5qxtn7k', 'miRNAtargets_dme_1.RData')
        miKeys[[2]]  <- c('unis5d46y8fo733', 'miRNAtargets_dme_2.RData')
        miKeys[[3]]  <- c('fd2xd4j4ph7wlvm', 'miRNAtargets_dme_3.RData')
        miKeys[[4]]  <- c('gf6jnhysfdiqwlu', 'miRNAtargets_dme_4.RData')
        miKeys[[5]]  <- c('d9tx12902ybfevc', 'miRNAtargets_dme_5.RData')
        miKeys[[6]]  <- c('h5kxhosghasg02d', 'miRNAtargets_dme_6.RData')
        miKeys[[7]]  <- c('qsrypy5xfywsd8t', 'miRNAtargets_dme_7.RData')
        miKeys[[8]]  <- c('yvtjpviu8995yp7', 'miRNAtargets_dme_8.RData')
        miKeys[[9]]  <- c('j59am4o0z7uaza5', 'miRNAtargets_dme_9.RData')
        miKeys[[10]] <- c('f0hrkye8aab0y2h', 'miRNAtargets_dme_10.RData')
        miKeys[[11]] <- c('86qt9wn5nbwnt6q', 'miRNAtargets_dme_11.RData')
    }
    if (org =='dre')
    {
        miKeys[[1]]  <- c('b20bwh7wddo328l', 'miRNAtargets_dre_1.RData')
        miKeys[[2]]  <- c('580gu82122yxm7k', 'miRNAtargets_dre_2.RData')
        miKeys[[3]]  <- c('t87acvrksydoqz7', 'miRNAtargets_dre_3.RData')
        miKeys[[4]]  <- c('gtxigfkrz7alav6', 'miRNAtargets_dre_4.RData')
        miKeys[[5]]  <- c('1cf0abkkhogfhi5', 'miRNAtargets_dre_5.RData')
        miKeys[[6]]  <- c('sch97a78dhtzl5j', 'miRNAtargets_dre_6.RData')
        miKeys[[7]]  <- c('vi98uvq80gi5ank', 'miRNAtargets_dre_7.RData')
        miKeys[[8]]  <- c('ik5w9mekugadnfj', 'miRNAtargets_dre_8.RData')
        miKeys[[9]]  <- c('d3c9dlnq0jy7223', 'miRNAtargets_dre_9.RData')
        miKeys[[10]] <- c('gye1kn4zjak6ewr', 'miRNAtargets_dre_10.RData')
        miKeys[[11]] <- c('rfazqoxw6pv7mqo', 'miRNAtargets_dre_11.RData')
    }
    if (org =='gga')
    {
        miKeys[[1]]  <- c('hcnmq8jdll2p3p3', 'miRNAtargets_gga_1.RData')
        miKeys[[2]]  <- c('avcysl4hvqj5cxp', 'miRNAtargets_gga_2.RData')
        miKeys[[3]]  <- c('qplewyvve72kvhj', 'miRNAtargets_gga_3.RData')
        miKeys[[4]]  <- c('0f2dqnj1ehagfgv', 'miRNAtargets_gga_4.RData')
        miKeys[[5]]  <- c('3sjyofi68uxjnpk', 'miRNAtargets_gga_5.RData')
        miKeys[[6]]  <- c('a01hqd5rs6lxwfj', 'miRNAtargets_gga_6.RData')
        miKeys[[7]]  <- c('dqqusm7sc8wtnfm', 'miRNAtargets_gga_7.RData')
        miKeys[[8]]  <- c('bcthqigqtbzqqi1', 'miRNAtargets_gga_8.RData')
        miKeys[[9]]  <- c('yqtx0w9lpcsvyrh', 'miRNAtargets_gga_9.RData')
        miKeys[[10]] <- c('c2xup5ccsfcidpi', 'miRNAtargets_gga_10.RData')
        miKeys[[11]] <- c('qjqc7yypf8tkf17', 'miRNAtargets_gga_11.RData')
    }
    if (org =='hsa')
    {
        miKeys[[1]]  <- c('zf7efg2f241sf7g', 'miRNAtargets_hsa_1.RData')
        miKeys[[2]]  <- c('2k82cbvm8zlwimb', 'miRNAtargets_hsa_2.RData')
        miKeys[[3]]  <- c('4ftt3afronh6zfz', 'miRNAtargets_hsa_3.RData')
        miKeys[[4]]  <- c('gsy14uirkwqknst', 'miRNAtargets_hsa_4.RData')
        miKeys[[5]]  <- c('nordc6iniokqfur', 'miRNAtargets_hsa_5.RData')
        miKeys[[6]]  <- c('wkda1rd9vf1ta0v', 'miRNAtargets_hsa_6.RData')
        miKeys[[7]]  <- c('9xvan8qky7dvjg7', 'miRNAtargets_hsa_7.RData')
        miKeys[[8]]  <- c('owlnyhz47972chw', 'miRNAtargets_hsa_8.RData')
        miKeys[[9]]  <- c('ec1yuhh9nnvj6mu', 'miRNAtargets_hsa_9.RData')
        miKeys[[10]] <- c('4ksggx0q2g97v2q', 'miRNAtargets_hsa_10.RData')
        miKeys[[11]] <- c('4uxnid3aeyxk85f', 'miRNAtargets_hsa_11.RData')
    }
    if (org =='mmu')
    {
        miKeys[[1]]  <- c('heaaqjfym8fhxli', 'miRNAtargets_mmu_1.RData')
        miKeys[[2]]  <- c('vjyogjecevvky5j', 'miRNAtargets_mmu_2.RData')
        miKeys[[3]]  <- c('p12p30cxdvia7qu', 'miRNAtargets_mmu_3.RData')
        miKeys[[4]]  <- c('3402paa4mda1thk', 'miRNAtargets_mmu_4.RData')
        miKeys[[5]]  <- c('qh0eapmat0vb5jv', 'miRNAtargets_mmu_5.RData')
        miKeys[[6]]  <- c('sxyclzohrrryd2i', 'miRNAtargets_mmu_6.RData')
        miKeys[[7]]  <- c('318ilma6p5v4hk6', 'miRNAtargets_mmu_7.RData')
        miKeys[[8]]  <- c('gg9jt8tcxzrswza', 'miRNAtargets_mmu_8.RData')
        miKeys[[9]]  <- c('i2qx257s3dw677e', 'miRNAtargets_mmu_9.RData')
        miKeys[[10]] <- c('0pr8gg5ksa7k94a', 'miRNAtargets_mmu_10.RData')
        miKeys[[11]] <- c('bf8jm388gvg3ssa', 'miRNAtargets_mmu_11.RData')
    }
    if (org =='oar')
    {
        miKeys[[1]]  <- c('h53xpu6u1o96rav', 'miRNAtargets_oar_1.RData')
        miKeys[[2]]  <- c('thvrfwsx3pggmrz', 'miRNAtargets_oar_2.RData')
        miKeys[[3]]  <- c('wc4eullh495ame7', 'miRNAtargets_oar_3.RData')
        miKeys[[4]]  <- c('ixtu0vwii6ncy6w', 'miRNAtargets_oar_4.RData')
        miKeys[[5]]  <- c('wf6paifn14l5ofz', 'miRNAtargets_oar_5.RData')
        miKeys[[6]]  <- c('ar9z5meiwdyxe08', 'miRNAtargets_oar_6.RData')
        miKeys[[7]]  <- c('a74qewl9pxujb3l', 'miRNAtargets_oar_7.RData')
        miKeys[[8]]  <- c('54viw7lha6oklvb', 'miRNAtargets_oar_8.RData')
        miKeys[[9]]  <- c('x9twcfnov3nuz3m', 'miRNAtargets_oar_9.RData')
        miKeys[[10]] <- c('tcvnienq5hh58bk', 'miRNAtargets_oar_10.RData')
        miKeys[[11]] <- c('qaxs98dvmhj508a', 'miRNAtargets_oar_11.RData')
    }
    if (org =='rno')
    {
        miKeys[[1]]  <- c('w5pk0lucthmijqf', 'miRNAtargets_rno_1.RData')
        miKeys[[2]]  <- c('fsq8c2699yedfcw', 'miRNAtargets_rno_2.RData')
        miKeys[[3]]  <- c('cf2oj6oy6n1u3lf', 'miRNAtargets_rno_3.RData')
        miKeys[[4]]  <- c('86er146w3ejv0ke', 'miRNAtargets_rno_4.RData')
        miKeys[[5]]  <- c('ecwmjmknf090mbb', 'miRNAtargets_rno_5.RData')
        miKeys[[6]]  <- c('ufh9hjk6anatmkj', 'miRNAtargets_rno_6.RData')
        miKeys[[7]]  <- c('2quj6pkjip1asup', 'miRNAtargets_rno_7.RData')
        miKeys[[8]]  <- c('t4h9v5y149xbhro', 'miRNAtargets_rno_8.RData')
        miKeys[[9]]  <- c('bvf0g9qm60y3sad', 'miRNAtargets_rno_9.RData')
        miKeys[[10]] <- c('9qocuom43ihl7um', 'miRNAtargets_rno_10.RData')
        miKeys[[11]] <- c('dmtek3fx529ky23', 'miRNAtargets_rno_11.RData')
    }
    if (!org %in% orgs)
    {
        miKeys[[1]]  <- c('vz774t044y2v4bq', 'miRNAtargets_null.RData')
        miKeys[[2]]  <- c('vz774t044y2v4bq', 'miRNAtargets_null.RData')
        miKeys[[3]]  <- c('vz774t044y2v4bq', 'miRNAtargets_null.RData')
        miKeys[[4]]  <- c('vz774t044y2v4bq', 'miRNAtargets_null.RData')
        miKeys[[5]]  <- c('vz774t044y2v4bq', 'miRNAtargets_null.RData')
        miKeys[[6]]  <- c('vz774t044y2v4bq', 'miRNAtargets_null.RData')
        miKeys[[7]]  <- c('vz774t044y2v4bq', 'miRNAtargets_null.RData')
        miKeys[[8]]  <- c('vz774t044y2v4bq', 'miRNAtargets_null.RData')
        miKeys[[9]]  <- c('vz774t044y2v4bq', 'miRNAtargets_null.RData')
        miKeys[[10]] <- c('vz774t044y2v4bq', 'miRNAtargets_null.RData')
        miKeys[[11]] <- c('vz774t044y2v4bq', 'miRNAtargets_null.RData')
    }

    key <- miKeys[[pn]][1]
    x   <- miKeys[[pn]][2]
    url <- paste0("https://dl.dropboxusercontent.com/s/", key, "/", x)

    return(url)
}

