describe( "makeExtData.Rscript", {
    
    script <- file.path( "data-raw", "makeExtData.Rscript" )
    wantFile <- file.path( "inst", "extdata", "Levine_32dim_10k.fcs" )

    oldDir <- getwd()
    on.exit( setwd( oldDir ), add = TRUE, after=FALSE )
    setwd( file.path( "..", ".." ))
    
    if ( file.exists( wantFile )) {
        unlink( wantFile )
    }
    
    describe( "Test setup", {
        it( "File does not exist before running.", {
            expect_false( file.exists( wantFile ))
        })
        it( "Have internet", {
            checksum= "ea8fac7c65fb589b0d53560f5251f74f9e9b243478dcb6b3ea79b5e36449c8d9"
            expect_silent( redownload("http://www.example.com", file= tempfile(), sha256= checksum))
        })
    })
    
    describe( "nominal run", {
        
        # 41.5 meg; 0.2 M/sec download; 1 min extra
        expectDownloadSeconds <- as.integer(( 41.5 / 0.2 ) + 60 )
        stdOutFile <- tempfile()
        stdErrFile <- tempfile()
        retVal <- system2(
            "/usr/bin/env",
            args= c("Rscript", script),
            stdout= stdOutFile,
            stderr= stdErrFile,
            wait= TRUE,
            timeout= expectDownloadSeconds
        )
        
        it( "Returns success (0) when run", {
            expect_equal(retVal, 0)
        })
        it( "Creates the expected file", {
            expect_true( file.exists( wantFile ))
        })
        it( "Is a valid FCS file", {
            expect_true(flowCore::isFCSfile( wantFile ))
        })
        it( "Has the expected subsampled data", {
            dat <- flowCore::read.FCS( wantFile )
            expect_equal( nrow(dat), 10L * 1000L )
            expect_equal( ncol(dat), 32L )
        })
        it( "No stderr generated while running", {
            expect_equal( readLines( stdErrFile, warn= FALSE ), character(0) )
        })
        it( "outputs file create to stdout.", {
            expect_equal( readLines( stdOutFile, warn= FALSE ), wantFile )
        })
    })
})