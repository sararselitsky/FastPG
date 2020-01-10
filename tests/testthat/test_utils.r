# Tests for the functions in utils.R

# Note: using describe() and try() mainly for nested context and
# expressability, not really trying to do behavior driven testing.

describe( "utils.R", {
    describe( "redownload()", {
        
        # redownload() is a state changing function that is implemented using
        # curl_download(). Will mock curl_download() with a spy. That means
        # there is no file to create a checksum for, so have to mock the
        # checksum-creating function openssl::sha256() with a spy also. Also
        # this means there is no file to check the existance of, so have to
        # mock that where it needs to return "TRUE".

        describe( "Minimal default run", {
            # Can't put setup inside the with-mock call.
            checksum <- paste0( rep.int( "feedfade", 8 ), collapse= "" )
            sha256_spy <- mockery::mock( hexStringToRaw( checksum ), cycle= TRUE )
            
            fileName <- "index.html"
            curl_download_spy <- mockery::mock( fileName, cycle= TRUE )
            
            with_mock(
                "openssl::sha256"= sha256_spy,
                "curl::curl_download"= curl_download_spy,

                describe( "Minimal default run (mocked)", {
                    url = paste0( "http://www.example.com/", fileName )
                    
                    got <- NULL
                    it( "Runs without error", {
                        expect_silent( got <<- redownload( url ))
                    })
                    it( "Returned the path to the file it should have downloaded", {
                        expect_equal( got, fileName )
                    })
                    it( "Defaulted the filename to the basename of the url", {
                        # Mocked: the file parameter to curl_download determines
                        # the filename as saved.
                        got <- mockery::mock_args(curl_download_spy)[[length(curl_download_spy)]][[2]]
                        expect_equal( got, fileName )
                        expect_equal( got, basename(url))
                    })
                    it( "Downloaded from the expected url", {
                        # Mocked: assumes curl_download() works as expected.
                        got <- mockery::mock_args(curl_download_spy)[[length(curl_download_spy)]][[1]]
                        expect_equal( got, url )
                    })
                    it( "Did not check sha256 by default", {
                        # Mocked: assumes curl_download() works as expected.
                        mockery::expect_called(sha256_spy, 0)
                    })
                    it( "No download, just warning if file exists and not checking sha256", {
                        file.create( fileName )
                        on.exit( unlink( fileName ))
                        wantWarnRE <- paste0(
                            "^Warning: Existing file assumed to match requested download.\n",
                            "\tNot checked as no sha256 was provided to test.\n$"
                        )
                        beforeCount <- length( curl_download_spy )
                        expect_warning( got <<- redownload( url ),
                                        wantWarnRE, perl= TRUE )
                        afterCount <- length( curl_download_spy )
                        expect_equal( got, fileName )
                        expect_equal( beforeCount, afterCount )
                        unlink( fileName )
                    })
                })
            )
        })
        describe( "Specifying the downloaded file name", {
            checksum <- paste0( rep.int( "feedfade", 8 ), collapse= "" )
            sha256_spy <- mockery::mock( hexStringToRaw( checksum ), cycle= TRUE )
            
            fileName <- tempfile()
            curl_download_spy <- mockery::mock( fileName, cycle= TRUE )

            with_mock(
                "openssl::sha256"= sha256_spy,
                "curl::curl_download"= curl_download_spy,
                
                describe( "Specifying the downloaded file name (mocked)", {
                    url = paste0( "http://www.example.com/", "index.html")
                    
                    it( "Runs without error", {
                        expect_silent( got <<- redownload( url, fileName ))
                    })
                    it( "Returned the path to the file it should have downloaded", {
                        expect_equal( got, fileName )
                    })
                    it( "Downloaded to the file as specified, not based on the url.", {
                        # Mocked: the file parameter to curl_download determines
                        # the filename as saved.
                        got <- mockery::mock_args(curl_download_spy)[[length(curl_download_spy)]][[2]]
                        expect_equal( got, fileName )
                        expect_false( got == basename(url))
                    })
                    it( "Downloaded from the expected url", {
                        # Mocked: assumes curl_download() works as expected.
                        got <- mockery::mock_args(curl_download_spy)[[length(curl_download_spy)]][[1]]
                        expect_equal( got, url )
                    })
                    it( "Did not check sha256 by default", {
                        # Mocked: assumes curl_download() works as expected.
                        mockery::expect_called(sha256_spy, 0)
                    })
                    it( "No download, just warning if file exists and not checking sha256", {
                        file.create( fileName )
                        on.exit( unlink( fileName ))
                        wantWarnRE <- paste0(
                            "^Warning: Existing file assumed to match requested download.\n",
                            "\tNot checked as no sha256 was provided to test.\n$"
                        )
                        beforeCount <- length( curl_download_spy )
                        expect_warning( got <<- redownload( url, fileName ),
                                        wantWarnRE, perl= TRUE )
                        afterCount <- length( curl_download_spy )
                        expect_equal( got, fileName )
                        expect_equal( beforeCount, afterCount )
                        unlink( fileName )
                    })
                    
                })
            )
        })
        describe( "Specifying the checksum and filename", {
            checksum <- paste0( rep.int( "feedfade", 8 ), collapse= "" )
            sha256_spy <-  mockery::mock( hexStringToRaw( checksum ), cycle= TRUE )
            
            fileName <- tempfile()
            curl_download_spy <- mockery::mock( fileName, cycle= TRUE )
            
            file_exists_spy <- mockery::mock( TRUE, cycle=TRUE )
            
            url = paste0( "http://www.example.com/", "index.html")
            
            with_mock(
                "openssl::sha256"= sha256_spy,
                "curl::curl_download"= curl_download_spy,
                
                describe( "When the file needs to be downloaded", {
                    
                    it( "Runs without error", {
                        expect_silent( got <<- redownload( url, fileName, checksum ))
                    })
                    it( "Returned the path to the file it should have downloaded", {
                        expect_equal( got, fileName )
                    })
                    it( "Downloaded to the file as specified, not based on the url.", {
                        # Mocked: the file parameter to curl_download determines
                        # the filename as saved.
                        got <- mockery::mock_args(curl_download_spy)[[length(curl_download_spy)]][[2]]
                        expect_equal( got, fileName )
                        expect_false( got == basename(url))
                    })
                    it( "Downloaded from the expected url", {
                        # Mocked: assumes curl_download() works as expected.
                        got <- mockery::mock_args(curl_download_spy)[[length(curl_download_spy)]][[1]]
                        expect_equal( got, url )
                    })
                    it( "Calculates the sha256 of the file expected", {
                        redownload( url, fileName, checksum )
                        
                        got <- mockery::mock_args(sha256_spy)[[length(sha256_spy)]][[1]]
                        expect_equal( summary(got)$description, fileName )
                    })
                    it( "Errors if checksum does not match when file does not exist", {
                        checksum <- paste0( rep.int( "baadbaad", 8 ), collapse= "" )
                        wantErrorRE= paste0(
                            "ERROR - sha256 does not match expected for downloaded file: '",
                            fileName, "'."
                        )
                        expect_error( redownload( url, fileName, checksum ),
                                      wantErrorRE )
                    })
                })
,               with_mock(
                    "file.exists" = file_exists_spy,
                    describe( "When file already exists", {
                            it( "Does no downloads", {
                                beforeCount <- length( curl_download_spy )
                                expect_silent( got <<- redownload( url, fileName, checksum ))
                                afterCount <- length( curl_download_spy )
                                expect_equal( got, fileName )
                                expect_equal( beforeCount, afterCount )
                            })
                            it( "Errors if checksum does not match", {
                                checksum <- paste0( rep.int( "baadbaad", 8 ), collapse= "" )
                                wantErrorRE= paste0(
                                    "^ERROR - File already exists but its sha256 does not match that provided.\n",
                                    "\tDelete or specify a different download name and try again.\n$"
                                )
                                expect_error( redownload( url, fileName, checksum ),
                                              wantErrorRE )
                            })
                        })
                )
            )
        })
        describe( "Specifying just the checksum", {
            checksum <- paste0( rep.int( "feedfade", 8 ), collapse= "" )
            sha256_spy <-  mockery::mock( hexStringToRaw( checksum ), cycle= TRUE )
            
            fileName <- "index.html"
            curl_download_spy <- mockery::mock( fileName, cycle= TRUE )
            
            file_exists_spy <- mockery::mock( TRUE, cycle=TRUE )
            
            url = paste0( "http://www.example.com/", "index.html")
            
            with_mock(
                "openssl::sha256"= sha256_spy,
                "curl::curl_download"= curl_download_spy,
                
                describe( "If the file needs to be downloaded", {
                    it( "Runs without error", {
                        expect_silent( got <<- redownload( url, sha256= checksum ))
                    })
                    it( "Returned the path to the file it should have downloaded", {
                        expect_equal( got, fileName )
                    })
                    it( "Downloaded to the file based on the url.", {
                        # Mocked: the file parameter to curl_download determines
                        # the filename as saved.
                        got <- mockery::mock_args(curl_download_spy)[[length(curl_download_spy)]][[2]]
                        expect_equal( got, fileName )
                        expect_equal( got, basename(url))
                    })
                    it( "Downloaded from the expected url", {
                        # Mocked: assumes curl_download() works as expected.
                        got <- mockery::mock_args(curl_download_spy)[[length(curl_download_spy)]][[1]]
                        expect_equal( got, url )
                    })
                    it( "Calculates the sha256 of the file expected", {
                        # Mocked: assumes curl_download() works as expected.
                        redownload( url, sha256= checksum )
                        
                        got <- mockery::mock_args(sha256_spy)[[length(sha256_spy)]][[1]]
                        expect_equal( summary(got)$description, fileName )
                    })
                    it( "Errors if checksum does not match when file does not exist", {
                        checksum <- paste0( rep.int( "baadbaad", 8 ), collapse= "" )
                        wantErrorRE= paste0(
                            "ERROR - sha256 does not match expected for downloaded file: '",
                            fileName, "'."
                        )
                        expect_error( redownload( url, sha256= checksum ),
                                      wantErrorRE )
                    })
                }),
                with_mock(
                    "file.exists" = file_exists_spy,
                    describe( "When file already exists", {
                        it( "Does not download file", {
                            beforeCount <- length( curl_download_spy )
                            expect_silent( got <<- redownload( url, sha256= checksum ))
                            afterCount <- length( curl_download_spy )
                            expect_equal( got, fileName )
                            expect_equal( beforeCount, afterCount )
                        })
                        it( "Errors if checksum does not match.", {
                            checksum <- paste0( rep.int( "baadbaad", 8 ), collapse= "" )
                            wantErrorRE= paste0(
                                "^ERROR - File already exists but its sha256 does not match that provided.\n",
                                "\tDelete or specify a different download name and try again.\n$"
                            )
                            expect_error( redownload( url, sha256= checksum ),
                                          wantErrorRE )
                        })
                    })
                )
            )
        })
        describe( "Parameter checking", {
            checksum <- paste0( rep.int( "feedfade", 8 ), collapse= "" )
            sha256_spy <-  mockery::mock( hexStringToRaw( checksum ), cycle= TRUE )
            
            fileName <- "index.html"
            curl_download_spy <- mockery::mock( fileName, cycle= TRUE )
            
            url = paste0( "http://www.example.com/", "index.html")
            
            with_mock(
                "openssl::sha256"= sha256_spy,
                "curl::curl_download"= curl_download_spy,
                
                describe( "Parameter checking (mocked)", {
                    it( "tests the url= parameter", {
                        badUrls <- list( null=NULL, missing=NA_character_,
                                    empty="", badformat= "http:/example.com/oops",
                                    vector=c("https://example.com", "https://example.com"),
                                    list=list("x"), func= sum )
                        for ( name in names( badUrls )) {
                            wantErrorRE <- "^Assertion on 'url' failed: "
                            expect_error( redownload( url= badUrls[[!!name]] ),
                                          wantErrorRE )
                        }
                    })
                    it( "tests the file= parameter", {
                        url = "http://www.example.com/"
                        badFiles <- list( null=NULL, missing=NA_character_,
                                       empty="",
                                       vector=c("https://example.com", "https://example.com"),
                                       list=list("x"), func= sum )
                        for ( name in names( badFiles )) {
                            wantErrorRE <- "^Assertion on 'file' failed: "
                            expect_error( redownload( url, file= badFiles[[!!name]] ),
                                          wantErrorRE )
                        }
                    })
                    it( "tests the checksum= parameter", {
                        url = "http://www.example.com/"
                        badChecksum <- list(
                            null= NULL, empty= "", badFormat= "123",
                            vector= c( checksum, checksum ),
                            list= list( x= checksum ), func= sum
                        )
                        for ( name in names( badChecksum )) {
                            wantErrorRE <- "^Assertion on 'sha256' failed: "
                            expect_error( redownload( url, sha256= badChecksum[[!!name]] ),
                                           wantErrorRE )
                            expect_error( redownload( url, fileName, sha256= badChecksum[[!!name]] ),
                                          wantErrorRE )
                        }
                    })
                    it( "parameter order is as expected", {
                        expect_silent( got <<- redownload( url, fileName, checksum ))
                        expect_equal( got, fileName )
                    })
                    it ( "parameter names are as expected", {
                        expect_silent( got <<- redownload( file= fileName, sha256= checksum, url= url ))
                        expect_equal( got, fileName )
                    })
                })
            )
        })
        
    })
    describe( "hexStringToRaw()", {
        describe( "Converting non-delimited hex strings", {
            it( "Converts hex strings to raw vectors case insensitively", {
                expect_equal( hexStringToRaw( "0001afff" ),
                              as.raw( c( 0, 1, 175, 255 )))
                expect_equal( hexStringToRaw( "0001AFff" ),
                              as.raw( c( 0, 1, 175, 255 )))
            })
            it( "Converts NULL, length 0, \"\" and single byte inputs to raw vectors", {
                expect_null( hexStringToRaw( NULL ))
                expect_equal( hexStringToRaw( character( 0 )), raw( 0 ))
                expect_equal( hexStringToRaw( "" ), raw( 0 ))
                expect_equal( hexStringToRaw( "a0" ), as.raw( 160 ))
            })
            it( "Errors if character count in input is not even", {
                wantErrorRE <- "Undelimited hex strings must have an even number of hex digits"
                expect_error( hexStringToRaw("abc"), wantErrorRE )
                expect_error( hexStringToRaw("a"), wantErrorRE )
            })
            it( "Errors if non-hex, non-delimiter character in input", {
                wantErrorRE <- "Undelimited hex strings may only contain hex digits."
                expect_error( hexStringToRaw("gg"), wantErrorRE )
                expect_error( hexStringToRaw("11 22"), wantErrorRE )
            })
        })
        describe( "Using non-default delimiters", {
            it( "Converts delimited hex strings to raw vectors, case insensitively", {
                expect_equal( hexStringToRaw( "00 01 af ff", sep= " " ),
                              as.raw( c( 0, 1, 175, 255 )))
                expect_equal( hexStringToRaw( "00::01::AF::ff", sep='::' ),
                              as.raw( c( 0, 1, 175, 255 )))
            })
            it( "Converts NULL, length 0, \"\" and single byte inputs to raw vectors", {
                expect_null( hexStringToRaw( NULL, sep=" " ))
                expect_equal( hexStringToRaw( character( 0 ), sep= " " ), raw( 0 ))
                expect_equal( hexStringToRaw( "" ), sep= "::", raw( 0 ))
                expect_equal( hexStringToRaw( "a0", '-' ), as.raw( 160 ))
                
                # Special corner case due to strsplit!
                expect_equal( hexStringToRaw( "a0-0a-", '-' ), as.raw( c( 160, 10 )))
            })
            it( "Errors if character count in input is not even", {
                wantErrorRE <- "All elements of a delimited hex string must contain two hex digits."
                expect_error( hexStringToRaw( "0 1 af ff", sep= " " ), wantErrorRE )
                expect_error( hexStringToRaw( "0000 1111 2222", sep= " " ), wantErrorRE )
                expect_error( hexStringToRaw( "::", sep= "::" ), wantErrorRE )
                expect_error( hexStringToRaw( "::01", sep= "::" ), wantErrorRE )
            })
            it( "Errors if non-hex, non-delimiter character in input", {
                wantErrorRE <- "Delimited hex strings may only contain hex digits and delimiters."
                expect_error( hexStringToRaw("00 01 af ff", sep="::"), wantErrorRE )
                expect_error( hexStringToRaw( "00 01::af ff", sep= " " ), wantErrorRE )
                expect_error( hexStringToRaw( "0 0g::af f", sep= " " ), wantErrorRE )
            })
        })
        describe( "Parameter checking", {
            it( "tests the x= data parameter", {
                badX <- list( vector=c("00", "01"),
                                 list=list("0001"), func= sum )
                for ( name in names( badX )) {
                    wantErrorRE <- "^Assertion on 'x' failed: "
                    expect_error( hexStringToRaw( x= badX[[!!name]] ),
                                  wantErrorRE )
                }
            })
            it( "tests the sep= parameter", {
                badSep <- list( null=NULL, empty="",
                                vector=c(" ", "::"),
                                list=list(" "), func= sum )
                for ( name in names( badSep )) {
                    wantErrorRE <- "^Assertion on 'sep' failed: "
                    expect_error( hexStringToRaw( "01 01", sep= badSep[[!!name]] ),
                                  wantErrorRE )
                }
            })
            it( "has the correct parameter order", {
                got <- hexStringToRaw( "00::00", '::' )
                want <- as.raw( c( 0, 0 ))
                expect_equal( got, want )
            })
            it( "has the correct parameter names", {
                got <- hexStringToRaw( sep= "::", x= "00::00" )
                want <- as.raw( c( 0, 0 ))
                expect_equal( got, want )
                
            })
        })
    })
})