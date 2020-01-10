# Utility functions needed by fastPG, vignettes, or scripts.

#' (Re)download a file with verification
#' 
#' Download a file from a URL and optionally verify the file matches a given
#' sha256 checksum. Will default to the url basename, but a different name can
#' be specified. Will only re-download a file that does not exist, but will
#' test a pre-existing file's checksum. If no checksum is provided when a file
#' pre-exists, a warning is given as the file present could be unrelated. It
#' is a fatal error if the checksum fails.
#'
#' @param url The URL specifying the file to download.
#' @param file The name (path) to download the file to.
#'   Defaults to saving into the current directory a file named the same as
#'   the basename of the url. Using `file= tempfile()` lets the system clean
#'   up the file at the end of the session.
#' @param sha256 The sha256 checksum expected from the downloaded file, as a
#'   character string with only hex digits. Case insensitive. Defaults to
#'   na_character_, which means don't check do checksum validation.
#'
#' @return The filename (invisibly)
#' 
#' Note: regardless of whether a checksum succeeds, fails, or is skipped, the
#' file is left as downloaded.
#'   
#' @examples
#' \dontrun{
#'   checksum <- "ea8fac7c65fb589b0d53560f5251f74f9e9b243478dcb6b3ea79b5e36449c8d9"
#'   toFile <- tempfile()
#'   redownload("http://www.example.com", file= toFile, sha256= checksum)
#' }
#' @export
redownload <- function( url, file= basename(url), sha256= NA_character_ ) {
    checkmate::assertString( url, fixed = "://" )
    checkmate::assertString( file, min.chars = 1 )
    checkmate::assertString( sha256, na.ok = TRUE, min.chars = 64 )
    
    alreadyExists <- file.exists( file )
    validating <- ! is.na( sha256 )

    # Get file    
    if ( ! alreadyExists ) {
        saved <- curl::curl_download( url, file )
    }
    else {
        saved <- file
        if ( ! validating ) {
            warning( paste0(
                "Warning: Existing file assumed to match requested download.\n",
                "\tNot checked as no sha256 was provided to test.\n"
            ))
        }
    }
    
    
    if ( validating ) {
        con <- file( saved )
        full_sha256 <- paste0(as.character(openssl::sha256( con )), collapse= "")
        valid <- tolower( sha256 ) == full_sha256
        if ( ! valid ) {
            if ( ! alreadyExists ) {
                stop( "ERROR - sha256 does not match expected for downloaded file: '", saved, "'." )
            }
            else {
                stop( paste0(
                    "ERROR - File already exists but its sha256 does not match that provided.\n",
                    "\tDelete or specify a different download name and try again.\n"
                ))
            }
        }
    }
    
    invisible( saved )
}


#' Convert a hex string to a raw vector
#'
#' @param x A hex string
#'   Besides delimiters, may only contain case insensitive hex digits, i.e.
#'   0-9, a-f, and A-F. Must have an even number of hex digits, and if
#'   delimiters are present, must split the string into 2 digit pairs.
#'   
#'   Converts `NULL`, `character(0)` and `""` to `NULL`, `raw(0)` and `raw(0)`
#'   respectively.
#' @param sep A delimiter separating hex digit pairs.
#'   By default is NA_character_, meaning no delimiter is present. The
#'   delimiter is ignored if input is NULL, character(0), "", or a single two
#'   digit value. Note that a string ends in delimiters, the last delimiter
#'   is ignored, but it is an error for a string to begin with a delimiter.
#'   This is an artifact of the way R implements strsplit().
#' @return A raw vector equivalent to the hex string pairs.
#' 
#' @examples
#' hexStringToRaw( "01ae03ff" )
#' #> 01 ae 03 ff
#' 
#' hexStringToRaw( "01 ae 03 FF", sep= " " )
#' #> 01 ae 03 ff
#' 
#' hexStringToRaw( "01::ae::03::", sep= "::" )
#' #> 01 ae 03
#' 
#' hexStringToRaw( NULL )
#' #> NULL
#' 
#' hexStringToRaw( "" )
#' #> raw(0)
#' 
#' hexStringToRaw( character(0) )
#' #> raw(0)
#' 
#' @export
hexStringToRaw <- function( x, sep= NA_character_ ) {
    checkmate::assertCharacter(x, null.ok= TRUE, max.len= 1, any.missing= FALSE)
    
    # Special case values
    if (is.null(x)) return( NULL )
    if (length(x) == 0 || nchar(x) == 0 ) return( raw() )
    if (nchar(x) == 2 && ! grepl( "[^0-9a-fA-F]", x )) {
        as.raw(strtoi( x, base= 16 ))
    }
    checkmate::assertString(sep, na.ok= TRUE, min.chars= 1 )
    
    if (! is.na( sep )) {
        vec <- strsplit( x, sep, fixed= TRUE )[[1]]
        notHex <- grepl( "[^0-9a-fA-F]", vec )
        if ( any( notHex )) {
            stop( "Delimited hex strings may only contain hex digits and delimiters.")
        }
        size <- sapply( vec, nchar )
        if (unique(size)[[1]] != 2) {
            stop( "All elements of a delimited hex string must contain two hex digits." )
        }
    }
    else {
        notHex <- grepl( "[^0-9a-fA-F]", x )
        if ( any( notHex )) {
            stop( "Undelimited hex strings may only contain hex digits.")
        }
        if ( nchar(x) %% 2 != 0 ) {
            stop( "Undelimited hex strings must have an even number of hex digits" )
        }
        
        # https://stackoverflow.com/questions/2247045
        vec <- strsplit(x, "")[[1]]
        vec <- paste0( vec[c(TRUE, FALSE)], vec[c(FALSE, TRUE)] )
    }
    
    
    as.raw(strtoi( vec, base= 16 ))
}