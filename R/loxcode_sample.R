#' S4 class to contain output of decode()
#' 
#' @slot data contains raw output of decode() as a data.frame
#' @slot read_ids list of FASTQ 
#' @export
setClass (
  "decode_output",

  # Defining slot type
  representation (
    data = "data.frame",
		saturation = "vector",
    read_ids = "list"
  ),

  # Initializing slots
  prototype = list(
		saturation = c(),
    data = data.frame(),
    read_ids = list()
  )
)


#' S4 class to represent loxcode experimental sample data
#' 
#' @slot decode A data.frame to contain raw decode data from loxcoder::decode()
#' @slot meta A data.frame for user-defined sample metadata
#' @export
setClass (
  "loxcode_sample",

  # Defining slot type
  representation (
    decode = "decode_output",
    meta = "data.frame",
    files = "vector" 
  ),

  # Initializing slots
  prototype = list(
    decode = new("decode_output"), 
    meta = data.frame(),
    files = c("", "")
  )  
)

setMethod("length", "loxcode_sample", function(x) nrow(x@decode@data))
setMethod("nrow", "loxcode_sample", function(x) length(x))

setGeneric("validate", function(x){ standardGeneric("validate") })

#' Add cassette validation column to decoded cassette data
#'
#' @export
setMethod("validate", "loxcode_sample", function(x){
  x@decode@data <- cbind(x@decode@data, data.frame(is_valid = is_valid(x@decode@data$code)))
  return(x)
})

setGeneric("get_data", function(x){ standardGeneric("get_data") })

#' Access decoded cassette data
#'
#' @export
setMethod("get_data", "loxcode_sample", function(x){ return(x@decode@data) })


