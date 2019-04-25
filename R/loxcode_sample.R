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
