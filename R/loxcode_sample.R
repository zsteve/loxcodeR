#' S4 class to contain output of decode()
#' 
#' @slot data contains raw output of decode() as a data.frame
#' @slot read_ids list of FASTQ  
setClass (
  "decode_output",

  # Defining slot type
  representation (
    data = "data.frame",
    read_ids = "list"
  ),

  # Initializing slots
  prototype = list(
    data = data.frame(),
    read_ids = list()
  )
)


#' S4 class to represent loxcode experimental sample data
#' 
#' @slot decode A data.frame to contain raw decode data from loxcoder::decode()
#' @slot meta A data.frame for user-defined sample metadata
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

