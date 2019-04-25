#' S4 class to represent loxcode experimental sample data
#' 
#' @slot decode A data.frame to contain raw decode data from loxcoder::decode()
#' @slot meta A data.frame for user-defined sample metadata
setClass (
  "loxcode_sample",

  # Defining slot type
  representation (
    decode = "data.frame",
    meta = "data.frame"
  ),

  # Initializing slots
  prototype = list(
    decode = data.frame(), 
    meta = data.frame()
  )
)
