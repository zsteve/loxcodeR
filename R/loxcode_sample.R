#' @export
load_origin_distmaps <- function(path){
  load_origin_files_wrapper(c(paste(path, 0, sep = '/'),
                            paste(path, 1, sep = '/'),
                            paste(path, 2, sep = '/'),
                            paste(path, 3, sep = '/'),
                            paste(path, 4, sep = '/')))
}

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

remove_existing <- function(x, n){
  if(!any(is.na(match(n, names(x))))){
    x <- x[, -match(n, names(x))]
  }
  return(x)
}

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
  x@decode@data <- remove_existing(x@decode@data, 'is_valid')
  x@decode@data <- cbind(x@decode@data, data.frame(is_valid = is_valid(x@decode@data$code)))
  return(x)
})

setGeneric("data", function(x){ standardGeneric("data") })

#' Access decoded cassette data
#'
#' @export
setMethod("data", "loxcode_sample", function(x){ return(x@decode@data) })

setGeneric("makeid", function(x){ standardGeneric("makeid") });

#' Append a colum of packed cassette IDs, or -1 if it cannot be packed
#'
#' @export
setMethod("makeid", "loxcode_sample", function(x){
  x@decode@data <- remove_existing(x@decode@data, 'id')
  x@decode@data <- cbind(x@decode@data, data.frame(id = pack(data(x)$code, data(x)$is_valid)))
  return(x)
})

setGeneric("get_origin_dist", function(x) { standardGeneric("get_origin_dist") })

#' Fetch distances from origin
#'
#' @export
setMethod("get_origin_dist", "loxcode_sample", function(x){
  x@decode@data <- remove_existing(x@decode@data, 'dist_orig')
  x@decode@data <- cbind(x@decode@data, data.frame(dist_orig = retrieve_dist_origin(data(x)$id, data(x)$size)))
  return(x)
})

setGeneric("valid", function(x){ standardGeneric("valid") })
#' Valid cassettes only
#'
#' @export
setMethod("valid", "loxcode_sample", function(x){
  v=data(x);
  return(v[v$is_valid == TRUE, ])
})
