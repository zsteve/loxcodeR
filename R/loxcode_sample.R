#' Load distance maps for distance-to-origin
#'
#' @param path path to directory containing maps named 0, 1, 2, 3, 4, corresponding to the size index of
#' 13, 9, 7, 5, 3 element cassettes respectively
#' @export
load_origin_distmaps <- function(path){
  load_origin_files_wrapper(c(paste(path, 0, sep = '/'),
                            paste(path, 1, sep = '/'),
                            paste(path, 2, sep = '/'),
                            paste(path, 3, sep = '/'),
                            paste(path, 4, sep = '/')))
}

#' Load distance maps for pairwise distances
#'
#' @param path path to directory containing the maps 0 (for 13-element distances) and 1 (for 9-element distances)
#' @export
load_pair_distmaps <- function(path){
  load_pair_files_wrapper(c(paste(path, 0, sep = '/'),
                            paste(path, 1, sep = '/')))
}

#' Load probability tables (Markov chain)
#'
#' @param path path to directory containing the Markov chain tables
#' @export
load_prob_files <- function(path){
  size_idx <- 0:4
  rec_dist <- 1:15
  l <- lapply(size_idx, function(s) sapply(rec_dist, function(r) paste(path , paste0('size_',s,'_rec',r), sep = '/')))
  load_prob_files_wrapper(l)
}

#' S4 class to contain output of decode()
#'
#' @slot data contains raw output of decode() as a data.frame
#' @slot read_ids list of FASTQ
#' @export
decode_output <- setClass (
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
loxcode_sample <- setClass (
  "loxcode_sample",

  # Defining slot type
  representation (
    decode = "decode_output",
    name = "character",
    meta = "data.frame",
    files = "vector",
    decode_stats = "list",
    consensus_filtered_data = "vector"
  ),
  # Initializing slots
  prototype = list(
    decode = new("decode_output"),
    name = '',
    meta = data.frame(),
    files = c("", ""),
    decode_stats = list(),
    consensus_filtered_data = c("")
  )
)

setMethod("length", "loxcode_sample", function(x) nrow(x@decode@data))

setMethod("nrow", "loxcode_sample", function(x) length(x))

#' Get loxcode sample name (description)
#'
#' @export
setGeneric("name", function(x){standardGeneric("name")})

setMethod("name", "loxcode_sample", function(x) x@name)

#' Set loxcode sample name (description)
#'
#' @export
setGeneric("name<-", function(x, v){standardGeneric("name<-")})

setMethod("name<-", "loxcode_sample", function(x, v){
  x@name <- v
})

setGeneric("validate", function(x){ standardGeneric("validate") })

#' Add cassette validation column to decoded cassette data
#'
#' @export
setMethod("validate", "loxcode_sample", function(x){
  x@decode@data <- remove_existing(x@decode@data, 'is_valid')
  x@decode@data <- cbind(x@decode@data, data.frame(is_valid = is_valid(x@decode@data$code)))
  return(x)
})

setGeneric("impute", function(x) {standardGeneric("impute")})

#' Impute missing code in 13-element cassettes
#'
#' For 13-element cassettes that are missing a single element, the
#' missing element is imputed to minimise the resulting dist_orig.
#' @export
setMethod("impute", "loxcode_sample", function(x){
  x@decode@data$code <- impute_13(x@decode@data$code, x@decode@data$size)
  return(x)
})

#' Get cassette IDs
#'
#' Appends a column of packed cassette IDs, or -1 if it cannot be packed
#' @export
setGeneric("makeid", function(x){ standardGeneric("makeid") });

setMethod("makeid", "loxcode_sample", function(x){
  x@decode@data <- remove_existing(x@decode@data, 'id')
  x@decode@data <- cbind(x@decode@data, data.frame(id = pack(data(x)$code, data(x)$is_valid)))
  return(x)
})

#' Fetch distances from origin (dist_orig)
#'
#' Appends a column of dist_orig values for valid cassettes
#' @export
setGeneric("get_origin_dist", function(x) { standardGeneric("get_origin_dist") })

setMethod("get_origin_dist", "loxcode_sample", function(x){
  x@decode@data <- remove_existing(x@decode@data, 'dist_orig')
  x@decode@data <- cbind(x@decode@data, data.frame(dist_orig = retrieve_dist_origin(data(x)$id, data(x)$size)))
  return(x)
})

#' Get recombination distance distribution
#'
#' Returns the proportion of valid cassettes for each valid recombination distance (dist_orig)
#' @export
setGeneric("get_rec_prob", function(x, size) {standardGeneric("get_rec_prob")})

setMethod("get_rec_prob", "loxcode_sample", function(x, size){
  r <- data.frame(table(valid(x)$dist_orig))
  names(r) <- c('rec', 'prob')
  r$prob <- r$prob/sum(r$prob)
  r$rec <- as.numeric(r$rec)
  return(r)
})

#' Get ensemble generation probability
#'
#' Retrieve ensemble probabilities as a weighted linear combination of Markov probabilities
#' where distances (dist_orig) are weighted using the sample distribution of recombination distances.
#'
#' Results are appended to readout data as a 'prob' column.
#'
#' @param x loxcode_sample object
#' @export
setGeneric("retrieve_prob_ensemble", function(x) {standardGeneric("retrieve_prob_ensemble")})

setMethod("retrieve_prob_ensemble", "loxcode_sample", function(x){
  sizes <- unique(loxcoder::valid(x)$size)
  x@decode@data$prob <- NA
  for(i in sizes){ # stratify by size
    r <- get_rec_prob(x, i)
    mask <- (x@decode@data$is_valid == T & x@decode@data$size == i)
    probs <- rep(0, sum(mask))
    for(j in 1:nrow(r)){
      print(paste('Weight: ', r$prob[j]))
      probs <- probs + r$prob[j]*loxcoder::retrieve_prob(x@decode@data[mask, ]$id, x@decode@data[mask, ]$size, rep(r$rec[j], sum(mask)))
    }
    x@decode@data[mask, ]$prob <- probs
  }
  return(x)
})

#' Access decoded cassette data, valid cassettes only
#'
#' @param x loxcode_sample object
#' @export
setGeneric("valid", function(x){ standardGeneric("valid") })

setMethod("valid", "loxcode_sample", function(x){
  v=data(x);
  return(v[v$is_valid == TRUE, ])
})

#' Access decoded cassette data
#'
#' @param x loxcode_sample object
#' @export
setGeneric("data", function(x){ standardGeneric("data") })

setMethod("data", "loxcode_sample", function(x){ return(x@decode@data) })

