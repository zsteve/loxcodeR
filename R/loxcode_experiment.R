#' S4 class for a single loxcode experiment, enabling handling of multiple samples
#'
#' @slot name string, name of the experiment
#' @slot suffix_R1 string, R1 suffix
#' @slot suffix_R2 string, R2 suffix
#' @slot dir string, directory containing R1, R2 fastq files
#' @slot samples list, contains loxcode_sample objects
#' @slot samp_table data.frame, user-specified table.

setClass(
  "loxcode_experiment",

  representation(
    name = "character",
    suffix_R1 = "character",
    suffix_R2 = "character",
    dir = "character",
    samples = "list",
    samp_table = "data.frame"
  ),

  prototype = list(
    name = '',
    suffix_R1 = '',
    suffix_R2 = '',
    dir = '',
    samples = list(),
    samp_table = data.frame()
  )
)

setGeneric("load_samples", function(x) {standardGeneric("load_samples")})
setMethod("load_samples", "loxcode_experiment", function(x){
    x@samples <- lapply(names(x@samples), function(z){
    print(z)
    samp_table_sliced <- x@samp_table[match(z, x@samp_table$sample), ]
    out <- loxcoder::decode(c(paste0(x@dir, x@samples[[z]], x@suffix_R1), paste0(x@dir, x@samples[[z]], x@suffix_R2)),
                       name = z, meta = samp_table_sliced$meta,
                     min_r1_len = samp_table_sliced$min_r1_len,
                     min_r2_len = samp_table_sliced$min_r2_len)
    out <- loxcoder::impute(out)
    out <- loxcoder::validate(out)
    out <- loxcoder::makeid(out)
    out <- loxcoder::get_origin_dist(out)
    return(out)
  })
  names(x@samples) <- x@samp_table$sample
  return(x)
})

#' Get loxcode_sample object
#'
#' @param x loxcode_experiment object
#' @param s string, sample name
#' @export
setGeneric("sample", function(x, s){standardGeneric("sample")})

setMethod("sample", "loxcode_experiment", function(x, s){
  return(x@samples[[s]])
})

#' @export
setGeneric("get_data", function(x, s){standardGeneric("get_data")})

setMethod("get_data", "loxcode_experiment", function(x, s){
  return(loxcoder::data(x@samples[[s]]))
})

#' @export
setGeneric("get_valid", function(x, s){standardGeneric("get_valid")})

setMethod("get_valid", "loxcode_experiment",function(x, s){
  return(loxcoder::valid(x@samples[[s]]))
})

#' Create experiment
#'
#' Create a loxcode_experiment from xlsx sample table.
#'
#' Must contain the following columns:
#'
#' * `sample` -- sample name. This is used as the unique identifier for the sample
#'
#' * `prefix` -- sample prefix. suffix_R1 and suffix_R2 are appended to this
#'
#' * `meta` -- sample metadata. This is user specified
#'
#' * `min_r1_len` -- minimum R1 length
#'
#' * `min_r2_len` -- minimum R2 length
#'
#' @param name string, name of the experiment
#' @param s string, path to user-specified Excel file containing experiment samples.
#' @param dir string, path to directory containing fastq files
#' @param suffix_R1 string, R1 suffix
#' @param suffix_R2 string, R2 suffix
#' @param load boolean, whether to load samples or not (default is TRUE)
#'
#' @export
load_from_xlsx <- function(name, s, dir, suffix_R1, suffix_R2, load = TRUE){
  x <- new("loxcode_experiment", name = name, dir = dir, suffix_R1 = suffix_R1, suffix_R2 = suffix_R2)
  x@samp_table <- data.frame(read_excel(s))
  x@samples = lapply(x@samp_table$prefix, function(z) z)
  names(x@samples) <- x@samp_table$sample
  if(load){
    x <- loxcoder::load_samples(x)
  }
  return(x)
}


