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

#' @export
setGeneric("venn_3way", function(x, a, b, c, size_range, dist_range, ...){standardGeneric("venn_3way")})

setMethod("venn_3way", "loxcode_experiment", function(x, a, b, c, size_range, dist_range, labels = NA){
  samp_names <- c(a, b, c)
  samp <- as.list(samp_names)
  names(samp) <- samp_names

  if(!is.na(labels)){
    samp_names <- labels
  }

  samp <- lapply(samp, function(s){ return(loxcoder::get_valid(x, s) %>%
                                             filter(size >= size_range[1] &
                                                      size <= size_range[2] &
                                                      dist_orig >= dist_range[1] &
                                                      dist_orig <= dist_range[2]))} )

  pairs <- list("n12" = c(1, 2), "n23" = c(2, 3), "n13" = c(1, 3));

  n_pairs <- lapply(pairs, function(x){
    length(intersect(samp[[x[1]]]$code, samp[[x[2]]]$code))
  })

  n_all <- length(intersect(samp[[1]]$code, intersect(samp[[2]]$code, samp[[3]]$code)))

  v <- draw.triple.venn(area1 = nrow(samp[[1]]),
                   area2 = nrow(samp[[2]]),
                   area3 = nrow(samp[[3]]),
                   n12 = n_pairs$n12,
                   n23 = n_pairs$n23,
                   n13 = n_pairs$n13,
                   n123 = n_all,
                   category = samp_names,
                   fontfamily = rep('sans', 7),
                   cat.fontfamily = rep('sans', 3),
                   fill = c('red', 'blue', 'green'))
  return(grobTree(children = v))
})
