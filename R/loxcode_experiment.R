#' Loxcode experiment object
#'
#' The `loxcode_experiment` object enables handling of data from multiple samples.
#'
#' @slot name string, name of the experiment
#' @slot suffix_R1 string, R1 suffix that is appended to the filename prefix, e.g. `_R1_001.fastq`
#' @slot suffix_R2 string, R2 suffix that is appended to the filename prefix, e.g. `_R2_001.fastq`
#' @slot dir string, directory containing R1, R2 *.fastq files. Must end with '/'
#' @slot samples list, contains loxcode_sample objects that can be accessed by the sample name
#' @slot samp_table data.frame, user-specified table that can be loaded from an Excel spreadsheet
loxcode_experiment <- setClass(
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
#' Load data from FASTQ for all samples in a loxcode_experiment object
#'
#' For each sample in the loxcode_experiment, barcodes readout is performed from FASTQ
#' followed by element imputation for the 13-element codes. Loxcodes are then validated, IDs fetched,
#' and distance-from-origin is retrieved.
#'
#' @param x loxcode_experiment object for which to load samples
#' @param full boolean, whether to produce full debugging output (default is FALSE, this uses significantly more memory)
#' @return loxc# ode_experiment object with sample data loaded
#' @export
setGeneric("load_samples", function(x, ...) {standardGeneric("load_samples")})

setMethod("load_samples", "loxcode_experiment", function(x, full = F){
    x@samples <- lapply(names(x@samples), function(z){
    print(z)
    samp_table_sliced <- x@samp_table[match(z, x@samp_table$sample), ]
    out <- loxcoder::decode(c(paste0(x@dir, x@samples[[z]], x@suffix_R1), paste0(x@dir, x@samples[[z]], x@suffix_R2)),
                       name = z, meta = samp_table_sliced$meta,
                     min_r1_len = samp_table_sliced$min_r1_len,
                     min_r2_len = samp_table_sliced$min_r2_len, full = full)
    out <- loxcoder::impute(out)
    out <- loxcoder::validate(out)
    out <- loxcoder::makeid(out)
    out <- loxcoder::get_origin_dist(out)
    return(out)
  })
  names(x@samples) <- x@samp_table$sample
  return(x)
})

#' Fetch names of samples
#'
#' The names returned by this function can be used to refer to each individual sample
#' @param x loxcode_experiment object
#' @return names (handles) of samples
#' @export
setGeneric("sampnames", function(x){standardGeneric("sampnames")})

setMethod("sampnames", "loxcode_experiment", function(x){
  return(names(x@samples))
})

#' Set names of samples
#'
#' @export
setGeneric("sampnames<-", function(x, v){ standardGeneric("sampnames<-")})

setMethod("sampnames<-", "loxcode_experiment", function(x, v){
  if(length(v) != length(unique(v))){
    stop("Sample names are not unique")
  }
  names(x@samples) <- v
})

#' Fetch experiment name
#'
#' @param x loxcode_experiment object
#' @return name of experiment
#' @export
setGeneric("name", function(x){ standardGeneric("name")})

setMethod("name", "loxcode_experiment", function(x){
  return(x@name)
})

#' Set experiment name
#'
#' @export
setGeneric("name<-", function(x, v) {standardGeneric("name<-")})

setMethod("name<-", "loxcode_experiment", function(x, v){
  x@name <- v
})

#' Get loxcode_sample object by sample name
#'
#' @param x loxcode_experiment object
#' @param s sample name
#' @return loxcode_sample object corresponding to s
#' @export
setGeneric("sample", function(x, s){standardGeneric("sample")})

setMethod("sample", "loxcode_experiment", function(x, s){
  return(x@samples[[s]])
})

#' Get unvalidated readout data for a given sample
#'
#' Equivalent to performing loxcoder::data() on the sample s
#' @param x loxcode_experiment object
#' @param s sample name
#' @return data.frame containing unvalidated readout data for sample s
#' @export
setGeneric("get_data", function(x, s){standardGeneric("get_data")})

setMethod("get_data", "loxcode_experiment", function(x, s){
  return(loxcoder::data(x@samples[[s]]))
})

#' Get validated readout data for a given sample
#'
#' Equivalent to performing loxcoder::valid() on the sample s
#' @param x loxcode_experiment object
#' @param s sample name
#' @export
setGeneric("get_valid", function(x, s){standardGeneric("get_valid")})

setMethod("get_valid", "loxcode_experiment",function(x, s){
  return(loxcoder::valid(x@samples[[s]]))
})

#' Get sample table
#'
#' @export
setGeneric("samptable", function(x){ standardGeneric("samptable")})

setMethod("samptable", "loxcode_experiment", function(x){
  return(x@samp_table)
})

#' Set sample table
#'
#' @export
setGeneric("samptable<-", function(x, value){ standardGeneric("samptable<-")})

setMethod("samptable<-", "loxcode_experiment", function(x, value){
  x@samp_table <- value
  return(x)
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
#' @param full boolean, whether to produce full debugging output (default is FALSE, this uses significantly more memory)
#' @return loxcode_experiment object
#' @export
load_from_xlsx <- function(name, s, dir, suffix_R1, suffix_R2, load = TRUE, full = FALSE){
  x <- new("loxcode_experiment", name = name, dir = dir, suffix_R1 = suffix_R1, suffix_R2 = suffix_R2)
  x@samp_table <- data.frame(read_excel(s))
  x@samples = lapply(x@samp_table$prefix, function(z) z)
  names(x@samples) <- x@samp_table$sample
  if(load){
    x <- loxcoder::load_samples(x, full = full)
  }
  return(x)
}

#' Merge two loxcode_sample objects
#'
#' Barcodes from s1 and s2 are merged. For barcodes present in both samples, read counts are
#' summed. Sample descriptors and other members are concatenated in order.
#' @param s1 sample 1
#' @param s2 sample 2
#' @return loxcode_sample object corresponding to merged sample
#' @export
merge_sample <- function(s1, s2){
  m <- merge(s1@decode@data, s2@decode@data, by = c('code', 'size', 'is_valid', 'id', 'dist_orig'), all = T)
  m$count.x[is.na(m$count.x)] <- 0
  m$count.y[is.na(m$count.y)] <- 0
  m$count <- m$count.x + m$count.y
  m <- m[, !(names(m) %in% c('count.x', 'count.y'))]
  m <- m[order(m$code), ]
  m <- m[, c(ncol(m), 1:(ncol(m)-1))]
  s <- new('loxcode_sample')
  s@decode <- new('decode_output')
  s@decode@data <- m
  # for the rest of the members, we just duplicate
  s@decode@saturation <- list(s1@decode@saturation, s2@decode@saturation)
  s@decode@read_ids <- list(s1@decode@read_ids, s2@decode@read_ids)
  s@meta <- cbind(s1@meta, s2@meta)
  s@name <- paste(s1@name, s2@name, sep = '_')
  s@files <- cbind(s1@files, s2@files)
  s@consensus_filtered_data <- list(s1@consensus_filtered_data,
                                    s2@consensus_filtered_data)
  s@decode_stats <- list(s1@decode_stats,
                         s2@decode_stats)
  return(s)
}

#' Merge samples by label
#'
#' Merge samples present in a loxcode_experiment object according to a specified label column in
#' samp_table. All samples with a common value for the label column are merged together. The resulting
#' loxcode_experiment object contains loxcode_samples named by the corresponding labels in by.
#'
#' @param x loxcode_experiment object
#' @param by name of column in samp_table to merge by
#' @return loxcode_experiment object containing merged samples
#' @export
setGeneric('merge_by', function(x, by){ standardGeneric('merge_by') })

setMethod('merge_by', 'loxcode_experiment', function(x, by){
  vals <- unique(x@samp_table[, by])
  x_merged <- new("loxcode_experiment", name = paste(x@name, 'merge_by', by, sep = '_'),
                  dir = x@dir,
                  suffix_R1 = x@suffix_R1,
                  suffix_R2 = x@suffix_R2)
  x_merged@samp_table <- x@samp_table
  x_merged@samples <- lapply(vals, function(z){
    print(z)
    index <- which(x@samp_table[, by] == z)
    if(length(index) > 0){
      out <- x@samples[[index[1]]]
      for(i in 2:length(index)){
        out <- merge_sample(out, x@samples[[index[i]]])
      }
      return(out)
    }else{
      return(NA)
    }
  })
  names(x_merged@samples) <- vals
  return(x_merged)
})

#' 2-way Venn diagram
#'
#' Plot 2-way Venn diagram showing common and distinct barcodes between two given samples in an experiment using VennDiagram package
#' @param x loxcode_experiment object
#' @param a sample A name
#' @param b sample B name
#' @param size_range range of loxcode sizes to show in Venn diagram. Both lower and upper bounds are inclusive
#' @param dist_range range of dist_orig values to show in Venn diagram. Both lower and upper bounds are inclusive
#' @param labels optional list of custom labels for samples A, B, otherwise sample names for A, B are used.
#' @return grobTree containing the 2-way Venn diagram graphic
#' @export
setGeneric("venn_2way", function(x, a, b, size_range, dist_range, ...){standardGeneric("venn_2way")})

setMethod("venn_2way", "loxcode_experiment", function(x, a, b, size_range, dist_range, labels = NA){
  samp_names <- c(a, b)
  samp <- as.list(samp_names)

  if(!is.na(labels)){
    samp_names <- labels
  }

  samp <- lapply(samp, function(s){
    return(loxcoder::get_valid(x, s) %>%
             filter(size >= size_range[1] &
                      size <= size_range[2] &
                      dist_orig >= dist_range[1] &
                      dist_orig <= dist_range[2]))
  })

  v <- draw.pairwise.venn(area1 = nrow(samp[[1]]),
                          area2 = nrow(samp[[2]]),
                          cross.area = length(intersect(samp[[1]]$code, samp[[2]]$code)),
                          category = samp_names,
                          fontfamily = rep('sans', 3),
                          cat.fontfamily = rep('sans', 2),
                          fill = c('red', 'blue'))
  return(grobTree(children = v))
})

#' 3-way Venn diagram
#'
#' Plot 3-way Venn diagram showing common and distinct barcodes between two given samples in an experiment using VennDiagram package
#' @param x loxcode_experiment object
#' @param a sample A name
#' @param b sample B name
#' @param c sample C name
#' @param size_range range of loxcode sizes to show in Venn diagram. Both lower and upper bounds are inclusive
#' @param dist_range range of dist_orig values to show in Venn diagram. Both lower and upper bounds are inclusive
#' @param labels optional list of custom labels for samples A, B, C, otherwise sample names for A, B, C are used.
#' @return grobTree containing the 3-way Venn diagram graphic
#' @export
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
