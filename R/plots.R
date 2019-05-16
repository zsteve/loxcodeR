setGeneric("size_plot", function(x) {standardGeneric("size_plot")})

#' Plots distribution of sizes found in samples
#'
#' @export
setMethod("size_plot", "loxcode_sample", function(x){
  g <- ggplot(data = data(x)) + geom_bar(aes(as.factor(size), fill = is_valid), position = "stack", stat = 'count') +
    xlab('size') + ylab('diversity') +
    ggtitle(loxcoder::name(x))
  return(g)
})

setGeneric("dist_orig_plot", function(x, size) {standardGeneric("dist_orig_plot")})
#' @export
setMethod("dist_orig_plot", "loxcode_sample", function(x, size){
  u <- valid(x)
  u <- u[u$size == size, ]
  g <- ggplot(data = u) + geom_bar(aes(dist_orig)) + ggtitle(sprintf("size = %d", size))+ scale_x_continuous(breaks = 0:10, limits = c(0, 10)) +
    xlab("Distance from origin") + ylab("Diversity") + ggtitle(loxcoder::name(x))
  return(g)
})

setGeneric("rank_count_plot", function(x, size) {standardGeneric("rank_count_plot")})
#' @export
setMethod("rank_count_plot", "loxcode_sample", function(x, size){
  u <- valid(x)
  u <- u[u$size == size, ]
  u <- u[order(u$count, decreasing = T), ]
  # u$count <- u$count/sum(u$count) # use frequency instead of raw count
  ggplot(data = u) + geom_point(aes(x = 1:nrow(u), y = log10(count))) +
    scale_x_log10(breaks = 1:10, labels = u$code[1:10]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    xlab('code') + ylab('log10(count)') + ggtitle(sprintf('%s; size = %d', loxcoder::name(x), size)) +
    geom_point(aes(x = 1:nrow(u), y = dist_orig*max(log10(count))/max(dist_orig)), color = 'red', alpha = 0.1) +
    scale_y_continuous(sec.axis = sec_axis(~.*1/(max(log10(u$count))/max(u$dist_orig)), name = 'distance'))
})

setGeneric("pair_comparison_plot", function(x1, x2) {standardGeneric("pair_comparison_plot")})
#' @export
setMethod("pair_comparison_plot", "loxcode_sample", function(x1, x2){
  u <- get_comparison_table(x1, x2)
  g <- ggplot(data = u) + geom_point(aes(y = 1 + rep1_count , x = 1 + rep2_count, color = as.factor(size))) +
    scale_x_log10() + scale_y_log10() +
    geom_abline(a = 1, b = 0) + ylab(loxcoder::name(x1)) + xlab(loxcoder::name(x2)) +
    ggtitle(paste(loxcoder::name(x1), 'vs', loxcoder::name(x2)))
  return(g)
})

barcode_union <- function(rep1, rep2){
  return(unique(c(loxcoder::valid(rep1)$code, loxcoder::valid(rep2)$code)))
}

get_barcode_stats_rep <- function(union_bc, rep){
  index <- match(union_bc, loxcoder::valid(rep)$code)
  u <- loxcoder::valid(rep)[index, ]
  u$count[is.na(u$count)] <- 0
  u$code <- union_bc
  return(u)
}

get_comparison_table <- function(rep1, rep2){
  bc_union <- barcode_union(rep1, rep2)
  u1 <- get_barcode_stats_rep(bc_union, rep1)
  u2 <- get_barcode_stats_rep(bc_union, rep2)
  u1$size <- ifelse(is.na(u1$size), u2$size, u1$size)
  u2$size <- u1$size
  # scale by total number of reads
  u1$count <- u1$count
  u2$count <- u2$count*(sum(loxcoder::data(rep1)$count)/sum(loxcoder::data(rep2)$count))
  return(data.frame(code = bc_union, size = u1$size, rep1_count = u1$count, rep2_count = u2$count))
}

