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
  fill_scale <- scale_fill_manual(breaks = 0:15, values = rep('blue', 16)) # use twice since gradient is hard to see otherwise
  g <- ggplot(data = u) + geom_bar(aes(x = dist_orig, fill = factor(dist_orig)), show.legend = F) + ggtitle(sprintf("size = %d", size))+ scale_x_continuous(breaks = 0:10, limits = c(0, 10)) +
    xlab("Distance from origin") + ylab("Diversity") + ggtitle(loxcoder::name(x)) + fill_scale
  return(g)
})

setGeneric("rank_count_plot", function(x, size, ymax, dmax) {standardGeneric("rank_count_plot")})
#' @export
setMethod("rank_count_plot", "loxcode_sample", function(x, size, ymax, dmax){
  u <- valid(x)
  u <- u[u$size == size, ]
  u <- u[order(u$count, decreasing = T), ]
  # u$count <- u$count/sum(u$count) # use frequency instead of raw count
  ggplot(data = u) + geom_point(aes(x = 1:nrow(u), y = log10(count))) +
    geom_point(aes(x = 1:nrow(u), y = dist_orig*ymax/dmax), color = 'red') +
    scale_x_log10(breaks = 1:10, labels = u$code[1:10]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    xlab('code') + ylab('log10(count)') +
    ggtitle(sprintf('%s; size = %d', loxcoder::name(x), size)) +
    scale_y_continuous(limits = c(0, ymax), sec.axis = sec_axis(~.*(dmax/ymax), name = 'distance'))
})

setGeneric("pair_comparison_plot", function(x1, x2, ...) {standardGeneric("pair_comparison_plot")})
#' @export
setMethod("pair_comparison_plot", "loxcode_sample", function(x1, x2, dist_range = NA, samp1_name= NA, samp2_name = NA){
  u <- get_comparison_table(x1, x2, dist_range)
  if(is.na(samp1_name) | is.na(samp2_name)){
    samp1_name <- loxcoder::name(x1)
    samp2_name <- loxcoder::name(x2)
  }
  nonzero_mask <- (u$rep1_count > 0 & u$rep2_count > 0)
  rep1_zero_mask <- (u$rep1_count == 0)
  rep2_zero_mask <- (u$rep2_count == 0)
  g <- ggplot() +
    geom_point(aes(y = log10(1 + u$rep1_count[nonzero_mask]) , x = log10(1 + u$rep2_count[nonzero_mask]),
              color = as.factor(u$size[nonzero_mask])), alpha = 0.25) +
    geom_quasirandom(aes(y = rep(-0.2, sum(rep1_zero_mask)), x = log10(1 + u$rep2_count[rep1_zero_mask]), color = factor(-1)), groupOnX = F, width = 0.2, alpha = 0.25) +
    geom_quasirandom(aes(y = log10(1 + u$rep1_count[rep2_zero_mask]), x = rep(-0.2, sum(rep2_zero_mask)), color = factor(-1)), groupOnX = T, width = 0.2, alpha = 0.25) +
    geom_abline() + ylab(paste('log10(reads in ', samp1_name, ')')) + xlab(paste('log10(reads in ', samp2_name, ')')) +
    ggtitle(paste(loxcoder::name(x1), 'vs', loxcoder::name(x2)))
  return(g)
})

barcode_union <- function(rep1, rep2, dist_range){
  return(unique(c(filter(loxcoder::valid(rep1), dist_orig >= dist_range[1], dist_orig <= dist_range[2])$code,
                  filter(loxcoder::valid(rep2), dist_orig >= dist_range[1], dist_orig <= dist_range[2])$code)))
}

get_barcode_stats_rep <- function(union_bc, rep){
  index <- match(union_bc, loxcoder::valid(rep)$code)
  u <- loxcoder::valid(rep)[index, ]
  u$count[is.na(u$count)] <- 0
  u$code <- union_bc
  return(u)
}

get_comparison_table <- function(rep1, rep2, dist_range){
  bc_union <- barcode_union(rep1, rep2, dist_range)
  u1 <- get_barcode_stats_rep(bc_union, rep1)
  u2 <- get_barcode_stats_rep(bc_union, rep2)
  u1$size <- ifelse(is.na(u1$size), u2$size, u1$size)
  u2$size <- u1$size
  # scale by total number of reads
  u1$count <- u1$count
  u2$count <- u2$count*(sum(loxcoder::data(rep1)$count)/sum(loxcoder::data(rep2)$count))
  return(data.frame(code = bc_union, size = u1$size, rep1_count = u1$count, rep2_count = u2$count))
}

setGeneric("dist_count_beeswarm_plot", function(x, count_threshold) {standardGeneric("dist_count_beeswarm_plot")})
#' @export
setMethod("dist_count_beeswarm_plot", "loxcode_sample", function(x, count_threshold){
  loxcoder::valid(x) %>% filter(count < count_threshold) -> y
  g <- ggplot(data = y) + geom_quasirandom(width = 0.9, groupOnX = F, aes(x = dist_orig, y = size, size = count, color = count), alpha = 0.2) +
    scale_size_area("num_reads") +
    scale_x_continuous("distance from origin", breaks = 0:15) +
    scale_y_continuous("size", breaks = c(3, 5, 7, 9, 13)) +
    scale_color_gradient(low = 'blue', high = 'red') +
    ggtitle(loxcoder::name(x))
  return(g)
})
