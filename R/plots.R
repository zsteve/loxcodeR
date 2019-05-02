setGeneric("size_plot", function(x) {standardGeneric("size_plot")})

#' Plots distribution of sizes found in samples
#'
#' @export
setMethod("size_plot", "loxcode_sample", function(x){
  g <- ggplot(data = data(x)) + geom_bar(aes(as.factor(size), fill = is_valid), position = "stack", stat = 'count') + xlab('size')
  return(g)
})

setGeneric("dist_orig_plot", function(x, size) {standardGeneric("dist_orig_plot")})
#' @export
setMethod("dist_orig_plot", "loxcode_sample", function(x, size){
  u <- valid(x)
  u <- u[u$size == size, ]
  g <- ggplot(data = u) + geom_bar(aes(dist_orig)) + ggtitle(sprintf("size = %d", size))+ scale_x_continuous(breaks = 0:10, limits = c(0, 10)) +
    xlab("Distance from origin")
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
    xlab('code') + ylab('log10(count)') + ggtitle(sprintf('size = %d', size)) +
    geom_point(aes(x = 1:nrow(u), y = dist_orig*max(log10(count))/max(dist_orig)), color = 'red', alpha = 0.1) +
    scale_y_continuous(sec.axis = sec_axis(~.*1/(max(log10(u$count))/max(u$dist_orig)), name = 'distance'))
})

# setGeneric("rank_dist_plot", function(x, size) {standardGeneric("rank_dist_plot")})
#
# #' @export
# setMethod("rank_dist_plot", "loxcode_sample", function(x, size){
#   u <- valid(x)
#   u <- u[u$size == size, ]
#   u <- u[order(u$count, decreasing = T), ]
#   u$count <- u$count/sum(u$count) # use frequency instead of raw count
#   ggplot(data = u) + geom_point(aes(x = 1:nrow(u), y = dist_orig)) +
#     scale_x_log10(breaks = 1:10, labels = u$code[1:10]) +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
#     xlab('code') + ylab('frequency') + ggtitle(sprintf('size = %d', size))
# })
