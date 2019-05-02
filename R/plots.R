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

setGeneric("count_plot", function(x, size) {standardGeneric("count_plot")})
#' @export
setMethod("count_plot", "loxcode_sample", function(x, size){
  u <- valid(x)
  u <- u[u$size == size, ]
  u <- u[order(u$count, decreasing = T), ]
  ggplot(data = u) + geom_point(aes(x = 1:nrow(u), y = count)) +
    scale_x_log10(breaks = 1:10, labels = u$code[1:10]) + scale_y_log10() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    xlab('code') + ggtitle(sprintf('size = %d', size))
})
