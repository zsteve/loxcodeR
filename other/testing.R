dyn.load('/home/zhang.s/project_2019/edlib/build/lib/libedlib.so')

loxcoder::decode(c('other/reconstruct_loxcode/fq/0_S47_R1_001.fastq', 'other/reconstruct_loxcode/fq/0_S47_R2_001.fastq'), meta = data.frame(),
                 min_r1_len = 300, min_r2_len = 280) -> t
# add cassette validity information
t <- loxcoder::validate(t)
table(data(t)$is_valid) # summary of validity
# we overload for multiple types of calling arguments
loxcoder::pack(loxcoder::get_cass_vec(data(t)$code), data(t)$is_valid)
loxcoder::pack(data(t)$code, data(t)$is_valid)

t <- loxcoder::makeid(t)x

# try loading some distance maps
loxcoder::load_origin_distmaps('maps/origin')
loxcoder::load_pair_distmaps('maps/pair')
loxcoder::retrieve_dist_origin(data(t)$id, data(t)$size)

t <- loxcoder::get_origin_dist(t)
loxcoder::size_plot(t)
loxcoder::dist_orig_plot(t, 9)
sizes <- list(3, 5, 7, 9)
a <- lapply(sizes, function(x) loxcoder::dist_orig_plot(t, x))

# plot some nice plots ...
grid.arrange(a[[1]], a[[2]], a[[3]], a[[4]], ncol = 1)

b <- lapply(sizes, function(x) loxcoder::rank_count_plot(t, x))
grid.arrange(b[[1]], b[[2]], b[[3]], b[[4]], ncol = 2)

u <- valid(t)
u <- u[u$size == 9, ]

loxcoder::get_cass_vec(u$code) -> x
loxcoder::retrieve_dist_pair(x[1:100], size = 9) -> m

Heatmap(m, cluster_rows = F, cluster_columns = F, col = colorRamp2(c(0, 10), c("white", "red")))
