# dyn.load('/home/zhang.s/project_2019/edlib/build/lib/libedlib.so')

loxcoder::load_origin_distmaps('/run/media/zsteve/ssd_ext/maps/origin')
loxcoder::load_pair_distmaps('/run/media/zsteve/ssd_ext/maps/pair')

loxcoder::decode(c('/run/media/zsteve/ssd_ext/renamed/50k_1_1__B_S10_R1_001.fastq', '/run/media/zsteve/ssd_ext/renamed/50k_1_1__B_S10_R2_001.fastq'), name = 'test', meta = data.frame(),
                 min_r1_len = 300, min_r2_len = 280) -> t
# add cassette validity information
t <- loxcoder::impute(t)
t <- loxcoder::validate(t)
table(data(t)$is_valid) # summary of validity
# we overload for multiple types of calling arguments
loxcoder::pack(loxcoder::get_cass_vec(data(t)$code), data(t)$is_valid)
loxcoder::pack(data(t)$code, data(t)$is_valid)

t <- loxcoder::makeid(t)

loxcoder::retrieve_dist_origin(data(t)$id, data(t)$size)

t <- loxcoder::get_origin_dist(t)
loxcoder::size_plot(t)
loxcoder::dist_orig_plot(t, 13)
sizes <- list(3, 5, 7, 9)
a <- lapply(sizes, function(x) loxcoder::dist_orig_plot(t, x))

# plot some nice plots ...
grid.arrange(a[[1]], a[[2]], a[[3]], a[[4]], ncol = 1)

b <- lapply(sizes, function(x) loxcoder::rank_count_plot(t, x))
grid.arrange(b[[1]], b[[2]], b[[3]], b[[4]], ncol = 2)

u <- valid(t)
u <- u[u$size == 9, ]

loxcoder::get_cass_vec(u$code) -> x
loxcoder::retrieve_dist_pair(x[1:10], size = 9) -> m

Heatmap(m, cluster_rows = F, cluster_columns = F, col = colorRamp2(c(0, 10), c("white", "red")))


# try to compare 2 samples
files <- list(c('50k_1_1__A_S9_R1_001.fastq', '50k_1_1__A_S9_R2_001.fastq'),
             c('50k_1_1__B_S10_R1_001.fastq', '50k_1_1__B_S10_R2_001.fastq'),
             c('50k_1_2__A_S3_R1_001.fastq', '50k_1_2__A_S3_R2_001.fastq'),
             c('50k_1_2__B_S4_R1_001.fastq', '50k_1_2__B_S4_R2_001.fastq'),
             c('12k_2_1__A_S35_R1_001.fastq', '12k_2_1__A_S35_R2_001.fastq'),
             c('12k_2_1__B_S36_R1_001.fastq', '12k_2_1__B_S36_R2_001.fastq'))

path <- '/run/media/zsteve/ssd_ext/renamed'
files <- lapply(files, function(x) { c(paste(path, x[1], sep = '/'), paste(path, x[2], sep = '/'))})
files_dec <- lapply(files, function(x) {
  t <- loxcoder::decode(x, data.frame(), min_r1_len = 300, min_r2_len = 290)
  t <- loxcoder::validate(t)
  t <- loxcoder::makeid(t)
  t <- loxcoder::get_origin_dist(t)
})

lapply(files_dec, function(x) x@decode_stats)

lapply(files_dec, function(x) table(data(x)$is_valid))

lapply(files_dec, function(f) loxcoder::size_plot(f))

sizes <- list(3, 5, 7, 9)
a <- lapply(files_dec, function(f) lapply(sizes, function(x) loxcoder::dist_orig_plot(f, x)))
plot_a <- function(x) grid.arrange(a[[x]][[1]], a[[x]][[2]], a[[x]][[3]], a[[x]][[4]], ncol = 1)

b <- lapply(files_dec, function(f) lapply(sizes, function(x) loxcoder::rank_count_plot(f, x)))
plot_b <- function(x) grid.arrange(b[[x]][[1]], b[[x]][[2]], b[[x]][[3]], b[[x]][[4]], ncol = 2)

# check replicates

u <- get_comparison_table(files_dec[[5]], files_dec[[6]])
loxcoder::pair_comparison_plot(files_dec[[6]], files_dec[[5]])


loxcoder::decode(c('/run/media/zsteve/ssd_ext/renamed/50k_1_2__B_S4_R1_001.fastq', '/run/media/zsteve/ssd_ext/renamed/50k_1_2__B_S4_R2_001.fastq'), name = 'test', meta = data.frame(),
                 min_r1_len = 300, min_r2_len = 280) -> t

t <- loxcoder::impute(t)
t <- loxcoder::validate(t)
table(data(t)$is_valid) # summary of validity
# we overload for multiple types of calling arguments
loxcoder::pack(loxcoder::get_cass_vec(data(t)$code), data(t)$is_valid)
loxcoder::pack(data(t)$code, data(t)$is_valid)

t <- loxcoder::makeid(t)

loxcoder::retrieve_dist_origin(data(t)$id, data(t)$size)

t <- loxcoder::get_origin_dist(t)

u <- valid(t)
u %>% filter(size == 9 & dist_orig >= 3 & count > 1) -> u9
u %>% filter(size == 13 & dist_orig >= 3) -> u13
rbind(u9, u13) -> u_both
a <- loxcoder::get_cass_vec(u_both$code)
loxcoder::retrieve_dist_pair(a, a) -> m
colnames(m) <- u_both$code
rownames(m) <- u_both$code
h <- Heatmap(m)
ro <- row_order(h)
co <- column_order(h)
Heatmap(m[ro, co], cluster_rows = F, cluster_columns = F)
# n <- melt(m[ro, co][120:270, 120:270])
n <- melt(m[ro, co][300:320, 300:320])
n$dist_orig1 <- u_both[match(n$Var1, u_both$code), c('dist_orig')]
n$dist_orig2 <- u_both[match(n$Var2, u_both$code), c('dist_orig')]
ggplot(data = n) + geom_tile(aes(x = Var1, y = Var2, fill = factor(value), text = paste(rownames(n), "info1 = ", dist_orig1, "info2 = ", dist_orig2))) + scale_color_discrete(guide = F)
ggplotly()


loxcoder::load_origin_distmaps('/stornext/HPCScratch/home/zhang.s/project_2019/loxcodeR/maps/origin')
loxcoder::load_pair_distmaps('/stornext/HPCScratch/home/zhang.s/project_2019/loxcodeR/maps/pair')
## Try out loxcode_experiment
x <- loxcoder::load_from_xlsx('Stephen_experiment',
                              "/stornext/HPCScratch/home/zhang.s/project_2019/loxcodeR/analysis/sample_test.xlsx",
                              '/stornext/HPCScratch/home/zhang.s/project_2019/NN128_renamed/', '_R1_001.fastq', '_R2_001.fastq')

loxcoder::sample(x, '12k_3_1__A') # can access sample (as S4 object) like this
loxcoder::get_valid(x, '50k_1_1__A') # or access valid codes directly as data.frame like this
t <- loxcoder::get_valid(x, '50k_1_1__A') %>% filter(size == 9 & dist_orig > 4 & count > 1) # get codes as a data.frame
m <- loxcoder::get_pair_dist(t, t)

loxcoder::load_prob_files('/stornext/HPCScratch/home/zhang.s/project_2019/wehi-project-19/markov/out')
s <- loxcoder::sample(x, '50k_1_1__A')
s <- loxcoder::retrieve_prob_ensemble(s)

ggplot(data = filter(valid(s), dist_orig > 0 & size == 5)) + geom_point(aes(x = 1/prob, y = count)) + geom_vline(aes(xintercept = 50000)) + scale_x_log10() + scale_y_log10()

