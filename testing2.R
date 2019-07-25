library(loxcoder)

loxcoder::load_origin_distmaps('/run/media/zsteve/ssd_ext/loxcode/maps/origin')
loxcoder::load_pair_distmaps('/run/media/zsteve/ssd_ext/loxcode/maps/pair')

# x_50k <- loxcoder::load_from_xlsx('NN128_renamed', '/home/zsteve/loxcoder/nn128_50k_samples.xlsx',
#                          dir='/home/zsteve/ssd_ext/loxcode/renamed/',
#                          suffix_R1='_R1_001.fastq', suffix_R2='_R2_001.fastq', load = T)
# save('x_50k', file = '/home/zsteve/loxcoder/NN128_50k.RData')

load('/home/zsteve/loxcoder/NN128_50k.RData')
load('/home/zsteve/loxcoder/NN128.RData')

merge_50k_by_pulse <- loxcoder::merge_by(x_50k, 'pulse')

size_plot_pulse1 <- loxcoder::dist_orig_plot(loxcoder::sample(merge_50k_by_pulse, '1'), 9)
size_plot_pulse2 <- loxcoder::dist_orig_plot(loxcoder::sample(merge_50k_by_pulse, '2'), 9)

grid.arrange(size_plot_pulse1, size_plot_pulse2, nrow = 2)

loxcoder::venn_2way(merge_50k_by_pulse, '1', '2', size_range = c(9, 9), dist_range = c(0, 10))
