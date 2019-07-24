library(loxcoder)

loxcoder::load_origin_distmaps('/run/media/zsteve/ssd_ext/loxcode/maps/origin')
loxcoder::load_pair_distmaps('/run/media/zsteve/ssd_ext/loxcode/maps/pair')

x_50k <- loxcoder::load_from_xlsx('NN128_renamed', '/home/zsteve/loxcoder/nn128_50k_samples.xlsx',
                         dir='/home/zsteve/ssd_ext/loxcode/renamed/',
                         suffix_R1='_R1_001.fastq', suffix_R2='_R2_001.fastq', load = T)

load('/home/zsteve/loxcoder/NN128.RData')

s_merged <- loxcoder::merge_sample(sample(x, '50k_1_1__A'), sample(x, '50k_1_1__B'))

x_50k@samp_table <- cbind(x_50k@samp_table, rep = unlist(lapply(x_50k@samp_table$sample, function(x) strsplit(x, split = '_')[[1]][5])))
merged_out <- loxcoder::merge_by(x_50k, 'rep')
unlist(lapply(samples(merged_out), function(x) nrow(sample(merged_out, x))))
unlist(lapply(samples(x_merged), function(x) nrow(sample(merged_out, x))))
