library(loxcoder)

loxcoder::load_origin_distmaps('/run/media/zsteve/ssd_ext/loxcode/maps/origin')
loxcoder::load_pair_distmaps('/run/media/zsteve/ssd_ext/loxcode/maps/pair')

x_50k <- loxcoder::load_from_xlsx('NN128_renamed', '/home/zsteve/loxcoder/nn128_50k_samples.xlsx',
                         dir='/home/zsteve/ssd_ext/loxcode/renamed/',
                         suffix_R1='_R1_001.fastq', suffix_R2='_R2_001.fastq', load = T)

load('/home/zsteve/loxcoder/NN128.RData')

s_merged <- loxcoder::merge_sample(sample(x_50k, '50k_1_2__A'), sample(x_50k, '50k_1_2__B'))
nrow(loxcoder::valid(s_merged))

merged_out <- loxcoder::merge_by(x_50k, 'merged_sample')
nrow(valid(sample(merged_out, '50k_1_2')))
