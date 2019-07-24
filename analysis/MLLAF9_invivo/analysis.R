library(loxcoder)
library(dplyr)
library(ComplexHeatmap)
library(plotly)

# loxcoder::load_origin_distmaps('/stornext/HPCScratch/home/zhang.s/project_2019/loxcodeR/maps/origin')
# loxcoder::load_pair_distmaps('/stornext/HPCScratch/home/zhang.s/project_2019/loxcodeR/maps/pair')
# loxcoder::load_prob_files('/stornext/HPCScratch/home/zhang.s/project_2019/wehi-project-19/markov/out')

loxcoder::load_origin_distmaps('/home/zsteve/ssd_ext/loxcode/maps/origin')
loxcoder::load_pair_distmaps('/home/zsteve/ssd_ext/loxcode/maps/pair')

x <- loxcoder::load_from_xlsx('MLLAF9_invivo_tibia',
                              s = '/stornext/HPCScratch/home/zhang.s/project_2019/loxcodeR/analysis/MLLAF9_invivo/sample_table_MLLAF9.xlsx',
                              dir = '/stornext/HPCScratch/home/zhang.s/project_2019/MLLAF9_invivo/',
                              suffix_R1 = '_R1_001.fastq',
                              suffix_R2 = '_R2_001.fastq',
                              load = T)


x_merged <- loxcoder::load_from_xlsx('MLLAF9_invivo_tibia_merged',
                              s = '/stornext/HPCScratch/home/zhang.s/project_2019/loxcodeR/analysis/MLLAF9_invivo/MLLAF9_merged.xlsx',
                              dir = '/stornext/HPCScratch/home/zhang.s/project_2019/MLLAF9_invivo/merged/',
                              suffix_R1 = '_R1.fastq',
                              suffix_R2 = '_R2.fastq',
                              load = T)


t <- loxcoder::get_valid(x_merged, "MLLAF9_tibia_L") %>% filter(size == 9 & dist_orig == 3)
p <- loxcoder::retrieve_prob(t$id, t$size, t$dist_orig)
t$prob <- p/sum(p)
t$empirical_prob <- t$count/sum(t$count)
ggplot(data = t) + geom_point(aes(x = prob, y = count, text = code, color = factor(dist_orig), size = size)) + scale_x_log10() + scale_y_log10()
ggplotly()

common <- intersect(filter(loxcoder::get_valid(x_merged, 'MLLAF9_tibia_L'), )$code,
                    filter(loxcoder::get_valid(x_merged, 'MLLAF9_tibia_R'), )$code)

t <- loxcoder::get_valid(x_merged, 'MLLAF9_tibia_R')
t$prob <- loxcoder::retrieve_prob(t$id, t$size, t$dist_orig)
t$is_common <- t$code %in% common

ggplot(data = filter(t, dist_orig >= 5)) + geom_boxplot(aes(x = is_common, y = -log10(prob))) +
  geom_beeswarm(aes(x = is_common, y = -log10(prob), size = count, color = factor(dist_orig), text = code))
ggplotly()

