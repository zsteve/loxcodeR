library(loxcoder)
library(dplyr)

load('analysis/MLLAF9_invivo/data.RData')

loxcoder::load_origin_distmaps('/stornext/HPCScratch/home/zhang.s/project_2019/loxcodeR/maps/origin')
loxcoder::load_pair_distmaps('/stornext/HPCScratch/home/zhang.s/project_2019/loxcodeR/maps/pair')

# x <- loxcoder::load_from_xlsx('MLLAF9_invivo_tibia',
#                               s = '/stornext/HPCScratch/home/zhang.s/project_2019/loxcodeR/analysis/MLLAF9_invivo/sample_table_MLLAF9.xlsx',
#                               dir = '/stornext/HPCScratch/home/zhang.s/project_2019/MLLAF9_invivo/',
#                               suffix_R1 = '_R1_001.fastq',
#                               suffix_R2 = '_R2_001.fastq',
#                               load = T)

t <- loxcoder::get_valid(x, "MLLAF9_tibia_L_A") %>% filter(size == 9)
loxcoder::get_pair_dist(t[1:10, ], t[1:10, ])
