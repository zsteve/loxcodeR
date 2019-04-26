loxcoder::decode(c('other/reconstruct_loxcode/fq/MLLAF9_tibia_L_A_S31_R1_001.fastq', 'other/reconstruct_loxcode/fq/MLLAF9_tibia_L_A_S31_R2_001.fastq'), meta = data.frame()) -> t
# add cassette validity information
t <- loxcoder::validate(t)
table(get_data(t)$is_valid) # summary of validity
