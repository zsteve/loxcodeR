loxcoder::decode(c('other/reconstruct_loxcode/fq/MLLAF9_tibia_L_A_S31_R1_001.fastq', 'other/reconstruct_loxcode/fq/MLLAF9_tibia_L_A_S31_R2_001.fastq'), meta = data.frame()) -> t
# add cassette validity information
t <- loxcoder::validate(t)
table(get_data(t)$is_valid) # summary of validity
# we overload for multiple types of calling arguments
loxcoder::pack(loxcoder::get_cass_vec(get_data(t)$code))
loxcoder::pack(get_data(t)$code)

t <- loxcoder::makeid(t)

# try loading some distance maps
load_origin_distmaps('other/tables/origin/0.gz', 'other/tables/origin/1.gz', 'other/tables/origin/2.gz', 'other/tables/origin/3.gz', 'other/tables/origin/4.gz')
