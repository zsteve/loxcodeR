dyn.load('/home/zhang.s/project_2019/edlib/build/lib/libedlib.so')

loxcoder::decode(c('other/reconstruct_loxcode/fq/MLLAF9_tibia_L_A_S31_R1_001.fastq', 'other/reconstruct_loxcode/fq/MLLAF9_tibia_L_A_S31_R2_001.fastq'), meta = data.frame()) -> t
# add cassette validity information
t <- loxcoder::validate(t)
table(data(t)$is_valid) # summary of validity
# we overload for multiple types of calling arguments
loxcoder::pack(loxcoder::get_cass_vec(data(t)$code), data(t)$is_valid)
loxcoder::pack(data(t)$code, data(t)$is_valid)

t <- loxcoder::makeid(t)

# try loading some distance maps
loxcoder::load_origin_distmaps('maps/origin')
loxcoder::retrieve_dist_origin(data(t)$id, data(t)$size)

t <- loxcoder::get_origin_dist(t)

