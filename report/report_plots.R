library(ggplot2)
library(gridExtra)
library(ggbeeswarm)
library(dplyr)
library(VennDiagram)
library(cowplot)
library(viridis)
library(heatmap3)


# also load merged files

loxcoder::load_origin_distmaps('/run/media/zsteve/ssd_ext/maps/origin')
loxcoder::load_pair_distmaps('/run/media/zsteve/ssd_ext/maps/pair')

# loxcoder::load_origin_distmaps('/stornext/HPCScratch/home/zhang.s/project_2019/loxcodeR/maps/origin')
# loxcoder::load_pair_distmaps('/stornext/HPCScratch/home/zhang.s/project_2019/loxcodeR/maps/pair')
loxcoder::load_prob_files('/stornext/HPCScratch/home/zhang.s/project_2019/wehi-project-19/markov/out')

# setwd('/stornext/HPCScratch/home/zhang.s/project_2019/loxcodeR/analysis')
load('50k_merged.RData')
load('Stephen_experiment_new.RData')

## Try out loxcode_experiment
# x_merged <- loxcoder::load_from_xlsx('Stephen_50k_merged',
#                               "/stornext/HPCScratch/home/zhang.s/project_2019/loxcodeR/analysis/50k_merged.xlsx",
#                               '/stornext/HPCScratch/home/zhang.s/project_2019/NN128_renamed/50k_merged_replicates/', '_R1.fastq', '_R2.fastq')

## Report plots
## 500 x 250

tot_reads <- list()
good_reads <- list()
good_reads_prop <- list()
for(i in x@samp_table$sample){
  tot_reads[i] <- x@samples[[i]]@decode_stats$tot_reads
  good_reads[i] <- sum(loxcoder::get_data(x, i)$count)
  good_reads_prop[i] <- good_reads[[i]]/tot_reads[[i]]
}

reads_info <- data.frame(sample = as.character(names(good_reads_prop)), prop = unlist(good_reads_prop))
reads_info$tot_reads <- unlist(tot_reads)
reads_info$good_reads <- unlist(good_reads)
reads_info <- reads_info[grep('50k', reads_info$sample, invert = F), ] # only look at 50k

reads_info$N <- sapply(as.character(reads_info$sample), function(x) strsplit(x, '_')[[1]][1])
reads_info$samp <- sapply(as.character(reads_info$sample), function(x) strsplit(x, '_')[[1]][2])
reads_info$pulse <- sapply(as.character(reads_info$sample), function(x) paste('Pulse', strsplit(x, '_')[[1]][3]))
reads_info$rep <- sapply(as.character(reads_info$sample), function(x) strsplit(x, '_')[[1]][5])
reads_info$diversity <- sapply(as.character(reads_info$sample), function(s) nrow(loxcoder::get_valid(x, s)))
reads_info$valid_count_sum <- sapply(as.character(reads_info$sample), function(s) sum(loxcoder::get_valid(x, s)$count))

ggplot(data = reads_info) + geom_bar(stat = 'identity', aes(x = factor(samp), y = prop, fill = factor(rep)), position = position_dodge()) +
  facet_wrap(~N + pulse, ncol = 5) +
  ylab('% informative reads') +
  xlab('sample') +
  labs(fill = 'replicate') +
  theme_classic() -> plot_read_prop

ggplot(data = reads_info) + geom_bar(stat = 'identity', aes(x = factor(samp), y = tot_reads, fill = factor(rep)), position = position_dodge()) +
  facet_wrap(~N + pulse, ncol = 5) +
  scale_y_log10(name = 'total reads') +
  xlab('sample') +
  labs(fill = 'replicate') +
  theme_classic() -> plot_read_tot

grid.arrange(plot_read_tot, plot_read_prop)

ggplot(data = reads_info) + geom_bar(stat = 'identity', aes(x = factor(samp), y = diversity, fill = factor(rep)), position = position_dodge()) +
  facet_wrap(~ pulse, ncol = 5) +
  xlab('sample') +
  labs(fill = 'replicate') +
  theme_classic()

pairwise_same <- loxcoder::pair_comparison_plot(loxcoder::sample(x, '50k_1_1__A'), loxcoder::sample(x, '50k_1_1__B')) +
  theme_classic() +
  scale_color_grey() +
  theme(legend.position = 'none') +
  ggtitle('Pulse 1 S1(A) vs S1(B)')
pairwise_diff <- loxcoder::pair_comparison_plot(loxcoder::sample(x, '50k_1_1__A'), loxcoder::sample(x, '50k_2_1__A')) +
  theme_classic() +
  scale_color_grey() +
  theme(legend.position = 'none') +
  ggtitle('Pulse 1 S1(A) vs S2(A)')
grid.arrange(pairwise_same, pairwise_diff, ncol = 2)

size_plot_pulse1 <- loxcoder::size_plot(loxcoder::sample(x_merged, '50k_1_1')) + ylim(0, 2500) + theme_classic() +
  ggtitle('50k_1_1 (Pulse 1)') +
  theme(legend.position = c(.85, .85)) +
  scale_fill_discrete(name = '', labels = c('invalid', 'valid'))

size_plot_pulse2 <- loxcoder::size_plot(loxcoder::sample(x_merged, '50k_1_2')) + ylim(0, 2500) + theme_classic() +
  ggtitle('50k_1_2 (Pulse 2)') +
  theme(legend.position = "none")

grid.arrange(size_plot_pulse1, size_plot_pulse2, nrow = 1)

count_plot_pulse1 <- ggplot(data = loxcoder::get_data(x_merged, '50k_1_1')) + geom_bar(stat = 'identity', aes(x = factor(size), y = count, fill = factor(is_valid))) +
  theme_classic() +
  ggtitle('50k_1_1 (Pulse 1)') +
  ylim(0, 1e5) +
  theme(legend.position = "none") +
  xlab('size') + ylab('reads')

count_plot_pulse2 <- ggplot(data = loxcoder::get_data(x_merged, '50k_1_2')) + geom_bar(stat = 'identity', aes(x = factor(size), y = count, fill = factor(is_valid))) +
  theme_classic() +
  ggtitle('50k_1_2 (Pulse 2)') +
  ylim(0, 1e5) +
  theme(legend.position = "none") +
  xlab('size') + ylab('reads')

grid.arrange(count_plot_pulse1, count_plot_pulse2, nrow = 1)

dist_plot_pulse1 <- loxcoder::dist_orig_plot(loxcoder::sample(x_merged, '50k_1_1'), 9) + theme_classic() + ylim(0, 800) +
  ggtitle('50k_1_1 (Pulse 1)')
dist_plot_pulse2 <- loxcoder::dist_orig_plot(loxcoder::sample(x_merged, '50k_1_2'), 9) + theme_classic() + ylim(0, 800) +
  ggtitle('50k_1_2 (Pulse 2)')

grid.arrange(dist_plot_pulse1, dist_plot_pulse2, ncol = 2)

# ggplot(data = loxcoder::get_valid(x_merged, '50k_1_1')) +
#   geom_bar(aes(x = factor(dist_orig), fill = factor(dist_orig))) +
#   scale_fill_viridis(discrete = T)
#
# ggplot(data = loxcoder::get_valid(x_merged, '50k_1_2')) +
#   geom_bar(aes(x = factor(dist_orig), fill = factor(dist_orig))) +
#   scale_fill_viridis(discrete = T)

rank_count_pulse1 <- loxcoder::rank_count_plot(loxcoder::sample(x_merged, '50k_1_1'), 9, 4, 10) +
  scale_x_log10() +
  xlab('barcode rank') +
  ggtitle('50k_1_1 (Pulse 1)') + theme_classic() + theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

rank_count_pulse2 <- loxcoder::rank_count_plot(loxcoder::sample(x_merged, '50k_1_2'), 9, 4, 10) +
  scale_x_log10() +
  xlab('barcode rank') +
  ggtitle('50k_1_2 (Pulse 2)') + theme_classic() + theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

grid.arrange(rank_count_pulse1,
             rank_count_pulse2, nrow = 1)

# Venn diagrams
# use 800x400
venn_pulse1_nofilter <- loxcoder::venn_3way(x_merged, '50k_1_1', '50k_2_1', '50k_3_1', size_range = c(0, 13), dist_range = c(0, 15), labels = c('S1', 'S2', 'S3'))
venn_pulse2_nofilter <- loxcoder::venn_3way(x_merged, '50k_1_2', '50k_2_2', '50k_3_2', size_range = c(0, 13), dist_range = c(0, 15), labels = c('S1', 'S2', 'S3'))

plot_grid(venn_pulse1_nofilter, venn_pulse2_nofilter, nrow = 1, labels = c('Pulse 1', 'Pulse 2'), scale = 0.8)

venn_pulse1_filter <- loxcoder::venn_3way(x_merged, '50k_1_1', '50k_2_1', '50k_3_1', size_range = c(7, 13), dist_range = c(4, 15), labels = c('S1', 'S2', 'S3'))
venn_pulse2_filter <- loxcoder::venn_3way(x_merged, '50k_1_2', '50k_2_2', '50k_3_2', size_range = c(7, 13), dist_range = c(4, 15), labels = c('S1', 'S2', 'S3'))

plot_grid(venn_pulse1_filter, venn_pulse2_filter, nrow = 1, labels = c('Pulse 1', 'Pulse 2'), scale = 0.8)

# Common codes
m <- loxcoder::get_valid(x_merged, x_merged@samp_table$sample[1])
for(i in x_merged@samp_table$sample[-1]){
  m <- merge(m, loxcoder::get_valid(x_merged, i), by = c('code', 'dist_orig', 'size', 'id', 'is_valid'), all = T)
}
rownames(m) <- m$code
names(m)[6:11] <- x_merged@samp_table$sample

m_pulse1 <- m[, c('50k_1_1', '50k_2_1', '50k_3_1')]
m_pulse1 <- m_pulse1[rowSums(!is.na(m_pulse1)) > 0,  ]
codes_nonuniq_pulse1 <- rownames(m_pulse1)[rowSums(is.na(m_pulse1)) < 2]
codes_uniq_pulse1 <- rownames(m_pulse1)[rowSums(is.na(m_pulse1)) == 2]

codes_nonuniq_pulse1_dist <- m[codes_nonuniq_pulse1, ]
codes_uniq_pulse1_dist <- m[codes_uniq_pulse1, ]

codes_nonuniq_pulse1_dist$class <- 'nonunique'
codes_uniq_pulse1_dist$class <- 'unique'

nonuniq_summary <- data.frame(table(codes_nonuniq_pulse1_dist$dist_orig))
nonuniq_summary$class <- 'nonunique'
nonuniq_summary$Var1 <- as.numeric(nonuniq_summary$Var1)
uniq_summary <- data.frame(table(codes_uniq_pulse1_dist$dist_orig))
uniq_summary$class <- 'unique'
uniq_summary$Var1 <- as.numeric(uniq_summary$Var1)

plot_uniq <- ggplot(data = codes_uniq_pulse1_dist) +
  geom_bar(stat = 'count', aes(x = dist_orig, fill = factor(dist_orig))) + ylim(0, 1200) + ylab('diversity') +
  scale_x_continuous(name = 'distance', breaks = 0:10, labels = 0:10, limits = c(0, 10)) +
  scale_fill_manual(breaks = 0:15, values = c(viridis(8), viridis(8))) +
  ggtitle('Pulse 1 unique') + theme_classic() +
  theme(legend.position = 'none')

plot_nonuniq <- ggplot(data = codes_nonuniq_pulse1_dist) +
  geom_bar(stat = 'count', aes(x = dist_orig, fill = factor(dist_orig))) + ylim(0, 1200) + ylab('diversity') +
  scale_x_continuous(name = 'distance', breaks = 0:10, labels = 0:10, limits = c(0, 10)) +
  scale_fill_manual(breaks = 0:15, values = c(viridis(8), viridis(8))) +
  ggtitle('Pulse 1 non-unique') + theme_classic() +
  theme(legend.position = 'none')

# 500x250 for this figure
grid.arrange(plot_uniq, plot_nonuniq, nrow = 1)

m_pulse2 <- m[, c('50k_1_2', '50k_2_2', '50k_3_2')]
m_pulse2 <- m_pulse2[rowSums(!is.na(m_pulse2)) > 0,  ]
codes_nonuniq_pulse2 <- rownames(m_pulse2)[rowSums(is.na(m_pulse2)) < 2]
codes_uniq_pulse2 <- rownames(m_pulse2)[rowSums(is.na(m_pulse2)) == 2]

codes_nonuniq_pulse2_dist <- m[codes_nonuniq_pulse2, ]
codes_uniq_pulse2_dist <- m[codes_uniq_pulse2, ]

codes_nonuniq_pulse2_dist$class <- 'nonunique'
codes_uniq_pulse2_dist$class <- 'unique'

nonuniq_summary <- data.frame(table(codes_nonuniq_pulse2_dist$dist_orig))
nonuniq_summary$class <- 'nonunique'
nonuniq_summary$Var1 <- as.numeric(nonuniq_summary$Var1)
uniq_summary <- data.frame(table(codes_uniq_pulse2_dist$dist_orig))
uniq_summary$class <- 'unique'
uniq_summary$Var1 <- as.numeric(uniq_summary$Var1)

plot_uniq <- ggplot(data = codes_uniq_pulse2_dist) +
  geom_bar(stat = 'count', aes(x = dist_orig, fill = factor(dist_orig))) + ylim(0, 3500) + ylab('diversity') +
  scale_x_continuous(name = 'distance', breaks = 0:10, labels = 0:10, limits = c(0, 10)) +
  scale_fill_manual(breaks = 0:15, values = c(viridis(8), viridis(8))) +
  ggtitle('Pulse 2 unique') + theme_classic() +
  theme(legend.position = 'none')

plot_nonuniq <- ggplot(data = codes_nonuniq_pulse2_dist) +
  geom_bar(stat = 'count', aes(x = dist_orig, fill = factor(dist_orig))) + ylim(0, 3500) + ylab('diversity') +
  scale_x_continuous(name = 'distance', breaks = 0:10, labels = 0:10, limits = c(0, 10)) +
  scale_fill_manual(breaks = 0:15, values = c(viridis(8), viridis(8))) +
  ggtitle('Pulse 2 non-unique') + theme_classic() +
  theme(legend.position = 'none')

# 250x400 for this figure
grid.arrange(plot_uniq, plot_nonuniq, nrow = 1)

# t <- loxcoder::get_valid(x_merged, '50k_3_1') %>% filter(size == 9)
# t$prob <- loxcoder::retrieve_prob(t$id, t$size, t$dist_orig)
# t$is_common <- t$code %in% common
#
# ggplot(data = filter(t, dist_orig == 5)) + geom_boxplot(aes(x = is_common, y = -log10(prob))) +
#   geom_beeswarm(aes(x = is_common, y = -log10(prob), size = count))
# ggplotly()

# Barcode vs sample heatmap
m <- filter(loxcoder::get_valid(x_merged, x_merged@samp_table$sample[1]), dist_orig > 4 & count > 1)
for(i in x_merged@samp_table$sample[-1]){
  m <- merge(m, filter(loxcoder::get_valid(x_merged, i), dist_orig > 4 & count > 1), by = c('code', 'dist_orig', 'size', 'id', 'is_valid'), all = T)
}
rownames(m) <- m$code
m <- m[, grepl('count', names(m))]
names(m) <- x_merged@samp_table$sample
m[is.na(m)] <- 0
m <- sweep(m, 2, colSums(m), '/')

library(heatmap3)
heatmap3(m, col = viridis(130))

# MDS plot
d <- dist(t(m))
fit <- cmdscale(d, eig = T, k = 2)
mds_points <- data.frame(fit$points)
names(mds_points) <- c('x', 'y')
mds_points$sample <- sapply(rownames(mds_points), function(x) strsplit(x, split = '_')[[1]][2])
ggplot(data = mds_points) + geom_point(aes(x = x, y = y, color = factor(sample))) +
  geom_text(aes(x = x, y = y, label = rownames(mds_points)), nudge_x = 0.01, nudge_y = 0.01)


# duplicate code from above but we do it for non-merged

samples <- x@samp_table$sample
samples <- samples[grep('50k', samples)]
m <- filter(loxcoder::get_valid(x, samples[1]), dist_orig >= 4 & count > 1)
for(i in samples[-1]){
  m <- merge(m, filter(loxcoder::get_valid(x, i), dist_orig >= 4 & count > 1), by = c('code', 'dist_orig', 'size', 'id', 'is_valid'), all = T)
}
rownames(m) <- m$code
m <- m[, grepl('count', names(m))]
names(m) <- samples
m[is.na(m)] <- 0
m <- sweep(m, 2, colSums(m), '/')
# m[is.na(m)] <- 0
library(heatmap3)

pdf('out.pdf', width = 6, height = 6)
heatmap3(t(m), col = viridis(130), labCol = NA, Rowv = NA, scale = 'column')
dev.off()

# MDS plot
d <- dist(t(m))
fit <- cmdscale(d, eig = T, k = 2)
mds_points <- data.frame(fit$points)
names(mds_points) <- c('x', 'y')
mds_points$sample <- sapply(rownames(mds_points), function(x) strsplit(x, split = '_')[[1]][2])
mds_points$rep <- sapply(rownames(mds_points), function(x) strsplit(x, split = '_')[[1]][5])
ggplot(data = mds_points) + geom_text(aes(x = x, y = y, label = rep, color = factor(sample)), nudge_x = 0.01, nudge_y = 0.01) +
  scale_color_discrete(name = 'sample', labels = c('1', '2', '3'), ) + theme_classic()


###

samples <- list('50k_1', '50k_2', '50k_3')
names(samples) <- samples

old_codes <- lapply(samples, function(x){
  pulse1 <- paste0(x, '_1')
  pulse2 <- paste0(x, '_2')
  pulse1_codes <- filter(loxcoder::get_valid(x_merged, pulse1), size == 9)$code
  pulse2_codes <- filter(loxcoder::get_valid(x_merged, pulse2), size == 9)$code
  return(intersect(pulse1_codes, pulse2_codes))
})

new_codes <- lapply(samples, function(x){
  pulse1 <- paste0(x, '_1')
  pulse2 <- paste0(x, '_2')
  pulse1_codes <- filter(loxcoder::get_valid(x_merged, pulse1), size == 9)$code
  pulse2_codes <- filter(loxcoder::get_valid(x_merged, pulse2), size == 9)$code
  return(pulse2_codes[!(pulse2_codes %in% intersect(pulse1_codes, pulse2_codes))])
})

# now filter codes that are repeated (common codes)
common_codes <- intersect(old_codes[[1]], intersect(old_codes[[2]], old_codes[[3]]))
common_codes <- c(common_codes, intersect(new_codes[[1]], intersect(new_codes[[2]], new_codes[[3]])))
old_codes <- lapply(old_codes, function(x){
  return(x[!(x %in% common_codes)])
})

new_codes <- lapply(new_codes, function(x){
  return(x[!(x %in% common_codes)])
})

old <- lapply(samples, function(x) {
  s <- paste0(x, '_2')
  t <- loxcoder::get_valid(x_merged, s) %>% filter(code %in% old_codes[[x]], dist_orig >= 3 & dist_orig <= 5, count > 10)
  t$sample <- x
  return(t)
})

new <- lapply(samples, function(x) {
  s <- paste0(x, '_2')
  t <- loxcoder::get_valid(x_merged, s) %>% filter(code %in% new_codes[[x]], dist_orig >= 4, count > 1)
  t$sample <- x
  return(t)
})

lapply(old, nrow)
lapply(new, nrow)

old_all <- rbind(old[[1]], old[[2]], old[[3]])
new_all <- rbind(new[[1]], new[[2]], new[[3]])

m <- loxcoder::get_pair_dist(old_all, new_all)
n <- m
n[n > 1] = NA

row_mask <- rowSums(!is.na(n)) > 0
col_mask <- colSums(!is.na(n)) > 0
n <- n[row_mask, col_mask]
n[is.na(n)] <- 3

row_cols <- recode(old_all$sample, `50k_1` = 'red', `50k_2` = 'green', `50k_3` = 'blue')
col_cols <- recode(new_all$sample, `50k_1` = 'red', `50k_2` = 'green', `50k_3` = 'blue')
heatmap3(n, Rowv = NA, Colv = NA, scale = 'none', col = viridis(130, direction = -1),
         RowSideColors = row_cols[row_mask], ColSideColors = col_cols[col_mask],
         RowSideLabs = c('old'), ColSideLabs = c('new'))

# old_new_data <- lapply(as.list(1:3), function(x){
#   return(lapply(as.list(1:3), function(y){
#     m <- loxcoder::get_pair_dist(old[[x]], new[[y]])
#     return(m)
#   }))
# })
#
# old_new_hm <- lapply(old_new_data, function(x) lapply(x, function(m){
#   m[m > 1] <- 2
#   return(m)
# }))
#
# grab_grob <- function(){
#   grid.echo()
#   grid.grab()
# }
#
# library(gplots)
# library(gridGraphics)
#
# gl <- lapply(old_new_hm, function(x) lapply(x, function(m){
#   heatmap3(m, scale = 'none', col = viridis(160, direction = -1), Rowv = NA, Colv = NA)
#   g <- grab_grob()
#   return(g)
# }))
#
# grid.arrange(grobs = c(gl[[1]], gl[[2]], gl[[3]]), ncol = 3)

t <- loxcoder::get_valid(x_merged, '50k_1_2') %>% filter(size == 9, dist_orig >= 3)
m <- loxcoder::get_pair_dist(t, t)
