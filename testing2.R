library(loxcoder)
library(plotly)

loxcoder::load_origin_distmaps('/run/media/zsteve/ssd_ext/loxcode/maps/origin')
loxcoder::load_pair_distmaps('/run/media/zsteve/ssd_ext/loxcode/maps/pair')

# x_50k <- loxcoder::load_from_xlsx('NN128_renamed', '/home/zsteve/loxcoder/nn128_50k_samples.xlsx',
#                          dir='/home/zsteve/ssd_ext/loxcode/renamed/',
#                          suffix_R1='_R1_001.fastq', suffix_R2='_R2_001.fastq', load = T, full = F)
# save('x_50k', file = '/home/zsteve/loxcoder/NN128_50k.RData')

load('/home/zsteve/loxcoder/nn128_50k.RData')

merge_50k_by_pulse <- loxcoder::merge_by(x_50k, 'pulse')

size_plot_pulse1 <- loxcoder::dist_orig_plot(loxcoder::sample(merge_50k_by_pulse, '1'), 9) + theme_classic()
size_plot_pulse2 <- loxcoder::dist_orig_plot(loxcoder::sample(merge_50k_by_pulse, '2'), 9) + theme_classic()
grid.arrange(size_plot_pulse1, size_plot_pulse2, nrow = 2)

loxcoder::venn_2way(merge_50k_by_pulse, '1', '2', size_range = c(9, 9), dist_range = c(0, 10))

merge_50k_by_rep <- loxcoder::merge_by(x_50k, 'rep')

loxcoder::venn_2way(merge_50k_by_rep, 'A', 'B', size_range = c(9, 9), dist_range = c(0, 10))

x_50k@samp_table$rep_pulse <- paste(x_50k@samp_table$rep, x_50k@samp_table$pulse, sep = ', pulse ')
merge_50k_by_pulse_and_rep <- loxcoder::merge_by(x_50k, 'rep_pulse')

dev.off()
loxcoder::venn_2way(merge_50k_by_pulse_and_rep, 'A, pulse 2', 'A, pulse 1', size_range = c(9, 9), dist_range = c(0, 10))

loxcoder::dist_count_beeswarm_plot(loxcoder::sample(merge_50k_by_pulse, '1'), 100)

## Flippy
t <- loxcoder::get_data(merge_50k_by_rep, 'A')
loxcoder::min_flip_dist(c = t$code, size = t$size, v = t$is_valid)

## Markov stuff
library(optimbase)
library(nloptr)

p_inv = function(x) 1
p_exc = function(x) 1
p_non = 0
R = lapply(c(13, 9, 7, 5, 3), function(x) get_tr_mat(x, p_inv, p_exc, p_non))

u <- loxcoder::run_markov(4, R)
lapply(u, function(x) sum(x$p))

# k=4
# ggplot() + geom_point(aes(x = u[[k]]$p, y = u_old[[k]]$p)) + scale_x_log10() + scale_y_log10() +
#   geom_abline(slope = 1, intercept = 0)

f <- function(x){

  p_inv <- function(r){
    if(r == 3){
      return(x[1]) # 1/100
    }else{
      return(x[2]) # 1
    }
  }
  p_exc <- function(r){
    if(r == 4){
      return(x[3]) # 5
    }else{
      return(x[4]) # 5
    }
  }
  p_non <- 0
  R = lapply(c(13, 9, 7, 5, 3), function(x) get_tr_mat(x, p_inv, p_exc, p_non))
  u_new <- loxcoder::run_markov(2, R)

  v <- filter(t, dist_orig <= 2 & size == 5)
  m <- merge(v, u[[3]], by = 'code', all = TRUE)
  m <- merge(m, u_new[[3]], by = 'code', all = TRUE)
  m$freq <- m$count/sum(m$count, na.rm = T)
#   return(mean((log(m$p.y) - log(m$freq))^2, na.rm = T))
  return(summary(lm(data = m, freq ~ p.y))$r.squared)
}

g <- ggplot(data = m) + geom_beeswarm(aes(x = log10(freq), y = log10(p.x), text = code), color = 'red', groupOnX = FALSE) +
  geom_beeswarm(aes(x = log10(freq), y = log10(p.y), text = code), color = 'green', groupOnX = FALSE) +
 geom_abline(slope = 1, intercept = 0)
ggplotly(g)


get_tr_mat <- function(n, p_inv, p_exc, p_non){
  M = zeros(n+1)
  for(i in 1:(n+1)){
    for(j in 1:(n+1)){
      p = 0
      if(abs(i-j) %% 2 == 1){
        # inversion
        if(abs(i-j) < 3){
          p = p_non
        }else{
          p = p_inv(abs(i-j))
        }
      }else{
        if(abs(i-j) < 4){
          p = p_non
        }else{
          p = p_exc(abs(i-j))
        }
      }
      M[i, j] = p
      M[j, i] = p
    }
  }
  return(M/(0.5*sum(M)))
}
