# diversity-begets-diversity
# Juan
# Jun 21, 2021

rm(list=ls())
library(dplyr)
library(ggplot2)
library(grid)
library(scales)
library(tidyr)

data <- read.table('dbd_corr_vs_params.txt',sep='\t',header=TRUE)
data <- data[!is.na(data$r),]
data_rpos <- read.table('dbd_rpos.txt',sep='\t',header=TRUE)
data_rneg <- read.table('dbd_rneg.txt',sep='\t',header=TRUE)

data_r <- rbind(cbind(data_rpos,run='pos'),
                cbind(data_rneg,run='neg'))
data_r$n <- 1
data_r <- aggregate(formula = n ~ n_families + n_species_focal + run,
                    data = data_r,
                    FUN = sum)
data_r$run <- factor(data_r$run, levels = c('pos','neg'))

data$signif <- data$p < 1e-2
sparsity_threshold <- 4
preference_threshold <- 0.5

myBootstrap <- function(x, conf = 0.05) {
  
  xmean <- mean(x)
  
  b <- rep(NA, 1000)
  for (i in 1:1000) {
    xi <- x[sample(1:length(x), size = length(x), replace = TRUE)]
    b[i] <- mean(xi)
  }
  
  b <- cumsum(table(b)/sum(table(b)))
  
  xmin <- as.numeric(names(b)[which(b>=conf)[1]])
  xmax <- as.numeric(names(b)[which(b>=(1-conf))[1] - 1])
  
  if(length(xmin) == 0) xmin <- xmean
  if(length(xmax) == 0) xmax <- xmean
  
  return(c(mean = xmean, err_min = xmin, err_max = xmax))
  
}

w <- 100
myplots <- vector(mode='list')
colors <- c('#dfc27d', '#a6611a',
            '#80cdc1', '#018571')

myplots[[1]] <-
  ggplot(data) +
    geom_jitter(aes(y=w/2,x=r,color=r),
                height=w/5,
                size=2) +
    geom_histogram(aes(x=r),
                   color='black',
                   fill='white',
                   bins = 20) +
    geom_vline(xintercept=0,
               lty=2) +
    scale_x_continuous(name='Correlation between number of non-focal families\nand number of species in focal family',
                       limits=c(-1,1)) +
    scale_y_continuous(name=NULL,
                       breaks=NULL,
                       limits=c(0,w)) +
    scale_color_gradient2(high = '#d7191c', low = '#2c7bb6', mid = 'gray') +
    coord_flip() +
    theme_bw() +
    theme(panel.grid.major=element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(),
          legend.title = element_blank(),
          aspect.ratio = 2,
          axis.text = element_text(size = 15),
          legend.position = 'none')

myplots[[2]] <-
  ggplot(data_r,aes(x=n_families,y=n_species_focal,group=run,alpha=n)) +
    geom_point() +
    geom_smooth(aes(color = run),
                method='lm',
                formula=y~x,
                se=F) +
    facet_wrap(~run,
               nrow=2,
               scales='free_y') +
    scale_x_continuous(name='Number of non-focal families') +
    scale_y_continuous(name='Number of species in focal family',
                       position='right',
                       breaks=pretty_breaks(n=4)) +
    scale_color_manual(values = c('#d7191c','#2c7bb6')) +
    expand_limits(y=0) +
    theme_bw() +
    theme(panel.grid.major=element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(),
          legend.title = element_blank(),
          aspect.ratio = 1/2,
          axis.text = element_text(size = 15),
          strip.background = element_blank(),
          strip.text = element_blank(),
          legend.position = 'none')

# categorize data according to model parameters
data_groups <- data[, c('metabolism', 'byproductSparsity', 'preferenceStrength', 'r')]
data_groups$metabolism <- gsub('common', 'Common core metabolism', data_groups$metabolism)
data_groups$metabolism <- gsub('specific', 'Family-specific metabolism', data_groups$metabolism)
data_groups$byproductSparsity <- c('Sparse secretions', 'Distributed secretions')[(data_groups$byproductSparsity < sparsity_threshold) + 1]
data_groups$preferenceStrength <- c('Strong specialists', 'Weak specialists/Generalists')[(data_groups$preferenceStrength < preference_threshold) + 1]
data_groups$group <- interaction(data_groups$metabolism,
                                 data_groups$byproductSparsity,
                                 data_groups$preferenceStrength)

data_groups$r <- as.numeric(data_groups$r > 0)

data_groups <- do.call(data.frame,
                       aggregate(formula = r ~ metabolism + byproductSparsity + preferenceStrength + group,
                                 data = data_groups,
                                 FUN = myBootstrap))
colnames(data_groups)[5:7] <- c('freq', 'freq_minus', 'freq_plus')

data_groups <- data_groups[order(data_groups$freq, decreasing = TRUE), ]
data_groups$group <- gsub('\\.', '\n', data_groups$group)
data_groups$group <- factor(data_groups$group, levels = data_groups$group)

myplots[[3]] <-
  ggplot(data_groups, aes(y = freq, x = group, group = group, fill = group)) +
    geom_errorbar(aes(ymin = freq, ymax = freq_plus),
                  width = 0.2) +
    geom_bar(stat = 'identity',
             color = 'black',
             width = 0.6) +
    scale_x_discrete(name = '') +
    scale_y_continuous(name = 'Recurrence of\npositive diversity slope',
                       limits = c(0, 1)) +
    scale_fill_brewer(palette = 'OrRd',
                      direction = -1) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = 'none',
          axis.text.x = element_text(angle = 45,
                                     hjust = 1),
          aspect.ratio = 0.5)

print(myplots[[1]])
print(myplots[[2]])
print(myplots[[3]])

ggsave(file.path('.','plots','simul_r_hist.pdf'),
       plot=myplots[[1]],
       device='pdf',
       height=90,
       width=150,
       units='mm',dpi=600)
ggsave(file.path('.','plots','simul_rpos_rneg.pdf'),
       plot=myplots[[2]],
       device='pdf',
       height=90,
       width=150,
       units='mm',dpi=600)
ggsave(file.path('.','plots','simul_r_common-vs-specific.pdf'),
       plot=myplots[[3]],
       device='pdf',
       height=90,
       width=150,
       units='mm',dpi=600)
