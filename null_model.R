rm(list=ls())

library(ggplot2)
library(scales)
library(reshape2)
library(RColorBrewer)
library(colorspace)

### Suppose we have a pool of species divided into $F$ families.
### There are $S_i$ species of family i in the pool.
### We define S as the total number of species in the pool:
### $ S = S_1 + S_2 + S_3 + ... + S_F $
### There are $N_{ij}$ individuals of species $j$ of family $i$.
### So there are $\sum_{i=1}^F \sum_{j=1}^{S_i} N_{ij}
### We assemble a 'random community' by sampling $n$ individuals
### from the pool.
### How does the number of species from family 1 scale with the number
### of families in the sample?

# general parameters
nF <- 50 # number of families
nS <- 100 # max. number of species per family
nN <- 1e6 # avg. number of individuals per species

###

# avg. population sizes of each family
logAvgPopSize <- sort(runif(nF,min=log10(nN)-2,max=log10(nN)+2), decreasing = TRUE)

# sample matrix N
N_ <- matrix (NA, nrow = nF, ncol = nS)
for (i in 1:nrow(N_)) {
  N_[i, ] <- floor(10^rnorm(nS,mean=logAvgPopSize[i],sd=logAvgPopSize[i]/10))
}
N_[sample(1:length(N_),size=floor(length(N_)/2))] <- 0 # randomly set to 0 half of the elements in matrix N
N_ <- t(apply(N_,1,sort,decreasing=T))
N_ <- N_[rowSums(N_)>0,colSums(N_)>0]
N_ <- rbind(N_[which.max(rowSums(N_)),],
            N_[rowSums(N_)<max(rowSums(N_)),])

N_df <- data.frame(n = as.numeric(N_[N_>0]),
                   fam = as.character(as.numeric(matrix(rep(1:nrow(N_), ncol(N_)), nrow = nrow(N_), byrow = FALSE)[N_>0])))
ggplot(N_df, aes(x = n, group = fam, fill = fam)) +
  geom_histogram(bins = 50,
                 position='stack') +
  scale_x_log10(name = 'Population size (# individuals)') +
  scale_y_continuous(name = 'Count') +
  theme_bw() +
  theme(panel.grid = element_blank())
ggplot(N_df, aes(x = 1, y = n, group = fam, fill = fam)) +
  geom_bar(position = 'stack',
           stat = 'identity')

# get pool composition and assign colors to each family
colnames(N_) <- paste('sp_',1:ncol(N_),sep='')
rownames(N_) <- paste('fam_',1:nrow(N_),sep='')

N <- sum(N_)
S <- rowSums(N_)

pool_composition <- melt(N_/sum(N_))
colnames(pool_composition) <- c('family','species','p')
pool_composition$id <- paste(#'Species ',
                             paste(gsub('fam_','',pool_composition$family),
                                   gsub('sp_','',pool_composition$species),
                                   sep='.'),
                             sep='')

mycols <- do.call('rbind',
                  lapply(1:nrow(N_),
                    FUN = function(i) {
                      
                      ncols <- sum(pool_composition$family==paste('fam_',i,sep='') & pool_composition$p>0)
                      base_color <- brewer.pal(n=nrow(N_),'Set1')[i]
                      ramp <- colorRampPalette(c(darken(base_color,amount=0.3),
                                                 lighten(base_color,amount=0.7)))
                      
                      return_this <- data.frame(family=paste('fam_',i,sep=''),
                                                species=paste('sp_',1:ncols,sep=''),
                                                color=ramp(ncols))
                      
                      return(return_this)
                      
                    }))
pool_composition <- merge(pool_composition,mycols,by=c('family','species'),all=T)
pool_composition$color <- as.character(pool_composition$color)
pool_composition$color[pool_composition$p == 0] <- 'white'

pool_composition <- pool_composition[order(as.numeric(gsub('sp_','',pool_composition$species))),]
pool_composition <- pool_composition[order(as.numeric(gsub('fam_','',pool_composition$family))),]
pool_composition$id <- factor(pool_composition$id,
                              levels=pool_composition$id)

### ANALYTIC PROBABILITIES ### -------------------------------------------------

# analytic probability and approximation of sampling a new focal species when sampling for the (m+1) time
p_up <- function(m) {
  sapply(m,
         FUN = function(m) {
           if(m>0) {
             sum(sapply(1:ncol(N_),
                        FUN = function(j) N_[1,j]/(sum(N_) - m) * prod(1 - N_[1,j]/(sum(N_) - 0:(m-1)))))
           } else {
             sum(N_[1,])/sum(N_)
           }
           })
}

# analytic probability and approximation of sampling a new family when sampling for the (m+1) time
p_right <- function(m) {
  sapply(m,
         FUN = function(m) {
           if(m>0) {
             sum(sapply(1:nrow(N_),
                        FUN = function(i) sum(N_[i,])/(sum(N_) - m) * prod(1 - sum(N_[i,])/(sum(N_) - 0:(m-1)))))
           } else {
             sum(N[2:nrow(N_)])/sum(N_)
           }
           })
}

# evaluate probabilities for various sample sizes
n <- 1:100

p <- rbind(data.frame(class='p_up',
                      color=pool_composition$color[1],
                      lty='solid',
                      n=n,
                      p=p_up(n)),
           data.frame(class='p_right',
                      color='black',
                      lty='solid',
                      n=n,
                      p=p_right(n)))

p <- p[p$p>0,]

# plot p vs m
myplot <-
  ggplot(data=p,aes(x=n,y=p,color=color,linetype=lty)) +
    geom_line() +
    scale_color_manual(values=c('red','black'),
                       name = NULL,
                       labels = c('P (step \'up\')',
                                  'P (step \'right\')')) +
    scale_linetype(guide = 'none') +
    scale_x_continuous(name=expression(paste('Step, ',italic('m')))) +
    scale_y_continuous(name='Probability') +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = c(0.6,0.85),
          aspect.ratio = 1,
          legend.background=element_rect(fill='transparent'),
          text=element_text(size=15),
          axis.text=element_text(size=15),
          axis.line=element_blank(),
          axis.ticks=element_line(size=0.25),
          panel.border=element_rect(size=0.25)) +
    guides(linetpye = 'none')
print(myplot)

# ggsave(file.path('.','plots','probs.pdf'),
#        plot=myplot,
#        device='pdf',
#        height=70,
#        width=70,
#        units='mm',dpi=600)

### SIMULATE DRAWS: SINGLE TRAJECTORY ### --------------------------------------

# function to sample a single community and keep track of nfamilies and nspeces_focal
sample_one_community <- function(n, N_) {
  
  # n     sample size
  # N_    matrix defining species pool
  
  species_pool <- expand.grid(rownames(N_),colnames(N_))
  species_pool <- matrix(paste(species_pool[,1],species_pool[,2],sep=' '),
                         nrow=nrow(N_))
  
  community <- setNames(rep(0,length(N_)),species_pool)
  nfamilies <- NULL
  nspecies_focal <- NULL
  
  for (i in 1:n) {
   p_i <- (N_ - matrix(community,nrow=nrow(N_)))/sum(N_)
   species_i <- sample(species_pool,
                       prob=p_i,
                       size=1)
   community[species_i] <- community[species_i] + 1
   community <- matrix(community,
                       nrow=nrow(N_))
   
   nfamilies <- c(nfamilies,sum(rowSums(community[2:nrow(community),])>0))
   nspecies_focal <- c(nspecies_focal,rowSums(community>0)[1])
   
   community <- setNames(c(community),species_pool)
  }
  
  return(data.frame(n=1:n,
                    nfamilies=nfamilies,
                    nspecies_focal=nspecies_focal))
  
}

# generate a single trajectory (to plot random walk), could also admit multiple ns
n_traj <- 5000
simul_traj <- data.frame(traj_id=numeric(),
                         n_final=numeric(),
                         n=numeric(),
                         nfamilies=numeric(),
                         nspecies_focal=numeric())
for (i in 1:length(n_traj)) {
  simul_traj <- rbind(simul_traj,
                      cbind(traj_id=as.character(i),
                            n_final=n_traj[i],
                            sample_one_community(n_traj[i],N_)))
}
simul_traj$n_final <- factor(as.character(simul_traj$n_final),
                             levels=as.character(sort(unique(n_traj),decreasing=T)))

# plot
myplot <-
  ggplot(data=simul_traj,
         aes(x=nfamilies,y=nspecies_focal,group=traj_id)) + # aes(x=nfamilies,y=nspecies_focal,color=n_final)) +
    geom_line(size=1,
              color='black') +
    scale_x_continuous(name='# non-focal families',
                       minor_breaks=seq(0,nrow(N_)-1,by=5),
                       limits=c(0,nrow(N_)-1)) +
    scale_y_continuous(name='# focal species',
                       minor_breaks=seq(0,ncol(N_),by=5),
                       limits=c(0,ncol(N_))) +
    theme_bw() +
    theme(aspect.ratio = 0.618,
          legend.position = 'none',
          legend.background=element_rect(fill='transparent'),
          text=element_text(size=15),
          axis.text=element_text(size=15),
          axis.line=element_blank(),
          axis.ticks=element_line(size=0.25),
          panel.border=element_rect(size=0.25))
print(myplot)

# ggsave(file.path('.','plots','traj.pdf'),
#        plot=myplot,
#        device='pdf',
#        height=90,
#        width=150,
#        units='mm',dpi=600)

### SIMULATE DRAWS: MULTIPLE DRAWS, DBD? ### -----------------------------------

# generate n runs of the random sampling process
n_communities <- 100
n <- setNames(c(sample(5:25,n_communities,replace=TRUE),
                sample(500:2500,n_communities,replace=TRUE),
                sample(50000:250000,n_communities,replace=TRUE)),
              c(rep('Small n', n_communities),
                rep('Intermediate n', n_communities),
                rep('Large n', n_communities)))

# function to sample a community
sample_communities <- function(n, N_) {
  
  # n     sample size
  # N_    matrix defining species pool
  
  species_pool <- expand.grid(rownames(N_),colnames(N_))
  species_pool <- matrix(paste(species_pool[,1],species_pool[,2],sep=' '),
                         nrow=nrow(N_))
  
  return_this <-
    sapply(n,
           FUN = function(n) {
             
             community <- setNames(rep(0,length(N_)),species_pool)
             
             for (i in 1:n) {
               p_i <- (N_ - matrix(community,nrow=nrow(N_)))/sum(N_)
               species_i <- sample(species_pool,
                                   prob=p_i,
                                   size=1)
               community[species_i] <- community[species_i] + 1
             }
             
             community <- matrix(community,
                                 nrow=nrow(N_))
             
             return(c(nfamilies=sum(rowSums(community[2:nrow(community),])>0),
                      nspecies_focal=rowSums(community>0)[1]))
             
           })
  
  return(return_this)
  
}

run <- sample_communities(n,N_)
simul_data <- data.frame(sample_size=n,
                         group=factor(names(n),levels=c('Small n','Intermediate n','Large n')),
                         nfamilies=run[1,],
                         nspecies_focal=run[2,])

# plot
myplot <-
  ggplot(data=simul_data,
         aes(x=nfamilies,y=nspecies_focal,group=group)) +
    geom_point(alpha=0.2,
               shape=16) +
    geom_smooth(method='lm',
                formula=y~x,
                size=0.5,
                color='black',
                se=FALSE) +
    facet_wrap(~group,
               nrow=1,
               scales='free') +
    scale_x_continuous(name='# non-focal families',
                       breaks=pretty_breaks()) +
    scale_y_continuous(name='# focal species',
                       breaks=pretty_breaks()) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          aspect.ratio = 0.618,
          legend.background=element_rect(fill='transparent'),
          text=element_text(size=15),
          axis.text=element_text(size=15),
          axis.line=element_blank(),
          axis.ticks=element_line(size=0.25),
          panel.border=element_rect(size=0.25),
          strip.background=element_blank(),
          strip.text=element_text(hjust=0))
print(myplot)

ggsave(file.path('.','plots','slopes-vs-n.pdf'),
       plot=myplot,
       device='pdf',
       height=90,
       width=450,
       units='mm',dpi=600)

