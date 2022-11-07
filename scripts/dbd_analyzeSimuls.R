rm(list = ls())

library(ggplot2)

files <- list.files('../simul_runs', full.names = T)
files <- files[grepl('.txt', files)]
data <- lapply(files, read.csv, sep = '\t')

for (i in 1:length(data)) colnames(data[[i]])[1:2] <- c('family', 'species')

# abundance threshold (below which a species is considered extinct)
extinction_threshold <- 1/10000 # note that experimental data are rarified to 10000 reads

# get relative abundances
for (i in 1:length(data)) {
  total_abundance <- colSums(data[[i]][, 3:ncol(data[[i]])])
  for (col in 3:ncol(data[[i]])) {
    data[[i]][, col] <- data[[i]][, col]/total_abundance[col-2]
  }
}

# rarify abundance data to 10000 "reads" (to mimic experimental protocol)
smart.round <- function(x) {
  # function to round array while preserving sum, from:
  # https://stackoverflow.com/questions/32544646/round-vector-of-numerics-to-integer-while-preserving-their-sum
  y <- floor(x)
  indices <- tail(order(x-y), round(sum(x)) - sum(y))
  y[indices] <- y[indices] + 1
  y
}
for (i in 1:length(data)) {
  for (col in 3:ncol(data[[i]])) {
    data[[i]][, col] <- smart.round(10000*data[[i]][, col])
  }
}

### LOOP THROUGH RUNS, get DBD info

#wrapper function
getDBDdata <- function(i) { # i is the run identifier, from 0 to 99
  
  which_file_t0 <- which(basename(files) == paste('Run', i, '_T0.txt', sep = ''))
  which_file_t20 <- which(basename(files) == paste('Run', i, '_T20.txt', sep = ''))
  
  data0_i <- data[[which_file_t0]]
  data_i <- data[[which_file_t20]]
  
  # identify focal (most abundant) family in each community
  family_abundances <- as.data.frame(matrix(NA, nrow = 0, ncol = 101))
  colnames(family_abundances) <- c('family', paste('W', 0:99, sep = ''))
  family_abundances_t0 <- family_abundances
  for (fam in unique(data_i$family)) {
    family_abundances <- rbind(family_abundances,
                               cbind(family = fam,
                                     as.data.frame(t(colSums(data_i[data_i$family == fam, c(3:ncol(data_i))])))))
    family_abundances_t0 <- rbind(family_abundances_t0,
                                  cbind(family = fam,
                                        as.data.frame(t(colSums(data0_i[data0_i$family == fam, c(3:ncol(data0_i))])))))
  }
  focal_family <- setNames(rep(NA, 100),
                           paste('W', 0:99, sep = ''))
  for (col in 0:99) {
    focal_family[paste('W', col, sep = '')] <- family_abundances$family[which.max(family_abundances[, col+2])]
  }
  focal_family <- setNames(rep('F0', 100),
                           paste('W', 0:99, sep = '')) ### FIXME: F0 should always be the focal family (it is the consumer of the main resource)
  
  # get number of non-focal families in each community at times 0 and final
  nFam_t0 <- NULL
  nFam <- NULL
  for (col in 2:101) {
    nFam <- c(nFam,
              sum(family_abundances[, col] > 0 & family_abundances$family != focal_family[col-1]))
    nFam_t0 <- c(nFam_t0,
                 sum(family_abundances_t0[, col] > 0 & family_abundances_t0$family != focal_family[col-1]))
  }
  
  # get number of focal species in each community at times 0 and final
  nSpecies_t0 <- NULL
  nSpecies <- NULL
  for (col in 3:102) {
    nSpecies <- c(nSpecies,
                  sum(data_i[data_i$family != focal_family[col-2], col] > 0))
    nSpecies_t0 <- c(nSpecies_t0,
                     sum(data0_i[data0_i$family != focal_family[col-2], col] > 0))
  }
  
  # output data frame
  output_data <- rbind(data.frame(run = paste('Run_', i, sep = ''),
                                  community = paste('Community_', i, '.', 0:99, '.T0', sep = ''),
                                  time = 'T0',
                                  n_nonFocalFamilies = nFam_t0,
                                  n_focalSpecies = nSpecies_t0),
                       data.frame(run = paste('Run_', i, sep = ''),
                                  community = paste('Community_', i, '.', 0:99, '.T20', sep = ''),
                                  time = 'T20',
                                  n_nonFocalFamilies = nFam,
                                  n_focalSpecies = nSpecies))
  
  return(output_data)
  
}

# get for all runs
dbd_data <- lapply(0:99, getDBDdata)
dbd_data <- do.call(rbind, dbd_data)

# add metadata for each run
metadata <- read.csv('../simul_params/dbd_corr_vs_params.txt', sep = '\t')
metadata <- cbind(run = paste('Run_', 0:99, sep = ''),
                  metadata[, c('preferenceStrength', 'byproductSparsity', 'metabolism')])
dbd_data <- merge(metadata, dbd_data, by = 'run', all = T)

# get DBD correlations
dbd_corr <- data.frame(run = character(0),
                       time = character(0),
                       slope = numeric(0),
                       pval = numeric(0))
for (i in unique(dbd_data$run)) {
  for (t in c('T0', 'T20')) {
  
    mylm <- lm(data = dbd_data[dbd_data$run == i & dbd_data$time == t, ],
               formula = n_focalSpecies ~ n_nonFocalFamilies)
    slope <- as.numeric(mylm$coefficients[2])
    if (!is.na(slope)) {
      pval <- as.numeric(pf(summary(mylm)$fstatistic[1],
                       summary(mylm)$fstatistic[2],
                       summary(mylm)$fstatistic[3],
                       lower.tail=F))
    } else {
      pval <- NA
  }
    
    dbd_corr <- rbind(dbd_corr,
                      data.frame(run = i,
                                 time = t,
                                 slope = slope,
                                 pval = pval))
  }
}

dbd_corr <- merge(metadata, dbd_corr, by = 'run', all = T)
nan_runs <- unique(dbd_corr$run[is.na(dbd_corr$slope)])
dbd_corr <- dbd_corr[!(dbd_corr$run %in% nan_runs), ]

# set thresholds for when we consider byproducts 'sparse' and consumers 'specialists'
sparsity_threshold <- 1 # 4
preference_threshold <- 0.5

dbd_corr$preferenceStrength <- sapply(dbd_corr$preferenceStrength,
                                      FUN = function(x) {
                                        if (x >= preference_threshold) return('specialists')
                                        else return('generalists')
                                      })
dbd_corr$byproductSparsity <- sapply(dbd_corr$byproductSparsity,
                                     FUN = function(x) {
                                       if (x >= sparsity_threshold) return('sparse')
                                       else return('dense')
                                     })

dbd_corr$group <- paste(dbd_corr$metabolism, dbd_corr$byproductSparsity, dbd_corr$preferenceStrength, sep = '\n')

groupMedians <- setNames(rep(NA, length(unique(dbd_corr$group))), unique(dbd_corr$group))
for (g in unique(dbd_corr$group)) {
  groupMedians[g] <- median(dbd_corr$slope[dbd_corr$group == g & dbd_corr$time == 'T20'])
}
groupMedians <- sort(groupMedians, decreasing = T)

groupMeans <- setNames(rep(NA, length(unique(dbd_corr$group))), unique(dbd_corr$group))
for (g in unique(dbd_corr$group)) {
  groupMeans[g] <- mean(dbd_corr$slope[dbd_corr$group == g & dbd_corr$time == 'T20'])
}
groupMeans <- sort(groupMeans, decreasing = T)

dbd_corr$group <- factor(dbd_corr$group, levels = names(groupMedians))

# what would have been the expected correlation if there had been no community assembly process, just random species filtering?
# this provides a statistical null expectation against which simulation results can be compared

# richness (num. of species) in stabilized communities
n_species <- lapply(0:99,
                    FUN = function(i) {
                      
                      which_file_t0 <- which(basename(files) == paste('Run', i, '_T0.txt', sep = ''))
                      which_file_t20 <- which(basename(files) == paste('Run', i, '_T20.txt', sep = ''))
                      
                      data0_i <- data[[which_file_t0]]
                      data_i <- data[[which_file_t20]]
                      
                      data_i[, 3:102] <- data_i[, 3:102]/colSums(data_i[, 3:102])
                      
                      n_species <- sapply(3:102, FUN = function(col) sum(data_i[, col] > 1/10000))
                      
                      return(n_species)
                      
                    })
n_species <- unlist(n_species)

if (F) { # uncomment this to check that the observed DBD slopes in the simulations are larger than the slope we would expect just from random species sampling
  
  # we now *randomly* assemble 10000 communities of specified richness (matching the richness observed in the simulations)
  inoc <- read.csv('./data/inoc.csv')
  
  makeRandomCommunity <- function(target_richness) {
    
    current_richness <- 0
    community <- rep(0, nrow(inoc))
    while (current_richness < target_richness) {
      sampled_species <- sample(1:nrow(inoc), 1, prob = inoc$Relative_Abundance.avg)
      community[sampled_species] <- community[sampled_species] + 1
      current_richness <- sum(community > 0)
    }
    
    # get # of families and # of focal species in randomly assembled community
    community <- cbind(inoc[, c('Family', 'ESV')], abundance = community)
    families <- data.frame(Family = unique(community$Family))
    families$abundance <- sapply(families$Family,
                                 FUN = function(fam) sum(community$abundance[community$Family == fam]))
    focal_family <- families$Family[which.max(families$abundance)][1]
    
    n_focalSpecies <- sum(community$abundance[community$Family == focal_family] > 0)
    n_nonFocalFamilies <- sum(families$abundance[families$Family != focal_family] > 0)
    
    return(c(n_focalSpecies = n_focalSpecies,
             n_nonFocalFamilies = n_nonFocalFamilies))
    
  }
  
  null_dbd <- lapply(n_species, FUN = makeRandomCommunity)
  null_dbd <- as.data.frame(do.call(rbind, null_dbd))
  
  # get null expectation for dbd slope
  slope_null <- lm(data = null_dbd, formula = n_focalSpecies ~ n_nonFocalFamilies)$coefficients[2]
  
  ggplot(null_dbd, aes(x = n_nonFocalFamilies, y = n_focalSpecies)) +
    geom_point(alpha = 0.25) +
    geom_smooth(method = 'lm')

}

# plot
ggplot(dbd_corr[dbd_corr$time == 'T20', ], aes(x = group, y = slope)) +
  geom_boxplot(fill = NA,
               width = 0.5,
               outlier.shape = NA) +
  geom_jitter(width = 0.15) +
  scale_y_continuous(name = 'DBD slope') +
  theme_bw() +
  theme(aspect.ratio = 0.6,
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16))

ggsave('../plots/simul_slopes.pdf',
       width = 150,
       height = 150,
       units = 'mm',
       dpi = 600)

