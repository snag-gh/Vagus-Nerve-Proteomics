#Proteomics analysis helper functions

density_plot <- function(qnt, title, log=FALSE, facet=TRUE) {
  datamelt <- melt(exprs(qnt))
  colnames(datamelt) <- c("PSM", "Sample", "abundance")
  datamelt$Region <- pData(qnt)[match(datamelt$Sample, qnt$Sample), "Region"]
  datamelt$RegionLR <- pData(qnt)[match(datamelt$Sample, qnt$Sample), "RegionLR"]
  datamelt$Region <- as.factor(datamelt$Region)
  if(facet) {
    if(log) {
      ggplot(datamelt, aes(x = log2(abundance), col = Sample)) + geom_density() + facet_grid(rows = vars(Region)) + geom_vline(xintercept = 3.5) + ggtitle(title) + theme_classic()
    } else {
      ggplot(datamelt, aes(x = abundance, col = Sample)) + geom_density() + facet_grid(rows = vars(Region)) +  geom_vline(xintercept = 3.5) + ggtitle(title) + theme_classic()
    }
  } else {
    if(log) {
      ggplot(datamelt, aes(x = log2(abundance), col = Sample)) + geom_density() + geom_vline(xintercept = 3.5) + ggtitle(title) + theme_classic()
    } else {
      ggplot(datamelt, aes(x = abundance, col = Sample)) + geom_density() + geom_vline(xintercept = 3.5) + ggtitle(title) + theme_classic()
    }
  }
  
}

## Data filtering function
filter_valids <- function(prot.obj, groups, min_count, at_least_one = FALSE) {
  # prot.obj = a MSnSet object
  # groups = a character vector dictating the grouping
  # min_count = a numeric vector of the same length as "groups" indicating the minimum 
  #     number of valid values for each group for retention
  # at_least_one = TRUE means to keep the row if min_count is met for at least one group
  #     FALSE means min_count must be met across all conditions for retention
  
  group.names <- lapply(groups, function(x) prot.obj$Sample[prot.obj$Region == x])
  
  group.filter = sapply(1:length(group.names), function(i) {
    df2 = exprs(prot.obj)[ ,group.names[[i]]]   # Extract columns of interest
    df2 = as.matrix(df2)   # Cast as matrix for the following command
    sums = rowSums(!is.na(df2)) # count the number of valid values for each condition
    sums >= min_count[i]   # Calculates whether min_count requirement is met
  })
  
  if (at_least_one) {
    fData(prot.obj)$KEEP = apply(group.filter, 1, any)
  } else {
    fData(prot.obj)$KEEP = apply(group.filter, 1, all)
  }
  
  return(prot.obj)  # No rows are omitted, filter rules are listed in the KEEP column
}

filter_valids_LR <- function(prot.obj, groups, min_count, at_least_one = FALSE) {
  # prot.obj = a MSnSet object
  # groups = a character vector dictating the grouping
  # min_count = a numeric vector of the same length as "groups" indicating the minimum 
  #     number of valid values for each group for retention
  # at_least_one = TRUE means to keep the row if min_count is met for at least one group
  #     FALSE means min_count must be met across all conditions for retention
  
  group.names <- lapply(groups, function(x) prot.obj$Sample[prot.obj$RegionLR == x])
  
  group.filter = sapply(1:length(group.names), function(i) {
    df2 = exprs(prot.obj)[ ,group.names[[i]]]   # Extract columns of interest
    df2 = as.matrix(df2)   # Cast as matrix for the following command
    sums = rowSums(!is.na(df2)) # count the number of valid values for each condition
    sums >= min_count[i]   # Calculates whether min_count requirement is met
  })
  
  if (at_least_one) {
    fData(prot.obj)$KEEP = apply(group.filter, 1, any)
  } else {
    fData(prot.obj)$KEEP = apply(group.filter, 1, all)
  }
  
  return(prot.obj)  # No rows are omitted, filter rules are listed in the KEEP column
}
