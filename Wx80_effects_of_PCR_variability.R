###############################################################################
# Analysis of variability among PCRs, and effects of multiple PCRs on 
# biodiversity estimates for different genes, extraction methods, and samples
###############################################################################

library(data.table)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(reshape2)
library(scales)
library(vegan)
library(vegetarian)

# Set ggplot theme
theme_set(theme_bw(base_size=8))

# Set working directory
setwd("./")

###############################################################################
# Calculate OTU occurrences and abundances among sets of ten PCRs,
# and make a barchart of the results
###############################################################################

# Compile data on occurrence and abundance of OTUs among sets of ten PCRs
genes <- c("16S","18S","26S","COI")

all_res <- data.frame() # Empty data frame for all results
for (gene in genes) {
  #f <- Sys.glob(paste0("Wx80_OTU_tables/Wx80_*", gene, "*_MEE1.0_min1_OTUtable.txt"))
  f <- Sys.glob(paste0("Wx80_OTU_tables/Wx80_*", gene, "*_MEE1.0_min2_OTUtable.txt"))
  OTUtable <- read.table(f, header = TRUE, row.names = 1, sep = "\t")
  print(paste("analysing", gene, "..."))
  
  # Calculate OTU occurrence and abundance data for each set of ten PCRs
  samples <- c("D_Psoil","D_Pmax","D_Ind","D_PO4","P_Psoil","P_Pmax","P_Ind","P_PO4")
  res <- data.frame() # Empty data frame for results for each set of PCRs
  
  for(s in samples){
    OTUtab <- OTUtable[, grep(s, colnames(OTUtable))]
    OTUtab <- OTUtab[rowSums(OTUtab) > 0, ] # Drop empty rows
    # Get a data frame of OTU occurrences among PCRs
    x <- as.data.frame(apply(OTUtab, 1, function(x) sum(x > 0, na.rm = TRUE)))
    colnames(x) <- "occs"
    # Add a column with abundances of each OTU
    x$sums <- rowSums(OTUtab)
    
    # Get a summary of OTU occurrences
    y <- aggregate(x$occs, by = list(x$occs), length)
    colnames(y) <- c("occurrence", "OTUs")
    y$prop_OTUs <- y$OTUs/sum(y$OTUs)
    y <- melt(y, id.vars = "occurrence")
    
    # Get a summary of OTU abundances
    z <- aggregate(x$sums, by = list(x$occs), sum)
    colnames(z) <- c("occurrence", "sequences")
    z$prop_sequences <- z$sequences/sum(z$sequences)
    z <- melt(z, id.vars = "occurrence")
    
    # Combine the occurrence and abundance data
    yz <- rbind(y, z)
    yz$sample <- s
    yz$gene <- gene
    res <- rbind(res, yz)
  }
  all_res <- rbind(all_res, res)
}

# Make a nice barchart of OTU abundances and occurrences among sets of ten PCRs
head(all_res)

# Fix method names and factors
all_res$subplot <- paste("Subplot", sapply(strsplit(all_res$sample, "_"), "[[", 1))
all_res$method <- gsub("D_|P_", "", all_res$sample)
all_res$method <- gsub("Psoil", "1.5 g\ndirect", all_res$method)
all_res$method <- gsub("Pmax", "7.5 g\ndirect", all_res$method)
all_res$method <- gsub("Ind", "15 g\nindirect", all_res$method)
all_res$method <- gsub("PO4", "15 g PO4\nbuffer", all_res$method)
all_res$method <- factor(all_res$method, ordered = TRUE, 
                         levels = c("1.5 g\ndirect", "7.5 g\ndirect", 
                                    "15 g\nindirect", "15 g PO4\nbuffer"))
all_res$variable <- gsub("prop_(OTUs|sequences)", "\\1 \\(proportion\\)", all_res$variable)
all_res$occurrence <- as.character(all_res$occurrence)
all_res$occurrence <- factor(all_res$occurrence, levels = c("1","2","3","4","5","6","7","8","9","10"), ordered = TRUE)

# Set up an RColorBrewer colour palette
colors <- brewer.pal(8, "Spectral")
pal <- colorRampPalette(colors) 

# Make a separate graph for each gene
p1 <- ggplot(all_res[all_res$gene == "16S", ], aes(x = method, y = value, fill = occurrence)) + geom_bar(stat = "identity") + 
  facet_grid(variable ~ subplot, scales = "free") +
  scale_y_continuous(labels = comma) + ggtitle("16S") +
  xlab("Subplot / DNA extraction method") + ylab(paste("OTUs / sequence reads")) +
  theme(strip.background = element_blank()) +
  scale_fill_manual(values = pal(10))

p2 <- ggplot(all_res[all_res$gene == "18S", ], aes(x = method, y = value, fill = occurrence)) + geom_bar(stat = "identity") + 
  facet_grid(variable ~ subplot, scales = "free") +
  scale_y_continuous(labels = comma) + ggtitle("18S") +
  xlab("Subplot / DNA extraction method") + ylab(paste("OTUs / sequence reads")) +
  theme(strip.background = element_blank()) + 
  scale_fill_manual(values = pal(10))

p3 <- ggplot(all_res[all_res$gene == "26S", ], aes(x = method, y = value, fill = occurrence)) + geom_bar(stat = "identity") + 
  facet_grid(variable ~ subplot, scales = "free") +
  scale_y_continuous(labels = comma) + ggtitle("26S") +
  xlab("Subplot / DNA extraction method") + ylab(paste("OTUs / sequence reads")) +
  theme(strip.background = element_blank()) + 
  scale_fill_manual(values = pal(10))

p4 <- ggplot(all_res[all_res$gene == "COI", ], aes(x = method, y = value, fill = occurrence)) + geom_bar(stat = "identity") + 
  facet_grid(variable ~ subplot, scales = "free") +
  scale_y_continuous(labels = comma) + ggtitle("COI") +
  xlab("Subplot / DNA extraction method") + ylab(paste("OTUs / sequence reads")) +
  theme(strip.background = element_blank()) + 
  scale_fill_manual(values = pal(10))

# Combine all four plots and output to pdf. Use Illustrator to remove redundant labels and legends. 
#pdf(file = "Wx80_OTU_occurrence_vs_abundance_boxplot_all_gene.pdf", width = 24/2.54, height = 21/2.54)
pdf(file = "Wx80_OTU_occurrence_vs_abundance_boxplot_all_gene_min2.pdf", width = 24/2.54, height = 21/2.54)
grid.arrange(p1, p2, p3, p4, ncol = 2)
dev.off()

# Output data to file
#write.table(all_res, file = "Wx80_OTU_occurrence_freqs_sums.txt", sep = "\t", quote = F, row.names = F)
write.table(all_res, file = "Wx80_OTU_occurrence_freqs_sums_min2.txt", sep = "\t", quote = F, row.names = F)

###############################################################################
# Get the ten most abundant OTUs in each PCR, for each gene, subplot and 
# DNA extraction method
###############################################################################

# Load a taxonomy reference file, used for setting the order of taxonomy factors   
tax_ref <- read.table("Wx80_metadata/Silva_taxa_order_list_v3.txt", header = TRUE, sep = "\t", comment.char = "", quote="\"")

# Function to determine the most abundant OTUs in each set of PCRs
get_top_OTUs <- function(OTUtable, taxa, gene, n = 10){
  samples <- c("D_Psoil","D_Pmax","D_Ind","D_PO4","P_Psoil","P_Pmax","P_Ind","P_PO4")  
  top_OTUs <- data.frame()
  for(s in samples){
    top10s <- data.frame()
    PCRs <- OTUtable[, grep(s, colnames(OTUtable))]
    for(i in 1:ncol(PCRs)){ 
      j <- as.data.frame(PCRs[i], drop = FALSE)
      j[, 2] <- rownames(j)
      lab <- gsub("sample_", "", colnames(j)[1])
      names(j) <- c("Abundance","OTU")
      j <- j[order(j$Abundance, decreasing = TRUE), ] # Order by abundance
      top10 <- j[1:n, ] # Take the top 10 abundances
      top10$rank <- seq(1:n)
      top10$sample <- s
      top10$lab <- lab
      top10$subplot <- paste("Subplot", sapply(strsplit(top10$lab, "\\.|_"), "[[", 2))
      top10$method <- sapply(strsplit(top10$lab, "\\.|_"), "[[", 3)
      top10$rep <- as.numeric(sapply(strsplit(top10$lab, "\\."), "[[", 3))

      top10s <- rbind(top10s, top10)
    }
    # Fix replicate numbers (rep 01 should be 10, all others shifted by 1 - only matters for Pmax samples)      
    if(grepl("Pmax", s)){
      top10s$rep <- as.numeric(top10s$rep) - 1
      top10s$rep[top10s$rep == 0] <- 10
    # And reverse Pmax reps 1-5 for improved plot clarity?
      top10s$rep[top10s$rep %in% c(1,5) & top10s$subplot == "D"] <- rev(top10s$rep[top10s$rep %in% c(1,5) & top10s$subplot == "D"])
      top10s$rep[top10s$rep %in% c(2,4) & top10s$subplot == "D"] <- rev(top10s$rep[top10s$rep %in% c(2,4) & top10s$subplot == "D"])
      top10s$rep[top10s$rep %in% c(1,5) & top10s$subplot == "P"] <- rev(top10s$rep[top10s$rep %in% c(1,5) & top10s$subplot == "P"])
      top10s$rep[top10s$rep %in% c(2,4) & top10s$subplot == "P"] <- rev(top10s$rep[top10s$rep %in% c(2,4) & top10s$subplot == "P"])
    }
    top_OTUs <- rbind(top_OTUs, top10s)
  }
  
  # Exclude any OTUs with excessively low abundances
  top_OTUs$Abundance[top_OTUs$Abundance < 25] <- NA
  # Fix extraction method labels
  top_OTUs$method <- gsub("Psoil", "1.5 g direct", top_OTUs$method)
  top_OTUs$method <- gsub("Pmax", "7.5 g direct", top_OTUs$method)
  top_OTUs$method <- gsub("Ind", "15 g indirect", top_OTUs$method)
  top_OTUs$method <- gsub("PO4", "15 g PO4 buffer", top_OTUs$method)
  top_OTUs$method <- factor(top_OTUs$method, ordered = TRUE, 
                                levels = c("1.5 g direct", "7.5 g direct", "15 g indirect", "15 g PO4 buffer"))
  # Add taxonomic data  
  top_OTUs$OTU <- gsub("=", "_", top_OTUs$OTU)
  taxa$OTU <- rownames(taxa)
  top_OTUs_tax <- merge(top_OTUs, taxa, by = "OTU")
  top_OTUs_tax$kingdom <- factor(top_OTUs_tax$kingdom, ordered = TRUE, levels = tax_ref$node[tax_ref$node %in% top_OTUs_tax$kingdom])
  top_OTUs_tax$phylum <- factor(top_OTUs_tax$phylum, ordered = TRUE, levels = tax_ref$node[tax_ref$node %in% top_OTUs_tax$phylum])
  top_OTUs_tax$class <- factor(top_OTUs_tax$class, ordered = TRUE, levels = tax_ref$node[tax_ref$node %in% top_OTUs_tax$class])
  top_OTUs_tax$order <- factor(top_OTUs_tax$order, ordered = TRUE, levels = tax_ref$node[tax_ref$node %in% top_OTUs_tax$order])
  top_OTUs_tax <- top_OTUs_tax[order(top_OTUs_tax$subplot, top_OTUs_tax$method, top_OTUs_tax$rep), ]
  top_OTUs_tax$gene <- gene
  
  return(top_OTUs_tax)
}

# Compile top OTUs data for each gene
genes <- c("16S","18S","26S","COI")
for (gene in genes) {
  #f <- Sys.glob(paste0("Wx80_OTU_tables/*", gene, "*MEE1.0_min1_OTUtable.txt"))
  f <- Sys.glob(paste0("Wx80_OTU_tables/*", gene, "*MEE1.0_min2_OTUtable.txt"))
  OTUtable <- read.table(f, header = TRUE, row.names = 1, sep = "\t")
  rownames(OTUtable) <- gsub("=", "_", rownames(OTUtable)) # Ensure rownames match between OTU table and taxonomy
  
  # Get taxonomy data
  #tx <- Sys.glob(paste0("Wx80_taxonomy/Wx80_*", gene, "*_MEE1.0_OTUs_min1_taxonomy_v2.txt"))
  tx <- Sys.glob(paste0("Wx80_taxonomy/Wx80_*",gene,"*_MEE1.0_OTUs_min2_taxonomy_v2.txt"))
  taxa <- read.table(tx, sep="\t", header=T, row.names=1, quote = "", comment.char = "")
  
  # Limit to OTUs with prokaryote or eukaryote phylum identifications
  taxa1 <- taxa[!grepl("root|cellular organisms|No hits|Not assigned", taxa$kingdom), ]
  if(gene == "16S"){
    taxa1 <- taxa1[grepl("Archaea|Bacteria", taxa1$superkingdom) & !grepl("Archaea|Bacteria", taxa1$phylum), ]
  }
  if(gene == "18S"){
    taxa1 <- taxa1[grepl("Eukaryota", taxa1$superkingdom) & !grepl("Eukaryota|Opisthokonta|Eumetazoa|Fungi", taxa1$phylum), ]
  }
  if(gene == "26S"){
    taxa1 <- taxa1[grepl("Fungi", taxa1$kingdom) & !grepl("Fungi", taxa1$phylum), ]
  }
  if(gene == "COI"){
    taxa1 <- taxa1[grepl("Metazoa", taxa1$kingdom) & !grepl("Eumetazoa|Metazoa", taxa1$phylum), ]
  }
  OTUtable1 <- OTUtable[rownames(OTUtable) %in% rownames(taxa1), ]
  
  # Get the most abundant OTUs
  top_OTUs_tax <- get_top_OTUs(OTUtable1, taxa1, gene)

  # Make a graph of top OTUs abundance per PCR, extraction method and subplot, and output to pdf
  p <- ggplot(top_OTUs_tax, aes(x = rep, y = Abundance)) + 
    geom_line(aes(group = OTU, colour = phylum), alpha = 0.4) +
    geom_point(aes(colour = phylum), size = 2, alpha = 1, shape = 5) +
    facet_grid(subplot ~ method) +
    scale_x_continuous(breaks = seq(1:10)) + 
    scale_y_log10(breaks = c(10, 100, 1000, 10000)) + 
    guides(alpha = FALSE, size = FALSE) +
    ggtitle(paste0(gene)) +
    xlab("PCR number") + ylab("OTU abundance (sequence reads)") +
    theme(strip.background = element_blank(), panel.grid = element_blank(), plot.title = element_text(size = 9))
  p
  #ggsave(p, filename = paste0("Wx80_", gene, "_top_10_OTUs_by_target_phylum.pdf"), height = 14, width = 21, units = "cm")
  ggsave(p, filename = paste0("Wx80_", gene, "_min2_top_10_OTUs_by_target_phylum.pdf"), height = 14, width = 21, units = "cm")
}

###############################################################################
# Analysis of effects of the number of PCR numbers on biodiversity estimates
# Specifically, calculation of gamma (aka total alpha) diversity for:
# (1) Combinations of one to ten PCRs
# (2) Combinations of one to ten PCRs, subsampled to equal sequence depth per combination;
# (3) Combinations of one to ten PCRs, limited to OTUs occurring in at least x PCRs
# (4) Combinations of one to ten PCRs, limited to OTUs with abundances above x
# (5) All PCRs combined and subsampled to 10 %, 20 %, 30 % etc. to 100 % 
# This is all a bit slow to calculate
###############################################################################

# Functions for biodiversity estimation based on variable combinations of PCRs
# No subsampling before biodiversity calculations:
gamma.calcs <- function(tab){
  gamma0 <- list()
  gamma1 <- list()
  gamma2 <- list()
  table1 <- tab[,colSums(tab) > 1000] # Exclude any low abundance (i.e. failed) PCRs
  for(i in 1:ncol(table1)){
    print(paste("analysing combos", i))
    combos <- combn(table1, i, simplify = FALSE) # Get all combinations
    if(length(combos) > 25){  # Limit to 25 combinations of data (to limit processing time)
      combos <- sample(combos, 25)
    }
    combos.t <- lapply(combos, t) # transpose
    gamma0[[i]] <- sapply(combos.t, function(x) d(x, lev = "gamma", q = 0))
    gamma1[[i]] <- sapply(combos.t, function(x) d(x, lev = "gamma", q = 1))
    gamma2[[i]] <- sapply(combos.t, function(x) d(x, lev = "gamma", q = 2))
  }
  gamma.list <- list("gamma0" = gamma0, "gamma1" = gamma1, "gamma2" = gamma2)
  return(gamma.list)
}

# Subsample combined PCRs to equal sequence depth before biodiversity calculations:
gamma.subs.calcs <- function(tab){
  gamma0.subs <- list()
  gamma1.subs <- list()
  gamma2.subs <- list()
  table1 <- tab[,colSums(tab) > 1000] # Exclude any low abundance (i.e. failed) PCRs
  m <- min(colSums(table1))
  for(i in 1:ncol(table1)){
    print(paste("analysing subsampled combos", i))
    combos <- combn(table1, i, simplify = FALSE) # Get all combinations
    if(length(combos) > 25){  # Limit to 25 combinations of data (to limit processing time)
      combos <- sample(combos, 25)
    }
    combos.t <- lapply(combos, t) # transpose
    combos.r <- lapply(combos.t, function(x) as.data.frame(rrarefy(x, m/i))) # Subsample to equal total number of reads
    gamma0.subs[[i]] <- sapply(combos.r, function(x) d(x, lev = "gamma", q = 0))
    gamma1.subs[[i]] <- sapply(combos.r, function(x) d(x, lev = "gamma", q = 1))
    gamma2.subs[[i]] <- sapply(combos.r, function(x) d(x, lev = "gamma", q = 2))
  }
  gamma.subs.list <- list("gamma0_subs" = gamma0.subs, "gamma1_subs" = gamma1.subs, "gamma2_subs" = gamma2.subs)
  return(gamma.subs.list)
}

# Functions to calculate biodiversity estimates for variable numbers of PCRs,
# but limited to OTUs occurring in multiple PCRs
restr <- function(tab, thr = 5){
  occs <- apply(tab, 1, function(x) sum(x>0))
  occs <- occs[occs >= thr]
  tab <- tab[which(rownames(tab) %in% names(occs)), ]
  print(nrow(tab))
  return(tab)
}

gamma.calcs.restr <- function(tab){
  gamma0 <- list()
  gamma1 <- list()
  gamma2 <- list()
  table1 <- tab[,colSums(tab) > 1000] # Exclude any low abundance (i.e. failed) PCRs
  table1 <- restr(table1, 3) # Limit to OTUs occurring in multiple PCRs  
  for(i in 1:ncol(table1)){
    print(paste("analysing combos", i))
    combos <- combn(table1, i, simplify = FALSE) # Get all combinations
    if(length(combos) > 25){  # Limit to 25 combinations of data (to limit processing time)
      combos <- sample(combos, 25)
    }
    combos.t <- lapply(combos, t) # transpose
    gamma0[[i]] <- sapply(combos.t, function(x) d(x, lev = "gamma", q = 0))
    gamma1[[i]] <- sapply(combos.t, function(x) d(x, lev = "gamma", q = 1))
    gamma2[[i]] <- sapply(combos.t, function(x) d(x, lev = "gamma", q = 2))
  }
  gamma.list <- list("gamma0" = gamma0, "gamma1" = gamma1, "gamma2" = gamma2)
  return(gamma.list)
}

# Function to calculate biodiversity estimates for percentage subsamples of total sequences:
gamma.pc.calcs <- function(tab){
  gamma0.pc <- list()
  gamma1.pc <- list()
  gamma2.pc <- list()
  table1 <- tab[,colSums(tab) > 1000] # Exclude any low abundance (i.e. failed) PCRs
  sums <- rowSums(table1)
  m <- sum(sums) # Total sequence count
  for(i in 1:10){
    x <- i/10*m # Subsample size (steps of 10 %)
    print(paste("analysing percent subsample", x))
    sums.r <- list()
    sums.t <- t(sums)
    for(j in 1:10){
      sums.r[[j]] <- rrarefy(sums.t, x)
    }
    gamma0.pc[[i]] <- sapply(sums.r, function(x) d(x, lev = "gamma", q = 0))
    gamma1.pc[[i]] <- sapply(sums.r, function(x) d(x, lev = "gamma", q = 1))
    gamma2.pc[[i]] <- sapply(sums.r, function(x) d(x, lev = "gamma", q = 2))
  }
  gamma.pc.list <- list("gamma0_pc" = gamma0.pc, "gamma1_pc" = gamma1.pc, "gamma2_pc" = gamma2.pc) 
  return(gamma.pc.list)
}

# Function to tidy diversity estimate results
tidy.gamma <- function(gamma.res, gene, s, pc = FALSE){
  gamma.long <- melt(gamma.res)
  if(pc == FALSE){
    colnames(gamma.long) <- c("gamma value", "rep columns", "measure")
  }else{
    colnames(gamma.pc.long) <- c("gamma value", "percent x 10", "measure")
  }
  gamma.long$sample <- s 
  gamma.long$gene <- gene
  return(gamma.long)
}

# Now calculate diversity estimates for all datasets (this is rather slow):
genelist <- c("16S","18S","26S","COI")

for(gene in genelist){
  f <- Sys.glob(paste0("Wx80_OTU_tables/Wx80_*", gene, "*_MEE1.0_min2_OTUtable.txt"))
  OTUtable <- read.table(f, header = TRUE, row.names = 1, sep = "\t")
  
  # For biodiversity estimates with low abundance OTUs excluded:
  OTUtable.limited <- OTUtable[which(rowSums(OTUtable) >= 5), ]

  # Set up empty lists to hold each set of results
  gamma.data <- list()
  gamma.subs.data <- list()
  gamma.restr.data <- list()
  gamma.abund.data <- list()
  gamma.pc.data <- list()

  n <- 1
  print(paste(n, gene))
  samples <- c("D_Psoil","D_Pmax","D_Ind","D_PO4","P_Psoil","P_Pmax","P_Ind","P_PO4")
    
  for(s in samples){
    print(s)
    # Subset to set of ten PCRs 
    tab <- OTUtable[, grepl(s, colnames(OTUtable))] 
    tab.lim <- OTUtable.limited[, grepl(s, colnames(OTUtable.limited))] # Abundance limited OTUtable
    
    # Calculate biodiversity estimates
    gamma <- gamma.calcs(tab) # Estimates for multiple PCRs without subsampling
    gamma.subs <- gamma.subs.calcs(tab) # Estimates for multiple PCRs subsampled to equal read depth
    gamma.restr <- gamma.calcs.restr(tab) # Estimates limited to OTUs detected in multiple PCRs
    gamma.abund <- gamma.calcs(tab.lim) # Estimates for multiple PCRs with low abundance OTUs excluded
    gamma.pc <- gamma.pc.calcs(tab) # Estimates by percentage subsets
    
    gamma.data[[n]] <- tidy.gamma(gamma, gene, s, pc = FALSE)
    gamma.subs.data[[n]] <- tidy.gamma(gamma.subs, gene, s, pc = FALSE)
    gamma.restr.data[[n]] <- tidy.gamma(gamma.restr, gene, s, pc = FALSE)
    gamma.abund.data[[n]] <- tidy.gamma(gamma.abund, gene, s, pc = FALSE)
    gamma.pc.data[[n]] <- tidy.gamma(gamma.pc, gene, s, pc = TRUE)
    
    n <- n +1
  }
  
  # Organise the data into a table, then output to file
  gamma.div <- melt(gamma.data, id.vars = c("gamma value", "rep columns", 
                                            "measure", "gene", "sample"))
  gamma.div.subs <- melt(gamma.subs.data, id.vars = c("gamma value", "rep columns", 
                                                      "measure", "gene", "sample"))
  gamma.div.restr <- melt(gamma.restr.data, id.vars = c("gamma value", "rep columns", 
                                                        "measure", "gene", "sample"))
  gamma.div.abund <- melt(gamma.abund.data, id.vars = c("gamma value", "rep columns", 
                                                        "measure", "gene", "sample"))
  gamma.div.pc <- melt(gamma.pc.data, id.vars = c("gamma value", "percent x 10", 
                                                  "measure", "gene", "sample"))
  
  write.table(gamma.div, file = paste0("Wx80_", gene, "_PCR_reps_min2_gamma_div.txt"), 
               sep = "\t", quote = F, row.names = F)
  write.table(gamma.div.subs, file = paste0("Wx80_", gene, "_PCR_reps_min2_gamma_div_subsampled.txt"), 
               sep = "\t", quote = F, row.names = F)
  write.table(gamma.div.restr, file = paste0("Wx80_", gene, "_PCR_reps_min2_gamma_div_restr3.txt"), 
              sep = "\t", quote = F, row.names = F)
  write.table(gamma.div.abund, file = paste0("Wx80_", gene, "_PCR_reps_min2_gamma_div_thr5.txt"), 
              sep = "\t", quote = F, row.names = F)
  write.table(gamma.div.pc, file = paste0("Wx80_", gene, "_PCR_reps_min2_gamma_div_percent.txt"), 
               sep = "\t", quote = F, row.names = F)
}

# Now plot alpha diversity estimate curves
# Load all the data
f1 <- Sys.glob("Wx80_alpha_diversity_per_reps/*replicates_min2_gamma_div.txt")
f2 <- Sys.glob("Wx80_alpha_diversity_per_reps/*replicates_min2_gamma_div_subsampled.txt") # Subsampled combinations data
f3 <- Sys.glob("Wx80_alpha_diversity_per_reps/*replicates_min2_gamma_div_restr3.txt") # Restricted data
f4 <- Sys.glob("Wx80_alpha_diversity_per_reps/*replicates_min2_gamma_div_thr5.txt") # Abundance limited data
f5 <- Sys.glob("Wx80_alpha_diversity_per_reps/*replicates_min2_gamma_div_percent.txt") # Percent subsampled data

# Load all the data
n <- 1
all.div <- data.frame()
for (f in f1) {
  df <- read.table(f, sep = "\t", header = TRUE)
  gene = sapply(strsplit(f, "per_reps/Wx80_"), "[[", 2)
  gene = sapply(strsplit(gene, "_PCR"), "[[", 1)
  df$gene <- gene
  df$var <- "combos"
  all.div <- rbind(all.div, df)
} 

for (f in f2) {  
  df <- read.table(f, sep = "\t", header = TRUE)
  gene = sapply(strsplit(f, "per_reps/Wx80_"), "[[", 2)
  gene = sapply(strsplit(gene, "_PCR"), "[[", 1)
  df$gene <- gene
  df$var <- "combos_subsampled"
  all.div <- rbind(all.div, df)
}

for (f in f3) {  
  df <- read.table(f, sep = "\t", header = TRUE)
  gene = sapply(strsplit(f, "per_reps/Wx80_"), "[[", 2)
  gene = sapply(strsplit(gene, "_PCR"), "[[", 1)
  df$gene <- gene
  df$var <- "combos_restricted"
  all.div <- rbind(all.div, df)
}

for (f in f4) {  
  df <- read.table(f, sep = "\t", header = TRUE)
  gene = sapply(strsplit(f, "per_reps/Wx80_"), "[[", 2)
  gene = sapply(strsplit(gene, "_PCR"), "[[", 1)
  df$gene <- gene
  df$var <- "abund_limited"
  all.div <- rbind(all.div, df)
}

for (f in f5) {  
  df <- read.table(f, sep = "\t", header = TRUE)
  gene = sapply(strsplit(f, "per_reps/Wx80_"), "[[", 2)
  gene = sapply(strsplit(gene, "_PCR"), "[[", 1)
  df$gene <- gene
  df$var <- "percent_subsamples"
  colnames(df) <- gsub("percent.x.10", "rep.columns", colnames(df)) # make column match the other datasets
  all.div <- rbind(all.div, df)
}

# Organise the data for plotting
all.div$subplot <- sapply(strsplit(as.character(all.div$sample), "_"), "[[", 1)
all.div$method <- sapply(strsplit(as.character(all.div$sample), "_"), "[[", 2)
all.div$method <- factor(all.div$method, ordered = TRUE, 
                         levels = c("Psoil", "Pmax", "Ind", "PO4"))

all.div$method2 <- gsub("Psoil", "10 x 1.5 g direct", all.div$method)
all.div$method2 <- gsub("Pmax", "2 x 7.5 g direct", all.div$method2)
all.div$method2 <- gsub("Ind", "1 x 15 g indirect", all.div$method2)
all.div$method2 <- gsub("PO4", "1 x 15 g PO4 buffer", all.div$method2)
all.div$method2 <- factor(all.div$method2, ordered = TRUE, 
                          levels = c("10 x 1.5 g direct", "2 x 7.5 g direct", "1 x 15 g indirect", "1 x 15 g PO4 buffer"))

all.div$measure2 <- gsub("gamma|_subs|_pc", "", all.div$measure)
all.div$method2.var <- paste(all.div$method2)

# Plot alpha 0/1 curves for combinations of PCRs and subsampled combinations of PCRs, for each subplot separately)
for(s in c("D", "P")){
  div1 <- all.div[all.div$subplot == s & all.div$measure2 != "2" & all.div$var == "combos", ]
  p1 <- ggplot(div1, aes(x = rep.columns, y = gamma.value, colour = measure2)) + 
    geom_jitter(alpha = 0.5) +
    scale_colour_manual(values = c("#E69F00", "#56B4E9")) +
    scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) +
    xlab("Number of pooled PCR replicates per subplot/method") + ylab("alpha diversity") +
    facet_grid(gene ~ method2, scales = "free") +
    theme(strip.background = element_blank(), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), legend.position = "None")
  p1
  
  div2 <- all.div[all.div$subplot == s & all.div$measure2 != "2" & all.div$var == "combos_subsampled", ]
  p2 <- ggplot(div2, aes(x = rep.columns, y = gamma.value, colour = measure2)) + 
    geom_jitter(alpha = 0.5) +
    scale_colour_manual(values = c("#E69F00", "#56B4E9")) +
    scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) +
    xlab("Number of pooled PCR replicates per subplot/method (subsampled)") + ylab(" ") +
    facet_grid(gene ~ method2, scales = "free") +
    theme(strip.background = element_blank(), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), legend.position = "None")
  p2
  pdf(file = paste0("Wx80_all_genes_PCR_reps_alpha-0-1_subplot", s, ".pdf"), height = 21/2.54, width = 29.7/2.54, useDingbats = FALSE)
  grid.arrange(p1, p2, ncol = 2)
  dev.off()
}

# Plot alpha 0/1 curves for data limited to OTUs with at least 5 sequences, for both subplots together
div <- all.div[all.div$measure2 != "2" & all.div$var == "combos", ]
div$method2.subplot <- paste0(div$subplot, ", ", div$method2)
div$method2.subplot <- factor(div$method2.subplot, levels = c("D, 10 x 1.5 g direct", "D, 2 x 7.5 g direct", "D, 1 x 15 g indirect",
                                                              "D, 1 x 15 g PO4 buffer", "P, 10 x 1.5 g direct", "P, 2 x 7.5 g direct",
                                                              "P, 1 x 15 g indirect","P, 1 x 15 g PO4 buffer"), ordered = TRUE)
p <- ggplot(div, aes(x = rep.columns, y = gamma.value, colour = measure2)) + 
  geom_jitter(alpha = 0.5) +
  scale_colour_manual(values = c("#E69F00", "#56B4E9")) +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) +
  xlab("Number of pooled PCR replicates per subplot/method (restricted to OTUs with at least five sequences)") + ylab("alpha diversity") +
  facet_grid(gene ~ method2_subplot, scales = "free") +
  theme(strip.background = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), legend.position = "None")
p
ggsave(p, "Wx80_all_genes_PCR_reps_alpha-0-1_thr5_subplots_D_P.pdf", height = 21, width = 29.7, units = "cm", useDingbats = FALSE)

### Plot alpha 0/1 curves for all data restricted to OTUs in multiple PCRs, for both subplots together
div <- all.div[all.div$measure2 != "2" & all.div$var == "combos_restricted", ]
div$method2.subplot <- paste0(div$subplot, ", ", div$method2)
div$method2.subplot <- factor(div$method2.subplot, levels = c("D, 10 x 1.5 g direct", "D, 2 x 7.5 g direct", "D, 1 x 15 g indirect",
                                                              "D, 1 x 15 g PO4 buffer", "P, 10 x 1.5 g direct", "P, 2 x 7.5 g direct",
                                                              "P, 1 x 15 g indirect","P, 1 x 15 g PO4 buffer"), ordered = TRUE)
p <- ggplot(div, aes(x = rep.columns, y = gamma.value, colour = measure2)) + 
  geom_jitter(alpha = 0.5) +
  scale_colour_manual(values = c("#E69F00", "#56B4E9")) +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) +
  xlab("Number of pooled PCR replicates per subplot/method (restricted to OTUs in at least 3 PCRs)") + ylab("alpha diversity") +
  facet_grid(gene ~ method2_subplot, scales = "free") +
  theme(strip.background = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), legend.position = "None")
p

ggsave(p, "Wx80_all_genes_PCR_reps_restr3_alpha-0-1.pdf", height = 21, width = 29.7, units = "cm", useDingbats = FALSE)

### Plot alpha 0/1 curves for all data subsampled by percentages of sequences, for both subplots together
div <- all.div[all.div$measure2 != "2" & all.div$var == "percent_subsamples", ]
div$method2.subplot <- paste0(div$subplot, ", ", div$method2)
div$rep.columns <- div$rep.columns * 10
div$method2.subplot <- factor(div$method2.subplot, levels = c("D, 10 x 1.5 g direct", "D, 2 x 7.5 g direct", "D, 1 x 15 g indirect",
                                                              "D, 1 x 15 g PO4 buffer", "P, 10 x 1.5 g direct", "P, 2 x 7.5 g direct",
                                                              "P, 1 x 15 g indirect","P, 1 x 15 g PO4 buffer"), ordered = TRUE)
p <- ggplot(div, aes(x = rep.columns, y = gamma.value, colour = measure2)) + 
  geom_jitter(alpha = 0.5) +
  scale_colour_manual(values = c("#E69F00", "#56B4E9")) +
  scale_x_continuous(breaks = c(20, 40, 60, 80, 100)) +
  xlab("Percentage of total sequence reads per subplot/method") + ylab("alpha diversity") +
  facet_grid(gene ~ method2_subplot, scales = "free") +
  theme(strip.background = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), legend.position = "None")
p
ggsave(p, "Wx80_all_genes_PCR_reps_percent_subs_alpha-0-1.pdf", height = 21, width = 29.7, units = "cm", useDingbats = FALSE)
