###############################################################################
# Comparisons of biodiversity composition and community structure among the 
# results of different DNA extraction methods
###############################################################################

library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(plyr)
library(RColorBrewer)
library(reshape2)
library(scales)
library(vegan)
library(vegetarian)

# Set ggplot theme
theme_set(theme_bw(base_size=8))

# Set up a colour palette
library(RColorBrewer)
colors <- brewer.pal(8, "Spectral")
pal <- colorRampPalette(colors) 

# Set working directory
setwd("./")

###############################################################################
# Comparison of biodiversity estimates by phylum among 
# different samples and DNA extraction methods
###############################################################################

# Compile biodiversity estimates for each gene, sample, DNA extraction method, and phylum
get.diversity.estimates <- function(s.tax, tax, i, s, m){
  tax.subset <- s.tax[grep(tax, s.tax$taxon), ]
  tax.subset$taxon <- NULL
  tax.subset$Row.names <- NULL
  a0 <- d(t(tax.subset), lev = "alpha", q = 0)
  a1 <- d(t(tax.subset), lev = "alpha", q = 1)
  b0 <- d(t(tax.subset), lev = "beta", q = 0)
  b1 <- d(t(tax.subset), lev = "beta", q = 1)
  g0 <- d(t(tax.subset), lev = "gamma", q = 0)
  g1 <- d(t(tax.subset), lev = "gamma", q = 1)
  div <- list("rep" = paste0("rep.", i), "sample" = s, "depth" = m, "taxon" = tax,
              "a0" = a0, "a1" = a1, "b0" = b0, "b1" = b1, "g0" = g0, "g1" = g1)
  return(div)
}

genes <- c("16S","18S","26S","COI")
for (gene in genes) {
  
  #f <- Sys.glob(paste0("Wx80_OTU_tables/Wx80_*", gene, "*_MEE1.0_min1_OTUtable.txt"))
  f <- Sys.glob(paste0("Wx80_OTU_tables/Wx80_*", gene, "*_MEE1.0_min2_OTUtable.txt"))
  OTUtable <- read.table(f, header = TRUE, row.names = 1, sep = "\t")
  rownames(OTUtable) <- gsub("=", "_", rownames(OTUtable)) # Ensure rownames match between OTU table and taxonomy
  print(paste("analysing", gene, "..."))
  
  # Get taxonomy data
  #tx <- Sys.glob(paste0("Wx80_taxonomy/Wx80_*", gene, "*_MEE1.0_OTUs_min1_taxonomy_v2.txt"))
  tx <- Sys.glob(paste0("Wx80_taxonomy/Wx80_*",gene,"*_MEE1.0_OTUs_min2_taxonomy_v2.txt"))
  taxa <- read.table(tx, sep="\t", header=T, row.names=1, quote = "", comment.char = "")
  
  #taxa <- taxa[,c("kingdom","phylum","class","order")]
  taxa <- taxa[,c("phylum"), drop = FALSE]
  colnames(taxa) <- "taxon"
  
  # Exclude excessively low abundance samples
  x <- mean(colSums(OTUtable))
  OTUtable.s <- as.data.frame(OTUtable[, colSums(OTUtable != 0) > 50])  # Exclude any samples with fewer than x OTUs
  OTUtable.s <- as.data.frame(OTUtable.s[, colSums(OTUtable.s) > x*0.25])  # Exclude any samples with excessively low abundance
  OTUtable.s <- as.data.frame(OTUtable.s[rowSums(OTUtable.s) > 0, ]) # Exclude any empty rows    
  print(paste("columns dropped:", ncol(OTUtable)-ncol(OTUtable.s)))
  
  # Get minimum for subsampling
  m <- min(colSums(OTUtable.s))
  print(paste("minimum colsums:", m))
  
  # Split table by subplot and extraction method
  samples <- c("D_Psoil","D_Pmax","D_Ind","D_PO4","P_Psoil","P_Pmax","P_Ind","P_PO4")  
  
  data.s = list()
  data.ns = list()
  
  # calculation of subsampled biodiversity estimates by phylum
  n <- 1
  for(s in samples){
    tab <- OTUtable.s[, grepl(s, colnames(OTUtable.s))]
    j <- 1
    dlist = list()
    for(i in 1:10){
      print(paste("Rep:", i))
      print(paste("Subsampling to", m))
      s.rare <- t(rrarefy(t(tab), m))
      s.tax <- merge(s.rare, taxa, by = "row.names")
      s.tax$taxon <- as.factor(s.tax$taxon)
      tax.levs <- levels(s.tax$taxon)   
      for(tax in tax.levs){
        print(paste("calculating diversity for", s, tax, "rep", i))
        div <- get.diversity.estimates(s.tax, tax, i, s, m)
        dlist[[j]] <- div
        j <- j + 1
      }
    }
    data.s[[n]] <- rbindlist(dlist, use.names = TRUE)
    n <- n + 1
  }
  
  # Calculation of non-subsampled biodiversity estimates by phylum
  n <- 1
  for(s in samples){
    j <- 1
    dlist = list()
    print(s)
    s.tax <- merge(s, taxa, by = "row.names")
    s.tax$taxon <- as.factor(s.tax$taxon)
    tax.levs <- levels(s.tax$taxon)   
    for(tax in tax.levs){
      print(paste("calculating diversity for", s, tax, "(unsubsampled)"))
      div <- get.diversity.estimates(s.tax, tax, i = "x", s, m = "x")
      dlist[[j]] <- div
      j <- j + 1
      
    }
    data.ns[[n]] <- rbindlist(dlist, use.names = TRUE)
    n <- n + 1
  }  
  
  # Output results to file
  data1 <- rbindlist(data.s, use.names = TRUE)
  write.table(data1, file = paste0("Wx80_", gene, "_effective_div_min2_by_method-phylum_subsampled.txt"), sep = "\t", 
              quote = FALSE, row.names = FALSE)
  data2 <- rbindlist(data.ns, use.names = TRUE)
  write.table(data2, file = paste0("Wx80_", gene, "_effective_div_min2_by_method-phylum_NOT-subsampled.txt"), sep = "\t", 
              quote = FALSE, row.names = FALSE)
}


# Now make plots of biodiversity estimates by DNA extraction method and phylum
# Load taxonomy reference file, used for setting order of taxonomic factors
tax_ref <- read.table("Wx80_metadata/Silva_taxa_order_list_v3.txt", sep = "\t", quote = "", comment.char = "", header = TRUE) 

genelist <- c("16S", "18S", "26S", "ShCO1")

for(g in genelist){
  
  f1 <- Sys.glob(paste0("Wx80_", g, "_effective_div_min2_by_method-phylum_subsampled.txt"))
  #f2 <- Sys.glob(paste0("Wx80_effective_div_analysis_2016/Wx80_", g, "_effective-a-b-g_div_min2_by_subplot-method-phylum_NOT-subsampled.txt"))

  df <- read.table(f1, header = TRUE, sep = "\t", check.names = FALSE)
  df <- melt(df, id.vars = c("rep", "sample", "depth", "taxon"))  

  # Get taxonomy data
  #tx <- Sys.glob(paste0("Wx80_taxonomy/Wx80_*", gene, "*_MEE1.0_OTUs_min1_taxonomy_v2.txt"))
  tx <- Sys.glob(paste0("Wx80_taxonomy/Wx80_*",gene,"*_MEE1.0_OTUs_min2_taxonomy_v2.txt"))
  taxa <- read.table(tx, sep="\t", header=T, row.names=1, quote = "", comment.char = "")
  df$kingdom <- taxa$kingdom[match(df$taxon, taxa$phylum)]
  df$phylum <- taxa$phylum[match(df$taxon, taxa$phylum)]
  
  ### Organise factor levels
  df$method <- paste(df$sample)
  df$method <- gsub("(D|P)_", "", df$method)
  df$method <- gsub("Psoil", "1.5 g direct", df$method)
  df$method <- gsub("Pmax", "7.5 g direct", df$method)
  df$method <- gsub("Ind", "15 g indirect", df$method)
  df$method <- gsub("PO4", "15 g PO4 buffer", df$method)
  df$method <- factor(df$method, levels = c("1.5 g direct", "7.5 g direct", "15 g indirect", 
                                                  "15 g PO4 buffer"), ordered = TRUE)
  df$taxon <- factor(df$taxon, levels = tax_ref$node, ordered = TRUE)
  df$variable <- gsub("a", "alpha ", df$variable)
  df$variable <- gsub("b", "beta ", df$variable)
  df$variable <- gsub("g", "gamma ", df$variable)
  
  df$var2 <- gsub("a", "alpha ", df$variable)
  df$var2 <- gsub("b", "beta ", df$variable)
  df$var2 <- gsub("g", "gamma ", df$variable)
  df$subplot <- sapply(strsplit(as.character(df$sample), "_"), "[[", 1)
  
  # Make a data summary?
  #outfile = "Wx80_effective_div_by_subplot-method-phylum_min2_subsampled_means.txt"
  #dt1 <- data.table(df)
  #means <- dt1[, .(min = min(value), max = max(value), mean = mean(value)), by=.(variable, sample, subplot, taxon, method)]
  #means$gene <- g
  #write.table(means, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE, append = TRUE) 
  
  # Drop non-phyla groups (i.e. those where phylum matches a kingdom), and other questionable groups
  df <- df[!(df$taxon %in% df$kingdom), ]
  df <- df[!(grepl("incertae sedis|Eumetazoa|Dikarya|Nucleariidae|Cnidaria|Echinodermata|Porifera|Brachiopoda|unclassified", 
                         df$taxon)), ]
  
  # And limit to prokaryotes (for 16S) or eukaryotes (for 18S, 26S, COI)
  if(g == "16S") {
    df.sub <- df[grepl("Archaea|Bacteria|Prokaryota", df$kingdom), ]
    df.sub$taxon <- gsub(" <phylum>", "", df.sub$taxon)
  } else {
    df.sub <- df[!(grepl("Archaea|Bacteria|Prokaryota", df$kingdom)), ]
  }
  
  # Drop low abundance groups (mean richness < 4)
  dt <- data.table(df.sub)
  means <- dt[, .(mean = mean(value)), by=.(variable, taxon)]
  exclude <- means[means$variable == "alpha 0" & means$mean < 4, ]
  df.sub <- df.sub[!(df.sub$taxon %in% exclude$taxon), ]
  
  # Drop alpha2/gamma2
  df.sub <- df.sub[!(grepl(" 2", df.sub$variable)), ]
  
  # Make a graph
  p <- ggplot(df.sub, aes(x = method, y = value)) + 
    geom_jitter(aes(colour = method, shape = subplot), width = 0.4, height = 0, size = 3, alpha = 0.75) + 
    scale_shape_manual(values = c(47,92)) + # slashes
    #scale_shape_manual(values = c(1,5)) + # circle and diamond
    facet_grid(variable ~ taxon, scales = "free") + 
    scale_colour_brewer(palette = "Set1") +
    ylab("Diversity estimate") + xlab("") +
    theme(strip.background = element_blank(), panel.border = element_rect(colour = "black"),
          legend.position = "none", axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))
  p
  # Output graph to file
  ggsave(p, filename = paste0("Wx80_", g, "_eff_div_by_subplot-phylum_min2_0-1_revised_slashes.pdf"), 
         width = 21, height = 18, units = "cm", useDingbats = FALSE)
}

###############################################################################
# Calculate the overlap of OTUs among DNA extraction methods
###############################################################################

# Compile OTU abundance and occurrence data per gene and method
genes <- c("16S","18S","26S","COI")
all.res <- data.frame()
for (gene in genes) {
  #f <- Sys.glob(paste0("Wx80_OTU_tables/Wx80_*", gene, "*_MEE1.0_min1_OTUtable.txt"))
  f <- Sys.glob(paste0("Wx80_OTU_tables/Wx80_*", gene, "*_MEE1.0_min2_OTUtable.txt"))
  OTUtable <- read.table(f, header = TRUE, row.names = 1, sep = "\t")
  rownames(OTUtable) <- gsub("=", "_", rownames(OTUtable)) # Ensure rownames match between OTU table and taxonomy
  print(paste("analysing", gene, "..."))
  
  # Get taxonomy data
  #tx <- Sys.glob(paste0("Wx80_taxonomy/Wx80_*", gene, "*_MEE1.0_OTUs_min1_taxonomy_v2.txt"))
  tx <- Sys.glob(paste0("Wx80_taxonomy/Wx80_*",gene,"*_MEE1.0_OTUs_min2_taxonomy_v2.txt"))
  taxa <- read.table(tx, sep="\t", header=T, row.names=1, quote = "", comment.char = "")
  taxa <- taxa[,c("kingdom","phylum","class","order")]
  
  # Calculate OTU abundances and occurrences among PCRs
  res <-data.frame(matrix(nrow = nrow(OTUtable), ncol = 0))
  for(m in c("Psoil","Pmax","Ind","PO4")){
    x <- OTUtable[, grep(m, colnames(OTUtable))]
    means <- rowMeans(x) # Mean abundance of OTUs per method
    sums <- rowSums(x) # Total abundance of OTUs per method
    occs <- apply(x, 1, function(x) sum(x > 0, na.rm = TRUE)) # OTU occurrence among PCRs per method
    y <- cbind(means, sums, occs)
    colnames(y) <- paste0(m, "_", colnames(y))
    res <- cbind(res, y)
  }
  res$gene <- gene
  res$total <- rowSums(res[, grepl("sums", colnames(res))])
  
  # Fix taxonomic group for protists
  taxa$group <- as.character(taxa$kingdom)
  taxa$group[grepl("Alveolata|Rhizaria|Stramenopiles", taxa$kingdom)]  <- "Protists"
  taxa$group[grepl("Eukaryota", taxa$group) & !grepl("Eukaryota", taxa$phylum)]  <- "Protists"
  taxa$group[grepl("Opisthokonta", taxa$kingdom) & grepl("Choanoflagellida", taxa$phylum)] <- "Protists"
  taxa$group[grepl("Opisthokonta", taxa$phylum)]  <- "Eukaryota"
  taxa$group[grepl("root|cellular organisms|Not assigned|No hits", taxa$kingdom)] <- "unknown"
  taxa$group <- factor(taxa$group, levels = c("Archaea","Bacteria","Eukaryota","Protists","Fungi","Metazoa","Viridiplantae","unknown"),
                       ordered = TRUE)
  # Combine data and taxa
  rownames(res) <- lapply(rownames(res), FUN = function(x) gsub("=", "_", x)) # Make OTU ids match between OTUtable and taxonomy
  res$group <- taxa$group[match(rownames(res), rownames(taxa))]
  res$group[is.na(res$group)] <- "unknown"
  head(res)
  all.res <- rbind(res, all.res)
}

summary(all.res)

# Now calculate overlap of OTUs between methods
get.overlap <- function(all.res , v = "sums", total.thr = 1, thr = 1){
  dt <- data.table(all.res[all.res$total >= total.thr, grepl(v, colnames(all.res))|grepl("gene|group", colnames(all.res))])
  colnames(dt) <- gsub(paste0("_", v), "", colnames(dt))
  summary(dt)
  
  # Determine OTU overlap between methods; summarise by "length" to get presence/absence
  a <- dt[(Psoil >= thr & Pmax < thr & Ind < thr & PO4 < thr), .(OTUs = length(Psoil)), by = c("gene", "group")]
  b <- dt[(Psoil < thr & Pmax >= thr & Ind < thr & PO4 < thr), .(OTUs = length(Psoil)), by = c("gene", "group")]
  c <- dt[(Psoil < thr & Pmax < thr & Ind >= thr & PO4 < thr), .(OTUs = length(Psoil)), by = c("gene", "group")]
  d <- dt[(Psoil < thr & Pmax < thr & Ind < thr & PO4 >= thr), .(OTUs = length(Psoil)), by = c("gene", "group")]
  a <- dt[(Psoil >= thr & Pmax < thr & Ind < thr & PO4 < thr), .(OTUs = length(Psoil)), by = c("gene", "group")]
  b <- dt[(Psoil < thr & Pmax >= thr & Ind < thr & PO4 < thr), .(OTUs = length(Psoil)), by = c("gene", "group")]
  c <- dt[(Psoil < thr & Pmax < thr & Ind >= thr & PO4 < thr), .(OTUs = length(Psoil)), by = c("gene", "group")]
  d <- dt[(Psoil < thr & Pmax < thr & Ind < thr & PO4 >= thr), .(OTUs = length(Psoil)), by = c("gene", "group")]
  ab <- dt[(Psoil >= thr & Pmax >= thr & Ind < thr & PO4 < thr), .(OTUs = length(Psoil)), by = c("gene", "group")]
  ac <- dt[(Psoil >= thr & Pmax < thr & Ind >= thr & PO4 < thr), .(OTUs = length(Psoil)), by = c("gene", "group")]
  ad <- dt[(Psoil >= thr & Pmax < thr & Ind < thr & PO4 >= thr), .(OTUs = length(Psoil)), by = c("gene", "group")]
  bc <- dt[(Psoil < thr & Pmax >= thr & Ind >= thr & PO4 < thr), .(OTUs = length(Psoil)), by = c("gene", "group")]
  bd <- dt[(Psoil < thr & Pmax >= thr & Ind < thr & PO4 >= thr), .(OTUs = length(Psoil)), by = c("gene", "group")]
  cd <- dt[(Psoil < thr & Pmax < thr & Ind >= thr & PO4 >= thr), .(OTUs = length(Psoil)), by = c("gene", "group")]
  abc <- dt[(Psoil >= thr & Pmax >= thr & Ind >= thr & PO4 < thr), .(OTUs = length(Psoil)), by = c("gene", "group")]
  abd <- dt[(Psoil >= thr & Pmax >= thr & Ind < thr & PO4 >= thr), .(OTUs = length(Psoil)), by = c("gene", "group")]
  acd <- dt[(Psoil >= thr & Pmax < thr & Ind >= thr & PO4 >= thr), .(OTUs = length(Psoil)), by = c("gene", "group")]
  bcd <- dt[(Psoil < thr & Pmax >= thr & Ind >= thr & PO4 >= thr), .(OTUs = length(Psoil)), by = c("gene", "group")]
  abcd <- dt[(Psoil >= thr & Pmax >= thr & Ind >= thr & PO4 >= thr), .(OTUs = length(Psoil)), by = c("gene", "group")]
  #total <- dt[,.(OTUs = length(Psoil)), by = c("gene", "group")]
  
  # Give each overlap category a simple label
  a$sample <- "A only" #"Ps_only"
  b$sample <- "B only" #"Pm_only"
  c$sample <- "C only" #"Ind_only"
  d$sample <- "D only" #"PO4_only"
  ab$sample <- "AB--" #"Ps_Pm"
  ac$sample <- "A-C-" #"Ps_Ind"
  ad$sample <- "A--D" #"Ps_PO4"
  bc$sample <- "-BC-" #"Pm_Ind"
  bd$sample <- "-B-D" #"Pm_PO4"
  cd$sample <- "--CD" #"Ind_PO4"
  abc$sample <- "ABC-" #"Ps_Pm_Ind"
  abd$sample <- "AB-D" #"Ps_Pm_PO4"
  acd$sample <- "A-CD" #"Ps_Ind_PO4"
  bcd$sample <- "-BCD" #"Pm_Ind_PO4"
  abcd$sample <- "All" #"All_methods"
  #total$sample <- "Total"
  
  overlap <- rbind(a, b, c, d, ab, ac, ad, bc, bd, cd, abc, abd, acd, bcd, abcd)
  overlap$total.thr <- total.thr
  overlap$method.thr <- thr
  overlap$proportion <- 0
  for(g in unique(overlap$gene)){ 
    overlap$proportion[overlap$gene == g] <- overlap$OTUs[overlap$gene == g]/(sum(overlap$OTUs[overlap$gene == g]))*100
  }
  return(overlap)
}

overlap.res <- list()
for(t in seq(1:10)){
  ov <- get.overlap(all.res = all.res, v = "sums", total.thr = t, thr = 1)
  overlap.res[[t]] <- ov
}

overlap.res <- do.call("rbind", overlap.res)
summary(overlap.res) 
unique(overlap.res$sample)
overlap.res$sample <- factor(overlap.res$sample, levels = rev(c("A only","B only","C only","D only", "AB--", "A-C-", "A--D",
                                                                "-BC-","-B-D","--CD","ABC-","AB-D","A-CD","-BCD","All")), ordered = TRUE)

# Output data to file
#write.table(overlap_res, file = "Wx80_OTUs_overlap_by_method_revised.txt", sep = "\t", quote = F, row.names = F)
write.table(overlap.res, file = "Wx80_OTUs_overlap_by_method_revised_min2.txt", sep = "\t", quote = F, row.names = F)

# Make a graph of overlap between methods by gene and abundance threshold
# For min2, threshold 1 and 2 results are the same
ggplot(overlap.res, aes(y = proportion, x = sample, fill = group, order = group)) + 
  geom_bar(stat = "identity") + scale_fill_brewer(palette = "Spectral") + 
  ylab("Proportion of OTUs") + xlab("DNA extraction methods") +
  coord_flip() + scale_y_continuous(breaks = pretty_breaks(3)) + 
  facet_grid(gene ~ total.thr, scales = "free") + 
  theme(panel.grid = element_blank(), strip.background = element_blank()) +
  guides(fill = guide_legend(title = ""))

#ggsave("Wx80_OTU_overlap_btw_methods_thr_1-10.pdf", width = 25, height = 15, units = "cm")
ggsave("Wx80_OTU_overlap_btw_methods_thr_1-10_min2.pdf", width = 25, height = 15, units = "cm")

# Only graph a subset of thresholds
ggplot(overlap.res[overlap.res$total.thr %in% c(2,3,5,10), ], 
       aes(y = proportion, x = sample, fill = group, order = group)) + 
  geom_bar(stat = "identity") + scale_fill_brewer(palette = "Spectral") + 
  ylab("Proportion of OTUs") + xlab("DNA extraction methods") +
  coord_flip() + scale_y_continuous(breaks = pretty_breaks(4)) + 
  facet_grid(gene ~ total.thr, scales = "free") + 
  theme(panel.grid = element_blank(), strip.background = element_blank()) +
  guides(fill = guide_legend(title = ""))

#ggsave("Wx80_OTU_overlap_btw_methods_thr_1-3_5_10.pdf", width = 20, height = 15, units = "cm")
ggsave("Wx80_OTU_overlap_btw_methods_thr_2-3_5_10_min2.pdf", width = 18, height = 15, units = "cm")

###############################################################################
# Multivariate comparisons of community structure between samples and 
# DNA extraction methods (MDS plots, dispersion and PERMANOVA tests)
###############################################################################

source("../R stuff/beta1_function.R")

# Function to generate an MDS ordination plot
MDS.plot <- function(df, n, metric = "jaccard", bin = TRUE, metadata){
  df.t <- t(df)
  
  if(metric == "beta1"){
    dist <- beta1(df.t, 1)
    attributes(dist)$Labels <- rownames(df.t)
  }else{
    dist <- vegdist(df.t, method = metric, binary = bin)
  }
  mds <- metaMDS(dist)
  
  # Extract data for ggplot
  pts <- as.data.frame(mds$points)
  str <- mds$stress

  # Combine with metadata
  rownames(pts) <- gsub("sample_", "", rownames(pts))
  pts <- merge(pts, metadata, by = "row.names")
  hulls <- ddply(pts, "Sample", find_hull)
  
  p <- ggplot(pts, aes(x = MDS1, y = MDS2)) + 
    geom_point(aes(colour = Subplot, shape = Method), size = 3, alpha = 0.7) + 
    scale_shape(solid = FALSE) +
    #geom_text(aes(label = pts$Rep), size = 4, vjust = 1.5) +
    geom_polygon(data = hulls, aes(mapping = hulls$Sample, colour = hulls$Subplot), fill = NA) +
    scale_colour_brewer(palette = "Set1") +
    theme(panel.grid = element_blank(), plot.title = element_text(size = 8),
          plot.margin = unit(c(0.1,0.25,0.25,0), "cm")) +
    xlab("") + ylab("") +
    ggtitle(paste0(letters[[n]], ". ", g,", ", taxon, ", ", metric, " distance (stress:", round(str, 3), ")"))
  
  return(p)
}

# Function to retrieve legend from plot
get.legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Function to find perimeter of points in MDS ordination plot
find.hull <- function(df) df[chull(df$MDS1, df$MDS2), ]

# Run beta-dispersion and PERMANOVA tests
bdisp.PERM.test <- function(df, metric = "jaccard", bin = TRUE, metadata = metadata){
  df.t <- t(df)
  if(metric == "beta1"){
    dist <- beta1(df.t, 1)
    attributes(dist)$Labels <- rownames(df.t)
  }else{
    dist <- vegdist(df.t, method = metric, binary = bin)
  }
  
  subplots <- as.factor(sapply(strsplit(rownames(df.t), "_"), "[[", 1))
  samples <- as.factor(sapply(strsplit(rownames(df.t), "-"), "[[", 2))
  methods <- as.factor(sapply(strsplit(rownames(df.t), "_|-"), "[[", 3))
  
  # Test for beta-dispersion differences between subplots   
  s.disp <- betadisper(dist, subplots)
  s.disp
  anova(s.disp)
  permutest(s.disp, pairwise = TRUE)
  plot(s.disp, main = paste(g, taxon, metric, ", subplot"))
  boxplot(s.disp, main = paste(g, taxon, metric, ", subplot"), cex.axis = 0.5, cex.main = 1)
  
  # Test for beta-dispersion differences between samples/methods 
  m.disp <- betadisper(dist, samples)
  m.disp
  anova(m.disp)
  permutest(m.disp, pairwise = TRUE)
  plot(m.disp, main = paste(g, taxon, metric, "sample/method"))
  boxplot(m.disp, main = paste(g, taxon, metric, "sample/method"), cex.axis = 0.5, cex.main = 1)
  
  # Permanova test
  #adonis(dist ~ subplots * methods, strata = subplots)
  print(adonis(dist ~ subplots * methods))
}

# Get within- and between-subplot multivariate distances
within.between <- function(df, metric = "jaccard", bin = TRUE){
  df.t <- t(df)
  if(metric == "1-beta1"){
    df.dist <- beta1(df.t, 1)
    attributes(df.dist)$Labels <- rownames(df.t)
  }else{
    df.dist <- vegdist(df.t, method = metric, binary = bin)
  }
  summary(df.dist)

  # Convert distance matrix to columns
  m.dist <- as.matrix(df.dist)
  m1 <- melt(m.dist)
  m2 <- m1[(m1$Var1 != m1$Var2), ]
  m2$plot1 <- sapply(strsplit(as.character(m2$Var1), "_"), "[[", 1)
  m2$plot2 <- sapply(strsplit(as.character(m2$Var2), "_"), "[[", 1)
  m2$method1 <- sapply(strsplit(as.character(m2$Var1), "_|-"), "[[", 3)
  m2$method2 <- sapply(strsplit(as.character(m2$Var2), "_|-"), "[[", 3)
  m2$within.between <- ifelse(m2$plot1 == m2$plot2, "within", "between")
  return(m2)
}

# Load and tidy up the metadata
metadata <- read.table("Wx80_metadata/Wx80_sample_data.txt", header = T, row.names = 1, sep = "\t")
metadata$Subplot <- gsub("S2-", "Subplot ", metadata$Subplot)
metadata$Method <- gsub("Psoil", "1.5 g direct", metadata$Method)
metadata$Method <- gsub("Pmax", "7.5 g direct", metadata$Method)
metadata$Method <- gsub("Ind", "15 g indirect", metadata$Method)
metadata$Method <- gsub("PO4", "15 g PO4 buffer", metadata$Method)
metadata$Method <- factor(metadata$Method, levels = c("1.5 g direct", "7.5 g direct", "15 g indirect", "15 g PO4 buffer"),
                             ordered = TRUE)

# Now run all the multivariate analyses
genes <- c("16S", "18S", "26S", "COI")
#metric <- "jaccard"
metric <- "beta1" # Calculating beta1 distance is rather slow
wb.dists <- data.frame()

for (gene in genes) {
  #f <- Sys.glob(paste0("Wx80_OTU_tables/Wx80*", gene, "*MEE1.0_min1_OTUtable.txt"))
  f <- Sys.glob(paste0("Wx80_OTU_tables/Wx80*", gene, "*MEE1.0_min2_OTUtable.txt"))
  df <- read.table(f, sep = "\t", header = TRUE, row.names = 1, check.names = F)
  colnames(df) <- gsub("^sample_", "", colnames(df))
  rownames(df) <- gsub("size=","size_", rownames(df))
  
  # Get taxonomy data
  #tx <- Sys.glob(paste0("Wx80_taxonomy/Wx80_*", gene, "*_MEE1.0_OTUs_min1_taxonomy_v2.txt"))
  tx <- Sys.glob(paste0("Wx80_taxonomy/Wx80_*",gene,"*_MEE1.0_OTUs_min2_taxonomy_v2.txt"))
  taxa <- read.table(tx, sep="\t", header=T, row.names=1, quote = "", comment.char = "")
  taxa <- taxa[,c("superkingdom","kingdom","phylum","class")]
  
  OTUtable <- merge(df, taxa, by = "row.names")
  rownames(OTUtable) <- OTUtable$Row.names
  OTUtable$Row.names <- NULL
  
  # Split OTU table by taxonomic groups
  OTUsubset <- subset(OTUtable, !(grepl("root|cellular organisms|No hits|Not assigned", OTUtable$superkingdom)))
  Not_ass <- subset(OTUtable, grepl("root|cellular organisms|No hits|Not assigned", OTUtable$superkingdom))
  Proks <- subset(OTUsubset, grepl("Archaea|Bacteria", OTUsubset$superkingdom))
  Euks <- subset(OTUsubset, grepl("Eukaryota", OTUsubset$superkingdom))
  Fungi <- subset(OTUsubset, grepl("Fungi", OTUsubset$kingdom))
  Metazoa <- subset(OTUsubset, grepl("Metazoa", OTUsubset$kingdom))
  Plants <- subset(OTUsubset, grepl("Viridiplantae", OTUsubset$kingdom))
  Protists <- subset(Euks, !(grepl("Opisthokonta|Fungi|Metazoa|Viridiplantae", Euks$kingdom)))
  Protists <- subset(Protists, !(grepl("Eukaryota", Protists$class)))
  
  taxa.subsets <- list("all OTUs" = OTUtable, "fungi" = Fungi, "protists" = Protists, 
                       "metazoans" = Metazoa, "plants" = Plants, "prokaryotes" = Proks, "not-assigned" = Not_ass)
  
  plist <- list() # List to hold the plots
  n <- 1 # To keep track of plot number

  for(i in 1:length(taxa.subsets)){
    taxon <- names(taxa.subsets)[[i]]
    df.s <- as.data.frame(taxa.subsets[[i]])
    df.s$superkingdom <- NULL
    df.s$kingdom <- NULL
    df.s$phylum <- NULL
    df.s$class <- NULL
    
    if(nrow(df.s) > 50){ ## Ignore tables with few OTUs
      print(paste(gene, taxon, ", ", nrow(df.s), "OTUs"))
      df.sub <- df.s[, colSums(df.s) > 50] ## Exclude samples with overly low abundance
      if(g == "26S" && metric == "Jaccard"){ ## Drop some samples with anomalous results for Jaccard dist plots
        df.sub <- dfsub[, colSums(dfsub) < 48000] # Exclude two samples over threshold for 26S Jaccard dist plots
        df.sub$"S2-D_Psoil-08" <- NULL
        df.sub$"S2-D_Psoil-09" <- NULL
      }
      if(g == "COI"){
        df.sub$"S2-P_Psoil-06" <- NULL
      }
      if(ncol(df.sub) > 20){ ## Ignore any tables with few samples remaining (should be ~80)
        
        # Generate MDS ordination plot
        p <- MDS.plot(df = df.sub, n = n, metric = metric, bin = TRUE, metadata = metadata)
        legend <- get.legend(p)  #3 Extract legend from plot
        p <- p + theme(legend.position = "none")
        print(p)
        plist[[n]] <- ggplotGrob(p)  
        n <- n + 1
        
        # Beta-dispersion and PERMANOVA tests
        bdisp.PERM.test(df = df.sub, metric = metric, bin = TRUE, metadata = metadata)
        
        # Get within-between subplot distances
        wb<- within.between(df = df.sub, metric = metric, bin = TRUE)
        wb$taxon <- taxon
        wb$row.names <- NULL
        wb$gene <- gene
        wb$metric <- metric
        wb.dists <- rbind(wb.dists, wb)
      }  
    }
  }
  
  # Output MDS plots to file
  pdf(file = paste0("Wx80_", g, "_min2_", metric, "_MDS_plots.pdf"),
      width = 20/2.54, height = 27/2.54, useDingbats = FALSE)
  args.list <- c(plist, list(ncol=2, nrow=3))
  print (grid.arrange(do.call(arrangeGrob, args.list), legend, ncol=2, widths=c(1, 0.2)))
  plist <- list() # reset plot list 
  dev.off()
}  

# Save within-between subplot distances to file
write.table(wb.dists, file = "Wx80_within_between_plot_distances.txt", sep = "\t", quote = FALSE, col.names = NA)

# Make a boxplot of within-between subplot distances
wb.dists <- read.table("Wx80_within_between_plot_distances.txt", sep = "\t", header = TRUE, row.names = 1)
summary(wb.dists)
unique(wb.dists$taxon)

# Fix DNA extraction method labels
fix.methods <- function(x){
  old <- c("Psoil","Pmax","Ind","PO4")
  new <- c("1.5 g direct", "7.5 g direct", "15 g indirect", "15 g PO4 buffer")
  y <- new[match(x, old)]
  return(y)
}

# Organise the data for plotting
wb.dists$gene <- factor(wb.dists$gene, levels = c("16S","18S","26S","COI"), ordered = TRUE)
wb.dists$within.between <- factor(wb.dists$within.between, levels = c("within", "between"), ordered = TRUE)
wb.dists$taxon <- factor(wb.dists$taxon, levels = c("all OTUs","prokaryotes","protists","fungi","metazoans","plants","not-assigned"), ordered = TRUE)
wb.dists$method1 <- sapply(wb.dists$method1, function(x) fix.methods(x))
wb.dists$method2 <- sapply(wb.dists$method2, function(x) fix.methods(x))
wb.dists$method1 <- factor(wb.dists$method1, levels = c("1.5 g direct", "7.5 g direct", "15 g indirect", "15 g PO4 buffer"), ordered = TRUE)
wb.dists$method2 <- factor(wb.dists$method2, levels = c("1.5 g direct", "7.5 g direct", "15 g indirect", "15 g PO4 buffer"), ordered = TRUE)
wb.dists$gene.taxon <- paste(wb.dists$gene, wb.dists$taxon)
wb.dists$m.wb <- paste(wb.dists$method1, wb.dists$within.between)
wb.dists$m.wb <- factor(wb.dists$m.wb, levels = c("1.5 g direct within", "1.5 g direct between", "7.5 g direct within", "7.5 g direct between",
                                                  "15 g indirect within", "15 g indirect between", 
                                                  "15 g PO4 buffer within", "15 g PO4 buffer between"), ordered = TRUE)
wb.dists <- wb.dists[wb.dists$method1 == wb.dists$method2, ]
wb.dists <- data.m[data.m$taxon != "not-assigned" & data.m$taxon != "plants" & data.m$taxon != "prokaryotes", ]

# Make the figure
ggplot(wb.dists, aes(x = m.wb, y = value, colour = within.between)) + 
  geom_boxplot(outlier.colour = alpha("black", 0.15)) + scale_shape(solid = FALSE) +
  xlab("") + ylab("Distance measure") +
  scale_colour_brewer(palette = "Set1") + facet_grid(taxon ~ gene) + 
  theme(legend.position = "none", panel.grid = element_blank(), strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))

