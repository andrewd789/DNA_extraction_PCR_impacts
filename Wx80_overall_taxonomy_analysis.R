###############################################################################
# Overall analyses of taxonomy composition per gene
###############################################################################

library(data.table)
library(ggplot2)
library(gridExtra)
library(plyr)
library(RColorBrewer)
library(reshape2)
library(scales)
library(vegan)
library(vegetarian)

theme_set(theme_bw(base_size=8))

# Set up an Colour palette
colors <- brewer.pal(8, "Spectral")
pal <- colorRampPalette(colors) 

setwd("./")

# Compile overall richness and abundance stats per gene and phylum
genes <- c("16S", "18S", "26S", "COI")
tx.all <- data.frame()
d1 <- data.frame()
for(gene in genes){
  for(m in c("min1","min2")){
    print(paste(gene, m))
    if(m == "min1"){
      f <- Sys.glob(paste0("Wx80_OTU_tables/Wx80_", gene, "*MEE1.0_min1_OTUtable.txt"))
      tx <- Sys.glob(paste0("Wx80_taxonomy/Wx80_*", gene, "*_MEE1.0_OTUs_min1_taxonomy_v2.txt"))
    }else{
      f <- Sys.glob(paste0("Wx80_OTU_tables/Wx80_", gene, "*MEE1.0_min2_OTUtable.txt"))
      tx <- Sys.glob(paste0("Wx80_taxonomy/Wx80_*",gene,"*_MEE1.0_OTUs_min2_taxonomy_v2.txt"))
    }
    OTUtable <- read.table(f, header = TRUE, row.names = 1, sep = "\t")
    taxa <- read.table(tx, sep="\t", header = TRUE, row.names = 1, quote = "", comment.char = "")
    taxa <- taxa[,c("kingdom","phylum","class","order")]

    # Fix OTU names so they match with taxonomy
    rownames(OTUtable) <- lapply(rownames(OTUtable), FUN = function(x) gsub("=", "_", x))
    
    # Get OTU counts and abundances by phylum
    OTUsums <- as.data.frame(rowSums(OTUtable)) # Abundances (sums)
    colnames(OTUsums)[[1]] <- "sum"
    OTUsums <- merge(OTUsums, taxa, by = "row.names")
    OTUsums$phylum <- factor(OTUsums$phylum) # Drop unneeded levels
    dt <- data.table(OTUsums)
    # OTUs per phylum
    counts <- dt[, length(sum), by = phylum] 
    colnames(counts)[[2]] <- "value"
    counts$var <- "OTUs"
    counts$gene <- gene
    counts$min <- m
    # Total abundances per phylum
    sums <- dt[, sum(sum), by = phylum]
    colnames(sums)[[2]] <- "value"
    sums$var <- "Seqs"
    sums$gene <- gene
    sums$min <- m
    d1 <- rbind(d1, counts, sums)
    tx.all <- rbind(tx.all, taxa)
    
    # Get effective alpha diversity by phylum?
    #OTUtable_tax <- merge(OTUtable, taxa, by = "row.names")
    #rownames(OTUtable_tax) <- OTUtable_tax$Row.names
    #OTUtable_tax$Row.names <- NULL
    #OTUtable_tax$phylum <- factor(OTUtable_tax$phylum) # Drop unneeded levels
    #div <- data.frame()
    #for(lev in levels(OTUtable_tax$phylum)){
    #  print(lev)
    #  tax_subset <- OTUtable_tax[grep(lev, OTUtable_tax$phylum), ]
    #  tax_subset <- tax_subset[, 1:80]
    #  g1 <- d(t(tax_subset), lev = "gamma", q = 1)
    #  d <- data.frame("phylum" = lev, "value" = g1)
    #  div <- rbind(div, d)
    #}
    #div$var <- "Alpha1"
    #div$gene <- gene
    #div$min <- m
    #d1 <- rbind(d1, div)
  }
}  

# Ensure results for all phyla are present for all genes, even if there are no values
d2 <- dcast(d1, phylum ~ var + gene + min) # Makes phyla consistent across all variables (empty values = NA)
d2 <- melt(d2, id.vars = "phylum")

# Replace NAs with zero?
d2[is.na(d2)] <- 0

# Get factors for plotting
d2$gene <- sapply(strsplit(as.character(d2$variable), "_"), "[[", 2)
d2$var <- sapply(strsplit(as.character(d2$variable), "_"), "[[", 1)
d2$min <- sapply(strsplit(as.character(d2$variable), "_"), "[[", 3)

# Organise the taxonomy data for plotting
# Load a reference taxonomy, used to determine plotting order
tax_ref <- read.table("Wx80_metadata/Silva_taxa_order_list_v3.txt", header = TRUE, sep = "\t", quote = "", comment.char = "")
tx.all <- unique(tx.all)
d2$kingdom <- tx.all$kingdom[match(d2$phylum, tx.all$phylum)]
d2$group1 <- tax_ref$group1[match(d2$phylum, tax_ref$node)] # Add a high-level group for summarising/plotting
d2$group1 <- factor(d2$group1, ordered = TRUE, 
                    levels = c("Archaea","Bacteria","Eukaryota","Fungi","Metazoa","Opisthokonta","Viridiplantae","unknown"))
tax_ref$node <- gsub("Bryophyta <mosses>", "Bryophyta", tax_ref$node) # Fix problem in tax_ref
d2$group1[d2$phylum == "Bryophyta"] <- "Viridiplantae"

# Sort out eukaryote taxonomic groups  - which are protists, and which are "others"?
protists <- as.character(unique(d2$phylum[grepl("Eukaryota|Opisthokonta", d2$group1)]))
protists <- protists[!grepl("Eukaryota|unclassified eukaryotes|Opisthokonta", protists)]
others <- c("Eukaryota", "Opisthokonta", "Opisthokonta incertae sedis", "unclassified eukaryotes")
d2$group2 <- as.character(d2$group1)
d2$group2[d2$phylum %in% protists] <- "Protists"
d2$group2[d2$phylum %in% others] <- "Other eukaryotes"
d2$group2 <- factor(d2$group2, ordered = TRUE, 
                    levels = c("Archaea","Bacteria","Protists","Fungi","Metazoa","Viridiplantae","Other eukaryotes","unknown"))

phy.ord <- tax_ref$node[tax_ref$node %in% d2$phylum]
d2$phy.ord <- match(d2$phylum, phy.ord)
d2 <- d2[order(d2$group2, d2$group1, d2$phy.ord), ]
d2$phylum <- factor(d2$phylum, ordered = TRUE, levels = rev(unique(d2$phylum)))

# Output to file
#write.table(d2, file = "Wx80_overall_otus_seqs_by_phylum.txt", sep = "\t", quote = F, row.names = F)

# Determine which "phyla" aren't really phylum-level groups, so they can be excluded from phylum counts,
# and coloured accordingly on the y-axis of figure below. This is tricky, as many protist and plant groups are "no rank".
nonphy <- as.character(unique(d2$phylum[!(d2$phylum %in% tax_ref$node[tax_ref$rank == "phylum"])])) # Includes various phylum-equivalent groups
nonphy <- nonphy[!grepl(paste0("Amoebozoa|Ciliophora|Embryophyta|Euphyllophyta|Magnoliophyta|Apusozoa|Centroheliozoa|",
                               "Cercozoa|Choanoflagellida|Cryptophyta|Dinophyceae|Euglenozoa|Heterolobosea|Oxymonadida|",
                               "Rhodophyta|Jakobida|Malawimonadidae|Bryophyta|Marchantiophyta"), nonphy)]
d2$phy <- ifelse(d2$phylum %in% nonphy, "non-phylum","phylum")

# Summarise OTU, sequence, and phylum counts as a table
d2.tab <- dcast(d2, var + group2 ~ min + gene, sum)
nphy.tab <- dcast(d2[d2$var == "OTUs" & !(d2$phylum %in% nonphy) & d2$value > 0, ], var + group2 ~ min + gene, length)
nphy.tab$var <- "phyla"
d2.tab <- rbind(d2.tab, nphy.tab)
#write.table(d2.tab, file = "Wx80_otus_seqs_nphy_by_group_gene_min.txt", sep = "\t", quote = F, row.names = F)  

# Make a figure showing overall taxonomy composition by gene
# Make a new variable for easy selection of data components
d2$var2 <- paste0(d2$var, "_", d2$min)

# Exclude Seqs min2 data
d2 <- d2[!d2$var2 == "Seqs_min2", ]

# Change any gene/phylum data points where all three of OTUs min1, Seqs min1, and OTUs min2 are zero to NA
# This prevents these points from being plotted at zero
for(x in unique(paste(d2$gene, d2$phylum))){
  #print(x)
  if(all(d2$value[paste(d2$gene, d2$phylum) == x & d2$var2 == "OTUs_min1"] == 0, 
         d2$value[paste(d2$gene, d2$phylum) == x & d2$var2 == "Seqs_min1"] == 0,
         d2$value[paste(d2$gene, d2$phylum) == x & d2$var2 == "OTUs_min2"] == 0)){
  d2$value[paste(d2$gene, d2$phylum) == x] <- NA  
  }
}

# Fix number format on the x-axis
scientific_10 <- function(x) {
  #parse(text=gsub("e", " %*% 10^", scientific_format()(x)))
  parse(text=gsub("1e\\+0", "10^", scientific_format()(x)))
}

# Now make the figure
p <- ggplot(d2) + 
  geom_point(data = d2[grepl("OTUs_min1", d2$var2), ], aes(x = phylum, y = value, colour = group2), shape = 1, size = 2) +
  geom_point(data = d2[grepl("Seqs_min1", d2$var2), ], aes(x = phylum, y = value, colour = group2), shape = 5, size = 2) +
  geom_point(data = d2[grepl("OTUs_min2", d2$var2), ], aes(x = phylum, y = value, colour = group2), shape = 16, size = 2) +
  #geom_point(data = d2[var == "Seqs_min2",], aes(x = phylum, y = value, colour = group2), shape = 17, size = 2) +
  geom_line(data = d2, aes(x = phylum, y = value, group=interaction(phylum, gene), 
                           colour = group2), size = 0.5, alpha = 0.5) +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000), label=scientific_10) +
  xlab("") + ylab("Number of OTUs/sequences") + 
  #scale_color_discrete(name = "Taxonomic group") + # Use default colour palette
  scale_colour_brewer(palette = "Set2", name = "Taxonomic group") + # Use a ColourBrewer palette
  facet_grid( ~ gene) + coord_flip() + 
  theme(strip.background = element_blank(), axis.text.y = element_text(
    colour=ifelse(levels(d2$phylum) %in% nonphy, "#7a7a7a", "black")), axis.text.x = element_text(size = 7))
p

# Output to pdf
ggsave(p, file = "Wx80_overall_taxonomy_OTUs-reads-by-phylum.pdf", width = 260, height = 160, units = "mm")
