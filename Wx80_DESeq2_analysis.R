###############################################################################
# Calculation of differential OTU abundances among samples and DNA extraction 
# methods using Phyloseq + DESeq2
# Based on http://joey711.github.io/phyloseq-extensions/DESeq2.html
###############################################################################

library(phyloseq)
library(DESeq2)
library(BiocParallel)
register(MulticoreParam(4))
library(ggplot2)
library(data.table)

theme_set(theme_bw(base_size = 9))
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}

# Set working directory
setwd("./")

# Load taxonomy reference file, used to set plotting order of taxa
tax.ref <- read.table("wx80_metadata/Silva_taxa_order_list_v3.txt", header = TRUE, sep = "\t", quote = "", comment.char = "")

# Generate DESeq2 comparison results for each gene
genes <- c("16S", "18S", "26S", "COI")

for(gene in genes){
  f <- Sys.glob(paste0("Wx80_OTU_tables/Wx80_", gene, "*MEE1.0_min1_OTUtable.txt"))
  tx <- Sys.glob(paste0("Wx80_taxonomy/Wx80_*",gene,"*_MEE1.0_OTUs_min1_taxonomy_v2.txt"))

  OTUtab <- read.table(f, header = TRUE, row.names = 1, sep = "\t")
  rownames(OTUtab) <- gsub("=", "_", rownames(OTUtab))
  taxa <- read.table(tx, sep="\t", header = TRUE, row.names = 1, quote = "", comment.char = "")
  taxa <- taxa[,c("kingdom","phylum","class","order")]

  print(paste("starting", label, "..."))

  samples <- read.table("Wx80_metadata/Wx80_sample_data.txt", header=T, row.names=1, sep="\t")

  OTUtab <- t(OTUtab)
  OTUtab <- OTUtab[order(rownames(OTUtab)),]
  samples <- samples[order(rownames(samples)), ]
  taxa1 <- taxa[, c("kingdom", "phylum", "class", "order", "family")]

  # Convert to phyloseq object
  OTUtable <- otu_table(OTUtab, taxa_are_rows=FALSE)
  sampledata <- sample_data(samples) 
  taxatable <- tax_table(as.matrix(taxa1))
  phylo <- phyloseq(OTUtable, sampledata, taxatable)
  phylo
  
  sampledata$Method <- factor(sampledata$Method, levels = c("Psoil", "Pmax", "Ind", "PO4"), ordered = T)
  sampledata$Sample <- factor(sampledata$Sample, levels = c("D_Psoil", "D_Pmax", "D_Ind", "D_PO4",
                                                            "P_Psoil", "P_Pmax", "P_Ind", "P_PO4"), ordered = T)
  #############################################################################
  # DESeq2 analysis starts here
  #############################################################################
  
  # Convert phyloseq object to DESeq format with dispersions estimated, and test formula
  # Then test for differences (default framework)
  dds = phyloseq_to_deseq2(phylo, ~ Subplot * Method)
  dds = DESeq(dds, test="Wald", fitType="parametric", parallel = T)
  resultsNames(dds)
  
  #norm_results <- counts(dds, normalized = T)
  #write.table(norm_results, file = paste0(label, "_normalized.txt"), sep = "\t", quote = F, row.names = T)
  #rlog_results <- rlog(dds, fast = T)
  #write.table(assay(rlog_results), file = paste0(label, "_rlog.txt"), sep = "\t", quote = F, row.names = T)
  
  comparisons = list('c("Subplot", "S2-P", "S2-D")', ## Subplot D vs. P
                     'c(\"Method\", \"Pmax\", \"Psoil\")', ## 1.5 g direct vs. 7.5 g direct
                     'c(\"Method\", \"Ind\", \"Psoil\")', ## 1.5 g direct vs. 15 g indirect
                     'c(\"Method\", \"PO4\", \"Psoil\")', ## 1.5 g direct vs. 15 g phosphate buffer 
                     'c(\"Method\", \"Ind\", \"Pmax\")', ## 7.5 g direct vs. 15 g indirect 
                     'c(\"Method\", \"PO4\", \"Pmax\")', ## 7.5 g direct vs. 15 g phosphate buffer 
                     'c(\"Method\", \"PO4\", \"Ind\")') ## 15 g indirect vs. 15 g phosphate buffer 
  
  pdf(file = paste0("Wx80_", gene, "_logFC_p05_plots_v2.pdf"), width = 29.7/2.54, height = 21/2.54,
      useDingbats = F)
  dlist <- list()
  summaries <- list()
  summary_tabs <-  list()
  
  n <- 1

  for(n in 1:length(comparisons)){
    print(n)
  
    # Tabulate difference results. 
    # Contrasts: (factor, A, B) means A (numerator) vs. B (denominator) for factor (A/B). 
    # e.g. contrasts = c("Method", "Psoil", "Pmax") means Psoil compared to Pmax
    res = results(dds, cooksCutoff = FALSE, contrast = eval(parse(text = comparisons[[n]])))
    alpha = 0.05
    sigtab = res[which(res$padj < alpha), ]
    sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phylo)[rownames(sigtab), ], "matrix"))
    head(sigtab)
  
    # Re-order results by max phylum log2FC value?
    #x = tapply(sigtab$log2FoldChange, sigtab$phylum, function(x) max(x))
    #x = sort(x, TRUE)
    #sigtab$phylum = factor(as.character(sigtab$phylum), levels=names(x))
    # Order by class
    #x = tapply(sigtab$log2FoldChange, sigtab$class, function(x) max(x))
    #x = sort(x, TRUE)
    #sigtab$class = factor(as.character(sigtab$class), levels=names(x))

    # Re-organise results by taxonomic group
    sigtab$kingdom <- factor(sigtab$kingdom, ordered = TRUE, levels = rev(tax.ref$node))
    sigtab$phylum <- factor(sigtab$phylum, ordered = TRUE, levels = rev(tax.ref$node))
    sigtab$class <- factor(sigtab$class, ordered = TRUE, levels = rev(tax.ref$node))
    sigtab$taxon_group <- tax.ref$group1[match(sigtab$phylum, tax.ref$node)]
    sigtab$taxon_group <- factor(sigtab$taxon_group, ordered = TRUE, levels = c("Archaea", "Bacteria", "Eukaryota", "Opisthokonta", 
                                                              "Fungi", "Metazoa", "Viridiplantae", "unknown"))
    
    label.bits <- gsub("(c\\(|\\)|\")", "", comparisons[[n]])
    x1 <- sapply(strsplit(label.bits, ","), "[[", 1)
    x2 <- sapply(strsplit(label.bits, ","), "[[", 2)
    x3 <- sapply(strsplit(label.bits, ","), "[[", 3)

    sigtab$numerator <- x2
    sigtab$denominator <- x3
    dlist[[n]] <- sigtab
    summaries[[n]] <- capture.output(summary(res))
 
    # Make a summary table of results
    dt <- data.table(sigtab)
    down <- dt[log2FoldChange < 0, .("nOTUs"=length(baseMean)), by=.(taxon_group,phylum,class)]
    up <- dt[log2FoldChange > 0, .("nOTUs"=length(baseMean)), by=.(taxon_group,phylum,class)]
    down$x_y <- paste("lower in", x2, "vs", x3) 
    up$x_y <- paste("higher in", x2, "vs", x3)
    summary_tab <- rbind(up, down)
    summary_tab$gene <- paste(gene)
    summary_tabs[[n]] <- summary_tab
       
    # Drop unwanted taxonomic groups for graph
    sigtab2 <- sigtab[!(grepl("No hits|Not assigned|cellular organisms|root", sigtab$kingdom)), ]
    if(gene != "16S"){ # Exclude any prokaryotes from non-16S results
      sigtab2 <- sigtab2[!(grepl("Archaea|Bacteria", sigtab2$kingdom)), ]
    }
    if(gene == "16S"){ # Only include prokaryotes in 16S results
      sigtab2 <- sigtab2[grepl("Archaea|Bacteria", sigtab2$kingdom), ]
    }
      
    ## Make a graph of OTUS with log-fold differences between subplots or methods
    title <- paste0(gene, ", ", x2, " compared to ", x3, ", adj.p < ", alpha)
    title <- gsub("S2-D", "subplot D", title)
    title <- gsub("S2-P", "subplot P", title)
    title <- gsub("Psoil", "1.5 g direct", title)
    title <- gsub("Pmax", "7.5 g direct", title)
    title <- gsub("Ind", "15 g indirect", title)
    title <- gsub("PO4", "15 g PO4 buffer", title)

    p <- ggplot(na.omit(sigtab2), aes(x=class, y=log2FoldChange, color=taxon_group)) + 
            geom_point(size=2, alpha = 0.75, shape = 5) + 
            scale_y_continuous(breaks = c(-10,-8,-6,-4,-2,0,2,4,6,8,10)) +
            guides(colour = guide_legend(override.aes = list(alpha = 1))) +
            ggtitle(title) + ylab("Log 2 fold-change") + xlab("") +
            scale_colour_brewer(palette="Set1") +
            coord_flip()
    print(p)
    
    print("printing plot")
    n <- n + 1
    
  }
  
  dev.off()
  
  ## Write data to file
  dlist.all <- rbindlist(dlist)
  write.table(dlist.all, file = paste0("Wx80_", gene, "_logFC_p05_data.txt"), 
              sep = "\t", quote = F, row.names = T)
  summary_tabs.all <- rbindlist(summary_tabs)
  write.table(summary_tabs.all, file = paste0("Wx80_", gene, "_logFC_p05_summary_table.txt"), 
              sep = "\t", quote = F, row.names = T)
}
