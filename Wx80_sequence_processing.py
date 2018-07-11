# -*- coding: utf-8 -*-

###############################################################################
# Bioinformatic processing of multi-gene MiSeq data.
# First, use PEAR paired end read merger to merge paired sequence files.
# PEAR read merger is called by "pear", requires sequence files to NOT be gzipped.
# This only works well for 18S and COI sequences, as very few 16S and 26S sequences overlap.
###############################################################################

import glob, os, subprocess, re

def run_PEAR(R1_file, R2_file, outfile, log):
    start_PEAR = "pear -f {0} -r {1} -v 50 -o {2}".format(R1_file, R2_file, outfile)
    print(start_PEAR)
    subprocess.call(start_PEAR.split(), stdout = log, stderr = subprocess.STDOUT)

os.chdir("NZGL01278_data") # Path to raw sequence data
log = open("PEAR_log.txt", "a")
R1_files = glob.glob("*R1_001.fastq")
for R1_file in R1_files:
    print("R1 file: {0}".format(R1_file))
    R2_file = glob.glob("{0}*R2_001.fastq".format(R1_file[:13]))[0]
    print("R2 file: {0}".format(R2_file))
    outfile = "{0}_pear_contigs".format(R1_file[:13])
    print("Outfile {0}".format(outfile))
    run_PEAR(R1_file, R2_file, outfile, log) 

log.close()

###############################################################################
# Next, split the merged (18S, COI) and/or unmerged (16S, 26S) sequence files 
# by primer, and relabel the sequences by sample.
###############################################################################

from Bio import SeqIO, SeqUtils # This part uses Biopython

def split_relabel(infile, gene, m, sample_IDs, output_handle):
    ## Set up logfile
    log = open("Primer_split_log.txt","a")
    log.write("Reads file\tPrimer\tRead count\n")
    if m == "merged":
        label = re.split('Wx80_|_pear', infile)[1]
    elif m == "unmerged":
        label = re.split('Wx80_|_L001', infile)[1]
    reads_handle = open(infile)
    #reads_handle = gzip.open(fastq_file) # Use if files are gzipped
    well = re.split('_S', label)[0] # NZGL well number
    sample_ID = sample_IDs.get(well) # Corresponding sample ID
    count = 0
    primer_f = primerlist.get(gene)[0]
    primer_r = primerlist.get(gene)[1]
    trimmed = ""
    for record in SeqIO.parse(reads_handle, "fastq"):
        primer_search = SeqUtils.nt_search(str(record.seq), primer_f)  # Searches record.seq for primer
        if len(primer_search) > 1 and (primer_search)[1] == 0: # Check if primer found (len > 1) at start of sequence ([1] == 0) 
            if m == "merged":
                trimmed = record[len(primer_f):-len(primer_r)]
            elif m == "unmerged":
                trimmed = record[len(primer_f):]
            trimmed.id = (("{0}|gene_{1}|sample_{2}").format(trimmed.id, gene, sample_ID)) # Adds gene and sample ID
            print(trimmed.id)
            print(len(trimmed.seq))
            SeqIO.write(trimmed, output_handle, "fastq")
            count += 1
            print("{0} {1} {2} {3}".format(gene, label, sample_ID, count))
    print("{0} {1}: Saved {2} reads".format(label, gene, count))
    log.write("{0} {1} {2}: {3}".format(gene, label, sample_ID, count))
    reads_handle.close()
    #output_handle.close()
    log.close()

# Load a tab-delimited list of gene names and primer sequences.
primerlist = {}
primerfile = open("../Wx80_metadata/Barcode_split_primers.txt", "r")
for row in primerfile:
    gene, primer_f, primer_r = row.strip().split("\t")
    primerlist[gene] = primer_f, primer_r
    print ("{0} {1} {2}".format(gene, primer_f, primer_r))

# Load a tab-delimited list of sample IDs
sample_IDs = {}
ID_file = open("../Wx80_metadata/NZGL01278_Sample_ID.txt", "r")
for row in ID_file:
    sample, name = row.strip().split("\t")
    if name.startswith("S2-"): # Exclude any non-Wx80 samples
        sample_IDs[sample] = name
        print ("{0} {1}".format(sample, name))

log = open("Primer_split_log.txt","a")
log.write("Reads file\tPrimer\tRead count\n")
log.close()

# Split and relabel all the sequences
for gene in "18S","COI":
    outfile = open("Wx80_{0}_all_merged.fastq".format(gene), "a")
    infiles = glob.glob("*_pear*.fastq")
    for infile in infiles:
        print ("Processing {0}".format(infile))
        split_relabel(infile, gene, m = "merged", sample_IDs = sample_IDs, output_handle = outfile)
    outfile.close()

for gene in "16S","26S":
    outfile = open("Wx80_{0}_all_R1_reads.fastq".format(gene), "a")
    infiles = glob.glob("*_R1_001.fastq")
    for infile in infiles:
        print ("Processing {0}".format(infile))
        split_relabel(infile, gene, m = "unmerged", sample_IDs = sample_IDs, output_handle = outfile)
    outfile.close()

primerfile.close()
sample_IDs.close()

###############################################################################
# Then process the each amplicon file into OTUs using USEARCH.
# For merged sequences, trim to length appropriate for each target.
# For unmerged R1 reads, trim to 250 bp.
###############################################################################

#usearch = "../usearch7.0.1090_i86linux64"
usearch = "../Usearch/usearch8.0.1623_i86linux64"

def filter_trim_merged(input_file, label, minlen, maxee, log):
    # Trim and quality filter merged sequences
    filter_trim = "{0} -fastq_filter {1} -fastq_minlen {2} -fastq_maxee {3} -fastaout {4}_filter_trim.fasta".format(usearch, input_file, minlen, maxee, label)
    subprocess.call(filter_trim.split(), stdout=log, stderr=subprocess.STDOUT)
    log.write('\nFilter and trim reads:\n' + str(filter_trim) + '\n')

def filter_trim_reads(input_file, label, trunclen, maxee, log):
    # Trim and quality filter unmerged R1 reads
    filter_trim = "{0} -fastq_filter {1} -fastq_trunclen {2} -fastq_maxee {3} -fastaout {4}_filter_trim.fasta".format(usearch, input_file, trunclen, maxee, label)
    subprocess.call(filter_trim.split(), stdout=log, stderr=subprocess.STDOUT)
    log.write('\nFilter and trim reads:\n' + str(filter_trim) + '\n')

def get_usearch_OTUs(label, minsize, log):
    # Remove duplicate sequences
    derep = "{0} -derep_fulllength {1}_filter_trim.fasta -fastaout {1}_uniques.fasta -sizeout".format(usearch, label)
    subprocess.call(derep.split(), stdout=log, stderr=subprocess.STDOUT)
    log.write("\nDereplicate filtered/trimmed reads:\n" + derep + "\n")
    # Sort sequences by size (abundance)
    sortbysize = "{0} -sortbysize {1}_uniques.fasta -minsize {2} -fastaout {1}_uniques_sorted.fasta".format(usearch, label, minsize)
    subprocess.call(sortbysize.split(), stdout=log, stderr=subprocess.STDOUT)
    log.write("\nSort reads by size, minsize {0}".format(minsize) + "\n" + sortbysize + "\n")
    # Cluster the filtered, dereplicated, and sorted sequences into OTUs
    clusterOTUs = "{0} -cluster_otus {1}_uniques_sorted.fasta -otus {1}_OTUs.fasta -relabel OTU -sizein -sizeout -fastaout {1}_OTUS_uparse.fasta".format(usearch, label)
    subprocess.call(clusterOTUs.split(), stdout=log, stderr=subprocess.STDOUT)
    log.write("\nCluster sequences into OTUs at 97%:\n" + clusterOTUs + "\n")
    # Filter the OTUs for de novo chimeras
    chimera_denovo = "{0} -uchime_denovo {1}_OTUs.fasta -uchimeout {1}_OTUs_filtered.uchime -nonchimeras {1}_OTUs_good.fasta".format(usearch, label)
    s = str(chimera_denovo)
    log.write("\n Chimera filter OTUs: \n" + s + "\n")
    subprocess.call(chimera_denovo.split(), stdout=log, stderr=subprocess.STDOUT)
    # Map the filtered sequence reads to the OTUs
    searchReads = "{0} -usearch_global {1}_filter_trim.fasta -db {1}_OTUs_good.fasta -strand plus -id 0.97 -uc {1}_OTUs_num_readmap.uc".format(usearch, label)
    subprocess.call(searchReads.split(), stdout=log, stderr=subprocess.STDOUT)
    log.write("\nMap the filtered/trimmed reads to the OTUs\n" + searchReads + "\n")
    # Make an OTU abundance table
    makeOTUtable = "python make_otutable_from_readmap.py {1}_OTUs_num_readmap".format(label)
    subprocess.call(makeOTUtable.split(), stdout=log, stderr=subprocess.STDOUT)

maxee_values = [2.5, 1.0, 0.5]
minsize_values = [1, 2]

genes = ["16S", "18S", "26S", "COI"]

for gene in genes:
    f = glob.glob("Wx80_{0}_all_*fastq".format(gene))[0]
    m = f.split("_")[3].split(".fastq")[0]
    if gene == "18S":
        minlen = 310
    elif gene == "COI":
        minlen = 300
    elif gene == "16S" or gene == "26S":
        trunclen = 250
    # To run with different maxee filtering values:
    #for maxee in maxee_values:
    # Or to run with different minsize values:
    for minsize in minsize_values:
        maxee = 1.0
        #minsize = 2
        label = ("Wx80_{0}_{1}_MEE{2}_min{3}".format(gene, m, maxee, minsize))
        log = open("{0}_usearch_log.txt".format(label), "a")
        if m == "merged":
            filter_trim_merged(f, label, minlen, maxee, log)
            get_usearch_OTUs(label, minsize, log)
        elif m == "R1":
            filter_trim_reads(f, label, trunclen, maxee, log)
            get_usearch_OTUs(label, minsize, log)
        os.remove("{0}_filter_trim.fasta".format(label))
        os.remove("{0}_uniques.fasta".format(label))
        os.remove("{0}_uniques_sorted.fasta".format(label))
