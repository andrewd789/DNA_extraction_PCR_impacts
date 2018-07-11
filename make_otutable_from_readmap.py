# -*- coding: utf-8 -*-
"""
May 2015

@author: Andrew Dopheide

Makes a table of OTU counts from USEARCH readmap .uc file.

"""

import glob
import sys

label = sys.argv[1]  # Label from usearch pipeline

#listing = os.listdir("./")
listing1 = glob.glob("*readmap.uc")  # Readmap files

for infile in listing1:
    file1 = [x for x in listing1 if x.startswith(label)][0]
    readmap = open(file1, "r")    
    print("Processing %s..." % label)

    samples = []
    otus = []
    unmatched = []

    # Get samples and otus from readmap
    print("Getting samples and otus from readmap...")
    for row in readmap:
        fields = row.strip().split("\t")
        query = fields[8]
        #print("query: %s" % query)
        sample = query.split("|")[2]
        if sample not in samples:
            samples.append(sample)
        if fields[0] is "H":  # OTU hit
            otu = fields[9]
            if otu not in otus:
                otus.append(otu)

    # Set up bins and tables for otu counts
    bins = [[0 for row in range(len(samples))] for col in range(len(otus))]
    otutable = open(("%s_otutable.txt" % label), "w")

    # Get counts from readmap
    print("Processing readmap...")
    readmap = open(file1, "r")  
    for row in readmap:
        fields = row.strip().split("\t")
        query = fields[8]
        sample = query.split("|")[2]
        if fields[0] is "H":  # OTU hit
            otu = fields[9]
            bins[otus.index(otu)][samples.index(sample)] += 1
        elif fields[0] is "N":  # unmatched
            hit = fields[8]
            match = fields[8]
            sample = hit.split("|")[2]   
            plot = hit.split("|")[3]
            bins[otulist.index(match)][samples.index(sample)] += 1

    # Write counts to table
    for item in samples:
        otutable.write("\t%s" % item)
    n = 0
    for row in bins:
        otutable.write("\n%s" % otus[n])
        for item in row:
            otutable.write("\t%s" % item)
        n += 1

    otutable.close()
    print("Finished %s" % label)
