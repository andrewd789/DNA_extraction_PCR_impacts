
**Sequence processing and analysis code used to investigate the effects of different DNA extraction methods and PCR variability on multi-taxon DNA meta-barcoding estimates of soil biodiversity.**

Citation:
Andrew Dopheide, Dong Xie, Thomas R. Buckley, Alexei J. Drummond, and Richard D. Newcomb, 2018. **Impacts of DNA extraction and PCR on DNA metabarcoding estimates of soil biodiversity.**  Methods in Ecology and Evolution.

Illumina sequence data for this project is available from the NCBI Sequence Read Archive under accession SRP148718. 
The data contains prokaryote 16S, eukaryote 18S, fungal 26S, and metazoan COI sequences derived from soil DNA extracted from two subplots using four different methods and amplified in sets of ten individually-indexed PCRs. 

The raw sequence data was merged (where possible), and processed into OTUs using Python with PEAR read merger (Zhang et al. 2014) and USEARCH (Edgar 2013): 
- Wx80_sequence_processing.py
- make_otutable_from_readmap.py 

Analysis of the OTUs and generation of results and figures were carried out using R:
- Wx80_overall_taxonomy_analysis.R (generation of overall taxonomy composition results)
- Wx80_effects_DNA_extraction_method.R (analysis of biodiversity results from different DNA extraction methods)
- Wx80_DESeq2_analysis.R (testing for differing OTU abundances between DNA extraction methods and subplots) 
- Wx80_effects_of_PCR_variability.R (analysis of effects of multiple PCRs on biodiversity estimates)
