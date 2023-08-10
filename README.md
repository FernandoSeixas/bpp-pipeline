# bpp-pipeline
scripts to prepare files and run BPP analysis

## Brief summary of the pipeline
This pipeline has three steps:
#1. Prepare bed files defining genomic regions (either coding or noncoding). These will be the loci that will be used for the BPP analysis.\
#2. Prepare alignments for the loci defined in the previous step, from VCF files.\
#3. Group loci alignments into blocks, prepare BPP .ctl files and run BPP.\
