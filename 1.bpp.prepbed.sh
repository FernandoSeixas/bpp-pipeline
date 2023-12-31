## BPP - part 1
## This script creates .bed files defining coding and non-coding regions for which the alignments will be created
## Software requirements: bedtools2


## Genomic coordinates *************************************************************************************************
mkdir bedFiles

#/ Rerence arguments
refnam=$1
genome="support/$refnam.genome.bed"                  # .bed file with coordinates of every scaffold/support/chromosome in the reference (e.g. Hmel201001o 1 17206585) 
exons="support/$refnam.exons.bed"                    # .bed file with coordinates of exons as annotated in the reference genome
scaln="support/$refnam.scaffoldsLength.noheader.txt" # same as the [genome] file but with only scaffold names and end position (e.g. Hmel201001o 17206585) 

#/ loci arguments
NCminlen=$2;  export NCminlen=$NCminlen           # noncoding minimum locus length
NCmaxlen=$3;  export NCmaxlen=$NCmaxlen           # noncoding maximum locus length
NCgaplen=$4;  export NCgaplen=$NCgaplen           # noncoding minimum gap between adjacent loci

CDminlen=$5;  export CDminlen=$CDminlen           # coding minimum locus length
CDmaxlen=$6;  export CDmaxlen=$CDmaxlen           # coding maximum locus length
CDgaplen=$7;  export CDgaplen=$CDgaplen           # coding minimum gap between adjacent loci


#/ Create buffer around exons
egrep -v "Hmel200" $exons > tmp; mv tmp bedFiles/hmelv25.exons.inchrom.bed 
~/software/bedtools2/bin/bedtools slop -b 2000 -i bedFiles/hmelv25.exons.inchrom.bed -g ${scaln} > bedFiles/hmelv25.exons.buffer_2kb.bed         # create a buffer of 2kb around exons to avoid linked selection
sort -k1,1 -k2,2n bedFiles/hmelv25.exons.buffer_2kb.bed > bedFiles/hmelv25.exons.buffer_2kb.sorted.bed                                           # sort bed file
~/software/bedtools2/bin/bedtools merge -i bedFiles/hmelv25.exons.buffer_2kb.sorted.bed > bedFiles/hmelv25.exons.buffer_2kb.merge.bed            # merge any overlapping blocks


## Create bed files for NONCOD (intergenic + intronic) and CODING (exonic) regions *************************************

## NONCOD //////
#/ create bed file non coding [NONCOD] regions
~/software/bedtools2/bin/bedtools complement -i bedFiles/hmelv25.exons.buffer_2kb.merge.bed -g $scaln > bedFiles/hmelv25.noncod.Step1.bed
egrep -v "Hmel200" bedFiles/hmelv25.noncod.Step1.bed > tmp; mv tmp bedFiles/hmelv25.ncod.Step1.bed
# run Rscript to get coordinates of valid regions to sample from
Rscript \
    ~/code/heliconius_seixas/2.elevatus_pardalinus/2.SpeciesTree/bpp-pipeline/DefineLocusToSelect.R \
    -d bedFiles/ \
    -i ${refnam}.noncod.Step1.bed \
    -l ${NCminlen} \
    -m ${NCmaxlen} \
    -g ${NCgaplen} \
    -o ${refnam}.noncod.minL${NCminlen}.maxL${NCmaxlen}.minG${NCgaplen}.bed

## EXONIC //////
# get exonic unique coordinates (sorted), and only in chromosomes
sort -k1,1 -k2,2n bedFiles/hmelv25.exons.inchrom.bed > bedFiles/${refnam}.exonic.Step1.bed
# run Rscript to get coordinates of valid regions to sample from
Rscript \
    ~/code/heliconius_seixas/2.elevatus_pardalinus/2.SpeciesTree/bpp-pipeline/DefineLocusToSelect.R \
    -d bedFiles/ \
    -i ${refnam}.exonic.Step1.bed \
    -l ${CDminlen} \
    -m ${CDmaxlen} \
    -g ${CDgaplen} \
    -o ${refnam}.exonic.minL${CDminlen}.maxL${CDmaxlen}.minG${CDgaplen}.bed
