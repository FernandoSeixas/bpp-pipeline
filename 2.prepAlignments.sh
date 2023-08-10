## BPP - part 2
## This script creates alignment files from vcf files for the genomic regions defined by the previous script (both exonic and noncod)
## Software requirements: bedtools2, bcftools


## variables & export
#/ parameters used in the previous scrip
locitp=$1; export locitp=$locitp    # type of loci [exonic, noncod]
minlen=$2; export minlen=$minlen    # minimim size of locus
maxlen=$3; export maxlen=$maxlen    # maximum size of locus
gaplen=$4; export gaplen=$gaplen    # minimim distance between loci
setloc=`echo "noncod.minL${minlen}.maxL${maxlen}.minG${gaplen}"`
export setloc=$setloc
#/ Other parameters
refnam=$5; export refnam=$refnam        # reference code
sgroup=$6; export sgroup=$sgroup        # species group
scalen="support/$refnam.scaffoldsLength.noheader.txt"
repeat="support/hmelv25_repeats.txt"    # repeat annotations can be found in lepbase
export scalen=$scalen
export repeat=$repeat


## Files that need to be prepared manually
examples/${sgroup}.inds.list             # file with the names of individuals to be analyzed (as in the vcf file). One individual per line  (e.g. indv1)
examples/${sgroup}.popmap.txt            # tab delimited file with individual names and population to assign (e.g indv1  pop1)
examples/${sgroup}.popmin.txt            # tab delimited file with population names and minimum number of individuals passing filters per population (e.g. pop1  2)


## Prepare loci alignments *********************************************************************************************
#/ create directories
mkdir $setloc.$sgroup
mkdir $setloc.$sgroup/1.vcfRegions
mkdir $setloc.$sgroup/2.vcf2fas
mkdir $setloc.$sgroup/3.fltFasta
mkdir $setloc.$sgroup/4.bppFormat

#/ Extract vcfs: extract regions from vcf [and mask repeat regions]
vcfdir="/n/mallet_lab/Lab/fseixas/1.projects/2.elevatus_pardalinus/0.data/2.3.merge/"; export vcfdir=$vcfdir

# make sure to remove any residual files that may have been created in previous runs 
rm $setloc.$sgroup/1.vcfRegions/*.vcf

cat bedFiles/$refnam.$setloc.bed | egrep Hmel201 | xargs -n 3 -P 8 sh -c '
    ~/software/samtools-1.17/bcftools-1.17/bcftools view \
    -r $0:$1-$2 \
    --samples-file examples/${sgroup}.inds.list \
    --targets-file ^$repeat \
    $vcfdir/$0.hmelv25.merge.vcf.gz \
    > $setloc.$sgroup/1.vcfRegions/$0.$1.$2.vcf
    # check if file as enough lines (here a minimum of 10 sites)
    nsnps=`~/software/samtools-1.17/bcftools-1.17/bcftools view -H $setloc.$sgroup/1.vcfRegions/$0.$1.$2.vcf | wc -l`
    if [ $nsnps -lt 10 ] 
    then 
        rm $setloc.$sgroup/1.vcfRegions/$0.$1.$2.vcf
    fi
'

#/ Convert vcf to fasta
rm $setloc.$sgroup/2.vcf2fas/*.fasta
ls $setloc.$sgroup/1.vcfRegions/ | sed 's,.vcf,,g' | xargs -n 1 -P 8 sh -c '
    echo $0;
    perl vcf2fas.iupac.pl \
    $setloc.$sgroup/1.vcfRegions/$0.vcf \
    $setloc.$sgroup/2.vcf2fas/$0.fasta
'

#/ Filter alignments both vertically [genotyping across individuals] and horizontally [genotyping in each individual]. These filtes are:
# 1. min number of individuals passing filters per population
# 2. max proportion of Ns per sequence [horizontal]. FLOAT
# 3. max proportion of Ns per position [vertical]. FLOAT
# 4. min length of alignment. INTEGER
rm $setloc.$sgroup/3.fltFasta/*.filter.fasta
for scaffold in `cut -f1 $scalen `; do
    echo ${scaffold}
    ls $setloc.$sgroup/2.vcf2fas/${scaffold}*.fasta | \
    sed "s,$setloc.$sgroup/2.vcf2fas/,,g" | \
    sed "s,.fasta,,g" | xargs -n 1 -P 1 sh -c '
    perl filterFastaAlignments.pl \
        $setloc.$sgroup/2.vcf2fas/$0.fasta \
        $setloc.$sgroup/3.fltFasta/$0.filter.fasta \
        examples/$sgroup.popmap.txt \
        examples/$sgroup.popmin.txt \
        0.50 0.20 50
    '
done

#/ Convert to BPP format
rm $setloc.$sgroup/4.bppFormat/*.fasta
ls $setloc.$sgroup/3.fltFasta/* | sed 's,.filter.fasta,,g' | sed "s,$setloc.$sgroup/3.fltFasta/,,g" | xargs -n 1 -P 8 sh -c '
    perl fasta2bpp.pl \
        $setloc.$sgroup/3.fltFasta/$0.filter.fasta \
        $setloc.$sgroup/4.bppFormat/$0.bpp.fasta \
        $0 \
        10
'