#!/bin/bash
#SBATCH -n 8                              # Number of cores requested
#SBATCH -N 1                              # Ensure that all cores are on one machine
#SBATCH -t 300                            # Runtime in minutes
#SBATCH -p serial_requeue,shared          # Partition to submit to
#SBATCH --mem-per-cpu=1000                # Memory per cpu in MB (see also --mem-per-cpu)
#SBATCH --open-mode=append               
#SBATCH -o BPP_%j.out                     # Standard out goes to this file
#SBATCH -e BPP_%j.err                     # Standard err goes to this filehostname


## Run BPP analysis in blocks ******************************************************************************************

# reference arguments
refnam="hmelv25";   export refnam=$refnam       #/ reference name
scalen="support/$refnam.scaffoldsLength.noheader.txt"
# loci arguments
nlocus="100";       export nlocus=$nlocus       #/ nlocus: number of loci per block (e.g. 100)
setloc=$1;          export setloc=$setloc       #/ setloc: same as in the previous script (e.g. noncod.minL100.maxL300.minG2000)
sgroup=$2;          export sgroup=$sgroup       #/ sgroup: same as in the previous script (e.g. example1)
# bpp arguments (these are things that need to be changed in BPP .ctl file and change from dataset to dataset). Read BPP manual for details
sline=$3;           export sline=$sline      #/ sline: species line (e.g. 4 parAnd parAma eleAma eleGui )
iline=$4;           export iline=$iline      #/ iline: number individuals/pop line (e.g. 2 2 2 2)
tline1=$5;          export tline1=$tline1    #/ tline1: initial tree rep1 - e.g  "((parAnd,parAma),(eleAma,eleGui));" 
tline2=$6;          export tline2=$tline2    #/ tline2: initial tree rep2 - e.g  "(parAnd,(parAma,(eleAma,eleGui)));" 
tline3=$7;          export tline3=$tline3    #/ tline3: initial tree rep3 - e.g  "(parAma,(parAnd,(eleAma,eleGui)));" 

## create directories
mkdir $setloc.$sgroup/5.0.blocks/
mkdir $setloc.$sgroup/5.1.blocksAlign/
mkdir $setloc.$sgroup/5.2.bpp
mkdir $setloc.$sgroup/logs


## Prepare blocks of loci **********************************************************************************************
cd ${setloc}.${sgroup}/
#/ make sure there are no files from previous analyses
rm 5.0.blocks/* 
rm 5.1.blocksAlign/* 
rm 5.2.bpp/*

#/ list loci to concatenate
cd 4.bppFormat
ls *.bpp.fasta > ../finalLoci.list
cd ..
sed 's,\.,\t,g' finalLoci.list | sort -k1,1 -k2,2n | sed 's,\t,.,g' > tmp; mv tmp finalLoci.list

#/ determine blocks of loci (for each chromosome or scaffold, create list of [nlocus] loci)
#/ after this step it's possible that blocks at the end of chromosomes/scaffolds have a reduced number of loci so you might want to merge them with the previous block
for scaff in `cut -f1 ../$scalen`; do
    grep ${scaff} finalLoci.list > ${scaff}.tmp
    split -l ${nlocus} --numeric-suffixes=1 $scaff.tmp 5.0.blocks/${scaff}.block\_;
    rm ${scaff}.tmp
done

#/ Concatenate alignments in the same block
for block in `ls 5.0.blocks/* | sed 's,5.0.blocks/,,g'  `; do
    echo $block;
    cat 5.0.blocks/${block} | xargs -n 1 sh -c 'cat 4.bppFormat/$0' >> 5.1.blocksAlign/${block}.seqs.txt;
done


## Prepare control files ***********************************************************************************************
rm BPP_A01.LociSummary.txt
ls 5.0.blocks/* | sed 's,5.0.blocks/,,g' | xargs -n 1 -P 1 sh -c '
    echo $0;
    mkdir 5.2.bpp/$0.rep1
    mkdir 5.2.bpp/$0.rep2
    mkdir 5.2.bpp/$0.rep3
    nlocus=`grep -v Hmel 5.1.blocksAlign/$0.seqs.txt | wc -l`
    # print number of loci
    printf "%s\t%s\n" $0 $nlocus >> BPP_A01.LociSummary.txt;
    # prepare control files 
    cp ../support/bpp-a01.ctl 5.2.bpp/$0.rep1/$0.rep1.ctl
    cp ../support/bpp-a01.ctl 5.2.bpp/$0.rep2/$0.rep2.ctl
    cp ../support/bpp-a01.ctl 5.2.bpp/$0.rep3/$0.rep3.ctl
    # specify number of loci
    sed -i "s,numberLoci,$nlocus,g" 5.2.bpp/$0.rep1/$0.rep1.ctl
    sed -i "s,numberLoci,$nlocus,g" 5.2.bpp/$0.rep2/$0.rep2.ctl
    sed -i "s,numberLoci,$nlocus,g" 5.2.bpp/$0.rep3/$0.rep3.ctl
    # specify alignments file
    sed -i "s,seqs.txt,5.1.blocksAlign/$0.seqs.txt,g" 5.2.bpp/$0.rep1/$0.rep1.ctl
    sed -i "s,seqs.txt,5.1.blocksAlign/$0.seqs.txt,g" 5.2.bpp/$0.rep2/$0.rep2.ctl
    sed -i "s,seqs.txt,5.1.blocksAlign/$0.seqs.txt,g" 5.2.bpp/$0.rep3/$0.rep3.ctl
    # specify output file
    sed -i "s,outp.txt,5.2.bpp/$0.rep1/$0.rep1.out.txt,g" 5.2.bpp/$0.rep1/$0.rep1.ctl
    sed -i "s,outp.txt,5.2.bpp/$0.rep2/$0.rep2.out.txt,g" 5.2.bpp/$0.rep2/$0.rep2.ctl
    sed -i "s,outp.txt,5.2.bpp/$0.rep3/$0.rep3.out.txt,g" 5.2.bpp/$0.rep3/$0.rep3.ctl
    # specify mcmc file
    sed -i "s,mcmc.txt,5.2.bpp/$0.rep1//$0.rep1.mcmc.list,g" 5.2.bpp/$0.rep1/$0.rep1.ctl
    sed -i "s,mcmc.txt,5.2.bpp/$0.rep2//$0.rep2.mcmc.list,g" 5.2.bpp/$0.rep2/$0.rep2.ctl
    sed -i "s,mcmc.txt,5.2.bpp/$0.rep3//$0.rep3.mcmc.list,g" 5.2.bpp/$0.rep3/$0.rep3.ctl
    # specify popmap file
    sed -i "s,imap.txt,../examples/${sgroup}.popmap.txt,g" 5.2.bpp/$0.rep1/$0.rep1.ctl
    sed -i "s,imap.txt,../examples/${sgroup}.popmap.txt,g" 5.2.bpp/$0.rep2/$0.rep2.ctl
    sed -i "s,imap.txt,../examples/${sgroup}.popmap.txt,g" 5.2.bpp/$0.rep3/$0.rep3.ctl
    # replace species line
    sed -i "s,sppline,${sline},g" 5.2.bpp/$0.rep1/$0.rep1.ctl
    sed -i "s,sppline,${sline},g" 5.2.bpp/$0.rep2/$0.rep2.ctl
    sed -i "s,sppline,${sline},g" 5.2.bpp/$0.rep3/$0.rep3.ctl
    # replace individual count line
    sed -i "s,indline,${iline},g" 5.2.bpp/$0.rep1/$0.rep1.ctl
    sed -i "s,indline,${iline},g" 5.2.bpp/$0.rep2/$0.rep2.ctl
    sed -i "s,indline,${iline},g" 5.2.bpp/$0.rep3/$0.rep3.ctl
    # replace tree line
    sed -i "s/treline/${tline1}/g" 5.2.bpp/$0.rep1/$0.rep1.ctl
    sed -i "s/treline/${tline2}/g" 5.2.bpp/$0.rep2/$0.rep2.ctl
    sed -i "s/treline/${tline3}/g" 5.2.bpp/$0.rep3/$0.rep3.ctl
'


## BPP *****************************************************************************************************************

#/ launch BPP
ls 5.1.blocksAlign/* | sed 's,5.1.blocksAlign/,,g' | sed 's,.seqs.txt,,g' | head -n 1 | xargs -n 1 -P 1 sh -c '
    echo $0;
    sbatch ../bpp_Launch.slurm $0 rep1
    # sbatch ../bpp_Launch.slurm $0 rep2
    # sbatch ../bpp_Launch.slurm $0 rep3
'
