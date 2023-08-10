## BPP - 

## 1. Prepare bed files
bash ~/code/heliconius_seixas/2.elevatus_pardalinus/2.SpeciesTree/bpp-pipeline/1.bpp.prepbed.sh \
    ~/code/heliconius_seixas/2.elevatus_pardalinus/2.SpeciesTree/bpp-pipeline/

## 2. Prepare alignment files
bash ~/code/heliconius_seixas/2.elevatus_pardalinus/2.SpeciesTree/bpp-pipeline/2.prepAlignments.sh \
    noncod \
    100 250 2000 \
    hmelv25 example1 
    
#/ 3. Run BPP
bash ~/code/heliconius_seixas/2.elevatus_pardalinus/2.SpeciesTree/bpp01/bpp-a01-3.runBPP.slurm \
    "noncod.minL100.minL250.minG2000" \
    "example1" \
    "4 parAnd parAma eleAma eleGui" \
    "2 2 2 2" \
    "((parAnd,parAma),(eleAma,eleGui));"  \
    "(parAnd,(parAma,(eleAma,eleGui)));"  \
    "(parAma,(parAnd,(eleAma,eleGui)));" 
