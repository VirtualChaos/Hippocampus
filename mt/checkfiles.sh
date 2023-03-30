#!/bin/bash

# No. of cell directories with mtcell.mat
#find . -name "mtcell.mat" | wc -l

# List of cell directories missing mtcell.mat
#find . -name "cell*" | grep -v -e cells_list | sort | cut -d "/" -f 2 > tmp1.txt
#find . -name "mtcell.mat" | sort | cut -d "/" -f 2 > tmp2.txt
#comm -23 tmp1.txt tmp2.txt

#cwd=`pwd`; for i in `comm -23 tmp1.txt tmp2.txt`; do echo $i; cd $i; qsub ~/matlab/Github/Hippocampus/mt/mtcellHPC_submit.pbs; cd $cwd; done

# No. of cell directories with mtcell.mat
find . -name "mtcellcombined.mat" | wc -l

# List of cell directories missing mtcell.mat
find . -name "cellBlock*" | grep -v -e cellBlockList | sort | cut -d "/" -f 2 > tmp1.txt
find . -name "mtcellcombined.mat" | sort | cut -d "/" -f 2 > tmp2.txt
comm -23 tmp1.txt tmp2.txt

cwd=`pwd`; for i in `comm -23 tmp1.txt tmp2.txt`; do echo $i; cd $i; qsub ~/matlab/Github/Hippocampus/mt/mtcellcombinedHPC_submit.pbs; cd $cwd; done
