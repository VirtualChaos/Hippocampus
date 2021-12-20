#!/bin/bash

# No. of cell directories with mtcell.mat
find . -name "mtcell.mat" | wc -l

# List of cell directories missing mtcell.mat
find . -name "cell*" | grep -v -e cells_list | sort | cut -d "/" -f 2 > tmp1.txt
find . -name "mtcell.mat" | sort | cut -d "/" -f 2 > tmp2.txt
comm -23 tmp1.txt tmp2.txt

#cwd=`pwd`; for i in `comm -23 tmp1.txt tmp2.txt`; do echo $i; cd $i; qsub ~/matlab/Github/Hippocampus/mt/mtcellHPC_submit_file.pbs; cd $cwd; done

