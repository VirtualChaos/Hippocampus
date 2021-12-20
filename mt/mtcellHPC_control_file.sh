#!/bin/bash 

#Submit in session directory
cmd1="cd cells"
eval=$cmd1

# filename="cells_list.txt"
# n=1
# while read line; do
# echo "Line No. $n : $line"
# n=$(n+1))
# done

cmd2="qsub ~/matlab/Github/Hippocampus/mt/mtcellHPC_submit_file.pbs"

while read cell_no; do
echo "$cell_no"
cmd1="cd $cell_no"
eval $cmd1
eval $cmd2
cmd1="cd .."
eval $cmd1
done < "cells_list.txt"
