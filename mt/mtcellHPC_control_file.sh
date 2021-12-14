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

while read cell_no; do
echo "$cell_no"
cmd1="cd "
eval $cmd1
"qsub "
done < "cells_list.txt"
