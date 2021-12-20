#!/bin/bash

#Submit in session directory
pushd cells

echo "$PWD"

cmd1="qsub ~/matlab/Github/Hippocampus/mt/mtcellHPC_submit_file.pbs"

while IFS= read -r cell_no; do
echo "$cell_no"
pushd $cell_no
eval $cmd1
popd
done < cells_list.txt
