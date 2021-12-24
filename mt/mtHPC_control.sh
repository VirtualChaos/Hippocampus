#!/bin/bash 

#Submit in session directory
JID1=$(qsub ~/matlab/Github/Hippocampus/mt/mtsessHPC_submit_file.pbs)

JID2=$(qsub -W depend=afterok:$JID1 ~/matlab/Github/Hippocampus/mt/mtcellHPC_control_file.pbs)
