#!/bin/bash 

#Submit in session directory
JID1=$(qsub mtsessHPC_submit_file.pbs)

JID2=$(qsub -W depend=afterok:$JID1 mtcellHPC_control_file.sh)
