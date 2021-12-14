#!/bin/bash

var1="auto"
var2="redo"
var3="save"

cmd1="matlab2019a -nojvm -nodisplay -nosplash -r \"addpath ~/matlab/Github/Hippocampus; addpath ~/matlab/Github/Hippocampus/mt; addpath ~/matlab/Github/DPV/miscellaneous; addpath ~/matlab/Github/npy-matlab-master/npy-matlab; mtsess($var1,$var2,$var3); cd cells; mt_createCellsList; exit\""
eval $cmd1
