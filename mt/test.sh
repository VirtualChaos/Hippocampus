#!/bin/bash

matlab -batch "addpath ~/matlab/Github/Hippocampus; addpath ~/matlab/Github/Hippocampus/mt; addpath(genpath('~/matlab/Github/DPV')); addpath ~/matlab/Github/npy-matlab-master/npy-matlab; mtsess('auto','redo','save'); ori = pwd; cd strcat(ori,'/','cells'); mt_createCellsList; exit"
