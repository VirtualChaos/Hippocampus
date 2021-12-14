#!/bin/bash

matlab -batch -r "addpath ~/matlab/Github/Hippocampus; addpath ~/matlab/Github/Hippocampus/mt; addpath(genpath('~/matlab/Github/DPV')); addpath ~/matlab/Github/npy-matlab-master/npy-matlab; mtsess('auto','redo','save'); cd cells; mt_createCellsList; exit"
