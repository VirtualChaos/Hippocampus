#!/bin/bash

#Submit in session directory

matlab -batch "addpath ~/matlab/Github/Hippocampus; addpath ~/matlab/Github/Hippocampus/mt; addpath(genpath('~/matlab/Github/DPV')) ; addpath ~/matlab/Github/npy-matlab/npy-matlab; addpath ~/matlab/Github/newNpt; mtcell('auto','redo','save'); exit;"
