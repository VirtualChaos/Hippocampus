/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * spikeBinner_matlab.h
 *
 * Code generation for function 'spikeBinner_matlab'
 *
 */

#pragma once

/* Include files */
#include "rtwtypes.h"
#include "spikeBinner_matlab_types.h"
#include "emlrt.h"
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Function Declarations */
void spikeBinner_matlab(const emlrtStack *sp,
                        const emxArray_real_T *flat_spiketimes,
                        const emxArray_real_T *stc,
                        const emxArray_real_T *conditions, real_T repeat,
                        emxArray_real_T *consol_arr, emxArray_real_T *binArr);

/* End of code generation (spikeBinner_matlab.h) */
