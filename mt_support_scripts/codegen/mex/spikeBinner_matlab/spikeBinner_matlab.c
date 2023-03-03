/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * spikeBinner_matlab.c
 *
 * Code generation for function 'spikeBinner_matlab'
 *
 */

/* Include files */
#include "spikeBinner_matlab.h"
#include "rt_nonfinite.h"
#include "spikeBinner_matlab_types.h"
#include "mwmathutil.h"

/* Custom Source Code */
#include "spikeBinner.h"

/* Function Definitions */
void spikeBinner_matlab(const emlrtStack *sp,
                        const emxArray_real_T *flat_spiketimes,
                        const emxArray_real_T *stc,
                        const emxArray_real_T *conditions, real_T repeat,
                        emxArray_real_T *consol_arr, emxArray_real_T *binArr)
{
  real_T d;
  int32_T i;
  (void)sp;
  /* coder.cinclude('spikeBinner.h'); */
  d = muDoubleScalarRound(repeat);
  if (d < 2.147483648E+9) {
    if (d >= -2.147483648E+9) {
      i = (int32_T)d;
    } else {
      i = MIN_int32_T;
    }
  } else if (d >= 2.147483648E+9) {
    i = MAX_int32_T;
  } else {
    i = 0;
  }
  spikeBinner(&flat_spiketimes->data[0], &stc->data[0], &conditions->data[0], i,
              &consol_arr->data[0], &binArr->data[0], flat_spiketimes->size[0],
              stc->size[0], binArr->size[0]);
}

/* End of code generation (spikeBinner_matlab.c) */
