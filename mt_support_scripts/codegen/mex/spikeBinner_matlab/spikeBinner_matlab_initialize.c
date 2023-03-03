/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * spikeBinner_matlab_initialize.c
 *
 * Code generation for function 'spikeBinner_matlab_initialize'
 *
 */

/* Include files */
#include "spikeBinner_matlab_initialize.h"
#include "_coder_spikeBinner_matlab_mex.h"
#include "rt_nonfinite.h"
#include "spikeBinner_matlab_data.h"

/* Custom Source Code */
#include "spikeBinner.h"

/* Variable Definitions */
static const volatile char_T *emlrtBreakCheckR2012bFlagVar = NULL;

/* Function Definitions */
void spikeBinner_matlab_initialize(void)
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  mex_InitInfAndNan();
  mexFunctionCreateRootTLS();
  emlrtBreakCheckR2012bFlagVar = emlrtGetBreakCheckFlagAddressR2012b();
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, NULL);
  emlrtEnterRtStackR2012b(&st);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (spikeBinner_matlab_initialize.c) */
