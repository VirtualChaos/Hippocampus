/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * spikeBinner_matlab_terminate.c
 *
 * Code generation for function 'spikeBinner_matlab_terminate'
 *
 */

/* Include files */
#include "spikeBinner_matlab_terminate.h"
#include "_coder_spikeBinner_matlab_mex.h"
#include "rt_nonfinite.h"
#include "spikeBinner_matlab_data.h"

/* Custom Source Code */
#include "spikeBinner.h"

/* Function Definitions */
void spikeBinner_matlab_atexit(void)
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

void spikeBinner_matlab_terminate(void)
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (spikeBinner_matlab_terminate.c) */
