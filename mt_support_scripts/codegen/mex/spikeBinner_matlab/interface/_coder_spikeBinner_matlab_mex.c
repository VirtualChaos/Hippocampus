/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_spikeBinner_matlab_mex.c
 *
 * Code generation for function '_coder_spikeBinner_matlab_mex'
 *
 */

/* Include files */
#include "_coder_spikeBinner_matlab_mex.h"
#include "_coder_spikeBinner_matlab_api.h"
#include "rt_nonfinite.h"
#include "spikeBinner_matlab_data.h"
#include "spikeBinner_matlab_initialize.h"
#include "spikeBinner_matlab_terminate.h"

/* Custom Source Code */
#include "spikeBinner.h"

/* Function Definitions */
void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs,
                 const mxArray *prhs[])
{
  mexAtExit(&spikeBinner_matlab_atexit);
  /* Module initialization. */
  spikeBinner_matlab_initialize();
  /* Dispatch the entry-point. */
  spikeBinner_matlab_mexFunction(nlhs, plhs, nrhs, prhs);
  /* Module termination. */
  spikeBinner_matlab_terminate();
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLSR2021a(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1,
                           NULL);
  return emlrtRootTLSGlobal;
}

void spikeBinner_matlab_mexFunction(int32_T nlhs, mxArray *plhs[2],
                                    int32_T nrhs, const mxArray *prhs[6])
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  const mxArray *outputs[2];
  int32_T b_nlhs;
  st.tls = emlrtRootTLSGlobal;
  /* Check for proper number of arguments. */
  if (nrhs != 6) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 6, 4,
                        18, "spikeBinner_matlab");
  }
  if (nlhs > 2) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 18,
                        "spikeBinner_matlab");
  }
  /* Call the function. */
  spikeBinner_matlab_api(prhs, nlhs, outputs);
  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    b_nlhs = 1;
  } else {
    b_nlhs = nlhs;
  }
  emlrtReturnArrays(b_nlhs, &plhs[0], &outputs[0]);
}

/* End of code generation (_coder_spikeBinner_matlab_mex.c) */
