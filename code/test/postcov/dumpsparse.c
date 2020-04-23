/*
 * EXTRACTDIAGBLOCKS - Extract diagonal blocks.
 *
 * C=EXTRACTDIAGBLOCKS(A,N), returns a sparse matrix C with the
 * N-by-N diagonal blocks of A. A must be real double with an integer
 * number of diagonal blocks. A can be sparse or full.
 *
 */

#include <math.h>
#include "mex.h"
#include "matrix.h"

#define A prhs[0]

void spprint(const mxArray *mx)
{
    mwSize n;
    mwIndex *ir, *jc;
    mwIndex j, k;
    double *pr;
    if( !mxIsSparse(mx) ) return;
    n = mxGetN(mx);
    pr = mxGetDoubles(mx);
    ir = mxGetIr(mx);
    jc = mxGetJc(mx);

    mexPrintf("colPtr = [ ");
    for (j=0; j<=n; j++) {
	mexPrintf("%ld",jc[j]);
	if (j<n) {
	    mexPrintf(",");
	}
	mexPrintf(" ");
    }
    mexPrintf("]\n");

    mexPrintf("rowPtr = [ ");
    for (k=0; k<jc[n]; k++) {
	mexPrintf("%ld",ir[k]);
	if (k<jc[n]-1) {
	    mexPrintf(",");
	}
	mexPrintf(" ");
    }
    mexPrintf("]\n");

    mexPrintf("val    = [ ");
    for (k=0; k<jc[n]; k++) {
	mexPrintf("%g",pr[k]);
	if (k<jc[n]-1) {
	    mexPrintf(",");
	}
	mexPrintf(" ");
    }
    mexPrintf("]\n");
}

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    // Check for proper number of arguments.
    if (nrhs!=1) {
	mexErrMsgTxt("One input required.");
    }
    if (nlhs!=0) {
	mexErrMsgTxt("Zero outputs required.");
    }

    spprint(A);
}
