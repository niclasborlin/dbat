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

//#define DEBUG(x) mexPrintf x; mexEvalString("pause");
#define DEBUG(x) 

#define IsNonZero(d) ((d)!=0.0)

/*
 * Extract elements of n-by-n diagonal blocks of the real sparse array
 * A and store in the sparse array C. A is square with an integer
 * number of diagonal blocks. C has the same size as A and has space
 * for nzMax elements.
 *
 * Return the number of copied elements, or nzMax+1 if error.
 */
mwIndex extractDiagBlksSparse(const mxArray *A,mxArray *C,mwSize n,
			mwSize nzMax)
{
	mwSize An=mxGetN(A);

	mwIndex *Ar=mxGetIr(A); // Row index vector
	mwIndex *Ac=mxGetJc(A); // Column index vector
	mxDouble *Av=mxGetDoubles(A); // Value vector

	DEBUG(("An = %ld, Ar = %p, Ac = %p, Av = %p.\n",An,Ar,Ac,Av));

	mwIndex *Cr=mxGetIr(C); // Row index vector
	mwIndex *Cc=mxGetJc(C); // Column index vector
	mxDouble *Cv=mxGetDoubles(C); // Value vector

	DEBUG(("Cr = %p, Cc = %p, Cv = %p.\n",Cr,Cc,Cv));

	mwIndex blk; // Column index of first block element.
	mwIndex k=0; // Index into sparse output vectors.
	mwIndex ix,j; // Indices within block.

	// For each block
	for (blk=0; blk<An; blk+=n) {
		DEBUG(("blk = %ld\n",blk));
		// For each column within the block.
		for (j=0; j<n; j++) {
			// Store output pointer into row,value vectors of first
			// non-zero element of this column.
			Cc[blk+j]=k;
			DEBUG(("j = %ld\n",j));

			// Start of non-zeros in this column.
			mwIndex cix=Ac[blk+j];

			// Number of non-zeros in this column.
			mwSize nzElems=Ac[blk+j+1]-Ac[blk+j];
			DEBUG(("cix = %ld, nzElems = %ld\n",cix, nzElems));

			// For each non-zero element in this column...
			for (ix=0; ix<nzElems; ix++) {
				DEBUG(("ix = %ld\n",ix));
				// Actual row index
				mwIndex i=Ar[cix+ix];
				DEBUG(("i = %ld\n",i));
				if (i<blk) {
						// Above block, ignore and continue.
						continue;
				}
				if (i>=blk+n) {
						// Below block, we're done in this column.
						break;
				}
				// Inside block, extract value and check if it is non-zero.
				mxDouble v=Av[cix+ix];

				if (IsNonZero(v)) {
					// Safety check to signal error if value vector is
					// already full.
					if (k>=nzMax) {
						// Close column indices and signal error.
						Cc[An]=k;
						return nzMax+1;
					}

					// Store row index and value.
					Cr[k]=i;
					Cv[k]=v;
					// Advance in row, value vectors.
					k++;
				}
			}
		}
	}
	// Close column indices.
	Cc[An]=k;

	return k;
}

/*
 * Extract elements of n-by-n diagonal blocks of the real dense array
 * A and store in the sparse array C. A is square with an integer
 * number of diagonal blocks. C has the same size as A and has space
 * for nzMax elements.
 *
 * Return the number of copied elements, or nzMax+1 if error.
 */
mwIndex extractDiagBlksDense(const mxArray *A,mxArray *C,mwSize n,
			mwSize nzMax)
{
	mwSize An=mxGetN(A);

	mxDouble *Av=mxGetDoubles(A); // Value vector

	DEBUG(("An = %ld, Av = %p.\n",An,Av));

	mwIndex *Cr=mxGetIr(C); // Row index vector
	mwIndex *Cc=mxGetJc(C); // Column index vector
	mxDouble *Cv=mxGetDoubles(C); // Value vector

	DEBUG(("Cr = %p, Cc = %p, Cv = %p.\n",Cr,Cc,Cv));

	mwIndex k=0; // Index into sparse output vectors.
	mwIndex i,j,ij[2]; // Indices into input array.
	mwIndex blk; // Column index of first block element.

	// For each block...
	for (blk=0; blk<An; blk+=n) {
		DEBUG(("blk = %ld\n",blk));
		// For each column within the block...
		for (j=0; j<n; j++) {
			// Store pointer into row,value vectors of first non-zero
			// element.
			Cc[blk+j]=k;
			DEBUG(("j = %ld\n",j));
			ij[1]=blk+j;
			// For each row within the block...
			for (i=0; i<n; i++) {
				DEBUG(("i = %ld\n",i));
				ij[0]=blk+i;

				DEBUG(("ij[] = [ %ld, %ld ]\n",ij[0],ij[1]));

				// Compute index and extract value.
				mwIndex ix=mxCalcSingleSubscript(A,2,ij);

				DEBUG(("ix = %ld\n",ix));

				double v=Av[ix];

				// Only store non-zero values.
				if (IsNonZero(v)) {
					// Safety check to signal error if value vector is
					// already full.
					if (k>=nzMax) {
						// Close column indices and signal error.
						Cc[An]=k;
						return nzMax+1;
					}

					// Store row index and value.
					Cr[k]=ij[0];
					Cv[k]=v;
					// Advance in row, value vectors.
					k++;
				}
			}
		}
	}
	// Close column indices.
	Cc[An]=k;

	return k;
}

#define A prhs[0]
#define N prhs[1]
#define C plhs[0]

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
	// Check for proper number of arguments.
	if (nrhs!=2) {
		mexErrMsgIdAndTxt("DBAT:extractdiagblocks:nrhs","Two inputs required.");
	}
	if (nlhs>1) {
		mexErrMsgIdAndTxt("DBAT:extractdiagblocks:nlhs","One output required.");
	}

	// Make sure N is real scalar, double.
	if (!mxIsDouble(N) || mxIsComplex(N) || mxGetNumberOfElements(N)!=1) {
		mexErrMsgIdAndTxt("DBAT:extractdiagblocks:notScalar",
				"N must be a real double scalar.");
	}

	// Make sure N is finite, positive integer.
	double nVal=mxGetScalar(N);

	if (!mxIsFinite(nVal) || nVal<=0 || nVal!=floor(nVal)) {
		mexErrMsgIdAndTxt("DBAT:extractdiagblocks:notInt",
						"N must be a positive integer.");
	}

	// Extract block size
	mwSize n=(mwSize)nVal;

	// Make sure A is double
	if (!mxIsDouble(A) || mxIsComplex(A)) {
		mexErrMsgIdAndTxt("DBAT:extractdiagblocks:notDouble",
						"A must be type double.");
	}

	mwSize Am=mxGetM(A);
	mwSize An=mxGetN(A);

	// Make sure A is square.
	if (Am!=An || mxGetNumberOfDimensions(A) != 2) {
		mexErrMsgIdAndTxt("DBAT:extractdiagblocks:notSquare",
				"A must be 2D square.");
	}

	// Make sure A contains has full blocks.
	if (Am % n != 0) {
		mexErrMsgIdAndTxt("DBAT:extractdiagblocks:badSize",
				"The size of A must be a multiple of N.");
	}

	// Pre-allocate sparse output array.
	mwSize nzmax=n*Am;

	if (mxIsSparse(A) && mxGetNzmax(A)<nzmax) {
			nzmax=mxGetNzmax(A);
	}

	DEBUG(("nzmax = %ld\n",nzmax));

	C=mxCreateSparse(Am,An,nzmax,false);
	mwIndex *Cr=mxGetIr(C); // Row index vector
	mwIndex *Cc=mxGetJc(C); // Column index vector
	mxDouble *Cv=mxGetDoubles(C); // Value vector

	if (C==NULL || Cr==NULL || Cc==NULL || Cv==NULL) {
		mexErrMsgIdAndTxt("DBAT:extractdiagblocks:memErr",
				"Memory allocation error for C.");
	}

	if ( mxIsSparse(A) ) {
		// Sparse stuff
		if (extractDiagBlksSparse(A,C,n,nzmax) > nzmax) {
			mexErrMsgIdAndTxt("DBAT:extractdiagblocks:internal",
					"Too many copied elements.");
		}
	} else {
		// Dense stuff
		if (extractDiagBlksDense(A,C,n,nzmax) > nzmax) {
			mexErrMsgIdAndTxt("DBAT:extractdiagblocks:internal",
					"Too many copied elements.");
		}
	}
}
