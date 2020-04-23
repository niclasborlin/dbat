/*
 * DIAGBLKOUTER - Perform an outer product, keeping only diagonal blocks.
 *
 * C=DIAGBLKOUTER(A,B,N) computes the N-by-N diagonal blocks of the
 * product B'*A*B, where the array A is real double, full and
 * symmetric. The array B is real of conforming size. The number of
 * columns of B is an integer multiple of N. The return matrix C is
 * sparse with N-by-N diagonal blocks.
 *
 * No test is performed as to whether A is symmetric.
 */

#include <math.h>
#include "mex.h"
#include "matrix.h"

//#define DEBUG(x) mexPrintf x; // mexEvalString("pause");
#define DEBUG(x)

#define IsNonZero(d) ((d)!=0.0)
#define IsZero(d) ((d)==0.0)

/* Multiplication

// 1-by-1

                   [ a11 a12 a13 ] [ b11 ]
   [ d11 d12 d13 ] [ a21 a22 a23 ] [ b21 ] =
                   [ a31 a32 a33 ] [ b31 ]

                   [ a11 b11 + a12 b21 + a13 b31 ]
   [ d11 d12 d13 ] [ a21 b11 + a22 b21 + a23 b31 ] =
                   [ a31 b11 + a32 b32 + a33 b31 ]

   d11 ( a11 b11 + a12 b21 + a13 b31 ) +
       d12 ( a21 b11 + a22 b21 + a23 b31 ) +
       d13 ( a31 b11 + a32 b32 + a33 b31 ) =

   d11 a11 b11 + d11 a12 b21 + d11 a13 b31 +
       d12 a21 b11 + d12 a22 b21 + d12 a23 b31 +
       d13 a31 b11 + d13 a32 b32 + d13 a33 b31 =

   a11 b11 d11 + a12 b21 d11 + a13 b31 d11 +
       a21 b11 d12 + a22 b21 d12 + a23 b31 d12 +
       a31 b11 d13 + a32 b32 d13 + a33 b31 d13 =


   c11 = sum (i=1..3) (j=1..3) aij bj1 d1i =
         sum (i=1..3) (j=1..3) aij bj1 bi1

// 2-by-2

                   [ a11 a12 a13 ] [ b11 b12 ]
   [ d11 d12 d13 ] [ a21 a22 a23 ] [ b21 b22 ] =
   [ d21 d22 d23 ] [ a31 a32 a33 ] [ b31 b23 ]

   k=1..2, l=1..2,
      ckl = sum (i=1..3) (j=1..3) aij bjl dki =
            sum (i=1..3) (j=1..3) aij bjl bik

// B is M-by-N, A is M-by-M

   k=1..M, l=1..N,
      ckl = sum (i=1..M) (j=1..M) aij bjl dki =
            sum (i=1..M) (j=1..M) aij bjl bik

*/

// Macro to compute linear index in buffer.
#define BUFIX(i,j,n) ((n)*(j)+(i))

// Compute linear index for an array
mwIndex lix(const mxArray *A,mwIndex i,mwIndex j)
{
    // Indices for mxCalcSingleSubscript
    mwIndex ij[2];

    ij[0]=i;
    ij[1]=j;
    return mxCalcSingleSubscript(A,2,ij);
}

/*
 * Compute n-by-n diagonal blocks of the product B'*A*B, where A and B
 * are full and of conforming sizes. A is assumed to be symmetric. C
 * has the appropriate size and has space for nzMax elements. The
 * buffer buf is preallocated to hold n^2 elements.
 *
 * Return the number of copied elements, or nzMax+1 if error.
 */

mwSize outerMultDense(const mxArray *A,const mxArray *B,mxArray *C,
		      mwSize n,mwSize nzMax,double *buf)
{
    mwSize An=mxGetN(A);
    mwSize Bn=mxGetN(B);

    mxDouble *Av=mxGetDoubles(A); // Value vector in A.
    mxDouble *Bv=mxGetDoubles(B); // Value vector in B.

    DEBUG(("An = %ld, Av = %p.\n",An,Av));
    DEBUG(("Bn = %ld, Bv = %p.\n",Bn,Bv));

    mwIndex *Cr=mxGetIr(C); // Row index vector
    mwIndex *Cc=mxGetJc(C); // Column index vector
    mxDouble *Cv=mxGetDoubles(C); // Value vector

    DEBUG(("Cr = %p, Cc = %p, Cv = %p.\n",Cr,Cc,Cv));

    mwIndex kk=0; // Index into sparse output vectors.

    // For each block. blk is first column in block.
    for (mwIndex blk=0; blk<Bn; blk+=n) {
	mexPrintf("pt=%d\n",blk/n);
	DEBUG(("blk = %ld\n",blk));

	// First compute the whole block in local buffer.

	// Result is symmetric. Thus, only subdiagonal part need to be
	// computed.

	// For each column within the block.
	for (mwIndex jj=0; jj<n; jj++) {
	    // Column index in C.
	    mwIndex l=blk+jj;

	    // For each subdiagonal row within the block.
	    for (mwIndex ii=jj; ii<n; ii++) {
		// Row index in C.
		mwIndex k=blk+ii;

		DEBUG(("ii,jj,k,l = %ld,%ld,%ld,%ld\n",ii,jj,k,l));

		double sum=0;
		// i is index of A(i,j) and B(i,k)
		for (mwIndex i=0; i<An; i++) {
		    double bik=Bv[lix(B,i,k)];

		    DEBUG(("ik = [%ld, %ld], bik=%g\n",i,k,bik));
		    if (IsZero(bik)) {
			continue;
		    }

		    // j is index of A(i,j) and B(j,l)
		    for (mwIndex j=0; j<An; j++) {
			double bjl=Bv[lix(B,j,l)];

			DEBUG(("jl = [%ld, %ld], bjl=%g\n",j,l,bjl));

			double aij=Av[lix(A,i,j)];

			DEBUG(("ij = [%ld, %ld], aij=%g\n",i,j,aij));

			sum+=aij*bjl*bik;
		    }
		}

		// Store computed value in buffer.
		buf[BUFIX(ii,jj,n)]=sum;

		DEBUG(("BUF ii,jj = [%ld, %ld], ix=%ld, val=%g\n",ii,jj,
		       BUFIX(ii,jj,n),sum));

		if (ii!=jj) {
		    buf[BUFIX(jj,ii,n)]=sum;
		    DEBUG(("BUF jj,ii = [%ld, %ld], v=%g\n",jj,ii,sum));
		}
	    }
	}

	// Copy all buffer elements to output

	// For each column within the block.
	for (mwIndex jj=0; jj<n; jj++) {
	    // Column of C.
	    mwIndex l=blk+jj;

	    // Store pointer into row,value vectors of first non-zero
	    // element.
	    Cc[l]=kk;

	    // For each row within the block.
	    for (mwIndex ii=0; ii<n; ii++) {
		// Row of C.
		mwIndex k=blk+ii;

		double v=buf[BUFIX(ii,jj,n)];

		// Only store non-zero values.
		if (IsNonZero(v)) {
		    // Safety check to signal error if value vector is
		    // already full.
		    if (kk>=nzMax) {
			// Close column indices and signal error.
			Cc[Bn]=kk;
			return nzMax+1;
		    }

		    // Store row index and value.
		    Cr[kk]=k;
		    Cv[kk]=v;
		    // Advance in row, value vectors.
		    kk++;
		}
	    }
	}
    }

    // Close column indices.
    Cc[Bn]=kk;

    return kk;
}

/*
 * Compute n-by-n diagonal blocks of the product B'*A*B, where A is
 * symmetric and full and B is sparse and of conforming size. C has
 * the appropriate size and has space for nzMax elements. The buffer
 * buf is preallocated to hold n^2 elements.
 *
 * Return the number of copied elements, or nzMax+1 if error.
 */

mwSize outerMultSparse(const mxArray *A,const mxArray *B,
		       mxArray *C,mwSize n,mwSize nzMax,double *buf)
{
    mwSize Bn=mxGetN(B);

    mxDouble *Av=mxGetDoubles(A); // Value vector in A.

    DEBUG(("Av = %p.\n",Av));

    mwIndex *Br=mxGetIr(B); // Row index vector
    mwIndex *Bc=mxGetJc(B); // Column index vector
    mxDouble *Bv=mxGetDoubles(B); // Value vector in B.

    DEBUG(("Br = %p, Bc = %p, Bv = %p.\n",Br,Bc,Bv));

    mwIndex *Cr=mxGetIr(C); // Row index vector
    mwIndex *Cc=mxGetJc(C); // Column index vector
    mxDouble *Cv=mxGetDoubles(C); // Value vector

    DEBUG(("Cr = %p, Cc = %p, Cv = %p.\n",Cr,Cc,Cv));

    mwIndex kk=0; // Index into sparse output vectors.

    // For each block. blk is first column in block.
    for (mwIndex blk=0; blk<Bn; blk+=n) {
	//mexPrintf("pt=%d\n",blk/n);
	DEBUG(("blk = %ld\n",blk));

	// First compute the whole block in local buffer.

	// Result is symmetric. Thus, only subdiagonal part need to be
	// computed.

	// For each column within the block.
	for (mwIndex jj=0; jj<n; jj++) {
	    // Column index in C.
	    mwIndex l=blk+jj;

	    // For each subdiagonal row within the block.
	    for (mwIndex ii=jj; ii<n; ii++) {
		// Row index in C.
		mwIndex k=blk+ii;

		DEBUG(("ii,jj,k,l = %ld,%ld,%ld,%ld\n",ii,jj,k,l));

		double sum=0;

		// i is index of A(i,j) and B(i,k)
		// Run down column B(:,k)

		// Non-zero elements in B(:,k)
		mwSize nnzBik=Bc[k+1]-Bc[k];
		// Index of first non-zero in B(:,k)
		mwIndex bikTop=Bc[k];
		// For each non-zero element of B(:,k)
		for (mwIndex bikIx=0; bikIx<nnzBik; bikIx++) {
		    mwIndex i=Br[bikTop+bikIx];
		    double bik=Bv[bikTop+bikIx];

		    DEBUG(("ik = [%ld, %ld], bik=%g\n",i,k,bik));

		    // Non-zero elements in B(:,l)
		    mwSize nnzBjl=Bc[l+1]-Bc[l];
		    // Index of first non-zero in B(:,l)
		    mwIndex bjlTop=Bc[l];
		    // For each non-zero element of B(:,l)
		    for (mwIndex bjlIx=0; bjlIx<nnzBjl; bjlIx++) {
			mwIndex j=Br[bjlTop+bjlIx];
			double bjl=Bv[bjlTop+bjlIx];

			DEBUG(("jl = [%ld, %ld], bjl=%g\n",j,l,bjl));

			double aij=Av[lix(A,i,j)];

			DEBUG(("ij = [%ld, %ld], aij=%g\n",i,j,aij));

			sum+=aij*bjl*bik;
		    }
		}

		// Store computed value in buffer.
		buf[BUFIX(ii,jj,n)]=sum;

		DEBUG(("BUF ii,jj = [%ld, %ld], ix=%ld, val=%g\n",ii,jj,
		       BUFIX(ii,jj,n),sum));

		if (ii!=jj) {
		    buf[BUFIX(jj,ii,n)]=sum;
		    DEBUG(("BUF jj,ii = [%ld, %ld], v=%g\n",jj,ii,sum));
		}
	    }
	}

	// Copy all buffer elements to output

	// For each column within the block.
	for (mwIndex jj=0; jj<n; jj++) {
	    // Column of C.
	    mwIndex l=blk+jj;

	    // Store pointer into row,value vectors of first non-zero
	    // element.
	    Cc[l]=kk;

	    // For each row within the block.
	    for (mwIndex ii=0; ii<n; ii++) {
		// Row of C.
		mwIndex k=blk+ii;

		double v=buf[BUFIX(ii,jj,n)];

		// Only store non-zero values.
		if (IsNonZero(v)) {
		    // Safety check to signal error if value vector is
		    // already full.
		    if (kk>=nzMax) {
			// Close column indices and signal error.
			Cc[Bn]=kk;
			return nzMax+1;
		    }

		    // Store row index and value.
		    Cr[kk]=k;
		    Cv[kk]=v;
		    // Advance in row, value vectors.
		    kk++;
		}
	    }
	}
    }

    // Close column indices.
    Cc[Bn]=kk;

    return kk;
}

#define C plhs[0]

// The gateway function
void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
    // Check for proper number of arguments.
    if (nrhs!=3) {
	mexErrMsgIdAndTxt("DBAT:diagblkouter:nrhs","Three inputs required.");
    }
    if (nlhs>1) {
	mexErrMsgIdAndTxt("DBAT:diagblkouter:nlhs","One output required.");
    }

    const mxArray *A=prhs[0];
    const mxArray *B=prhs[1];
    const mxArray *N=prhs[2];

    // Make sure N is real scalar, double.
    if (!mxIsDouble(N) || mxIsComplex(N) || mxGetNumberOfElements(N)!=1) {
	mexErrMsgIdAndTxt("DBAT:diagblkouter:notScalar",
			  "N must be a real double scalar.");
    }

    // Make sure N is finite, positive integer.
    double nVal=mxGetScalar(N);

    if (!mxIsFinite(nVal) || nVal<=0 || nVal!=floor(nVal)) {
	mexErrMsgIdAndTxt("DBAT:diagblkouter:notInt",
			  "N must be a positive integer.");
    }

    // Extract block size
    mwSize n=(mwSize)nVal;

    // Make sure A is double real
    if (!mxIsDouble(A) || mxIsComplex(A)) {
	mexErrMsgIdAndTxt("DBAT:diagblkouter:notDouble",
			  "A must be type real double.");
    }

    mwSize Am=mxGetM(A);
    mwSize An=mxGetN(A);

    // Make sure A is square.
    if (Am!=An || mxGetNumberOfDimensions(A) != 2) {
	mexErrMsgIdAndTxt("DBAT:diagblkouter:notSquare",
			  "A must be 2D square.");
    }

    // Make sure A is full.
    if ( mxIsSparse(A) ) {
	mexErrMsgIdAndTxt("DBAT:diagblkouter:notFull",
			  "A must be full.");
    }

    // Make sure B is double real
    if (!mxIsDouble(B) || mxIsComplex(B)) {
	mexErrMsgIdAndTxt("DBAT:diagblkouter:notDouble",
			  "B must be type real double.");
    }

    mwSize Bm=mxGetM(B);
    mwSize Bn=mxGetN(B);

    // Make sure A, B has conforming dimensions
    if (An!=Bm || mxGetNumberOfDimensions(B) != 2) {
	mexErrMsgIdAndTxt("DBAT:diagblkouter:notSquare",
			  "Incorrect dimensions.");
    }

    // Make sure B correspond to full blocks
    if (Bn % n != 0) {
	mexErrMsgIdAndTxt("DBAT:diagblkouter:badSize",
			  "The column size of B must be a multiple of N.");
    }

    // Pre-allocate sparse output array.
    mwSize nzmax=n*Bn;

    DEBUG(("nzmax = %ld\n",nzmax));

    C=mxCreateSparse(Bn,Bn,nzmax,false);
    mwIndex *Cr=mxGetIr(C); // Row index vector
    mwIndex *Cc=mxGetJc(C); // Column index vector
    mxDouble *Cv=mxGetDoubles(C); // Value vector

    if (C==NULL || Cr==NULL || Cc==NULL || Cv==NULL) {
	mexErrMsgIdAndTxt("DBAT:diagblkouter:memErr",
			  "Memory allocation error for C.");
    }

    // Local buffer
    double *buf=mxCalloc(n*n,sizeof(*buf));

    if (buf==NULL) {
	mexErrMsgIdAndTxt("DBAT:diagblkouter:outOfMemory",
			  "Out of memory for local buffer.");
    }

    if ( mxIsSparse(B) ) {
	// Sparse products
	if (outerMultSparse(A,B,C,n,nzmax,buf) > nzmax) {
	    mexErrMsgIdAndTxt("DBAT:diagblkouter:internal",
			      "Too many copied elements.");
	}
    } else {
	// Dense products
	if (outerMultDense(A,B,C,n,nzmax,buf) > nzmax) {
	    mexErrMsgIdAndTxt("DBAT:diagblkouter:internal",
			      "Too many copied elements.");
	}
    }

    if (buf) {
	mxFree(buf);
	buf=NULL;
    }
}
