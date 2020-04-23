/*
 * ICPC_MEX Inverse Cholesky Post Cov computation
 *
 * C=ICPC_MEX(LA,LB,LC,N) computes a subset of the posterior
 * covariance matrix C of the object points from the Cholesky
 * factorisation L=[LA, 0; LB, LC] of the normal matrix. The matrix
 * LA is 3M-by-3M sparse, block diagonal with 3-by-3 lower triangular
 * blocks. The matrix LB is K-by-3M and sparse. The K-by-K matrix LC
 * is full lower triangular. All matrices are real.
 *
 * The returned 3M-by-3M matrix C is symmetric block diagonal with
 * N-by-N blocks, where N is 1 or 3.
 *
 * No tests are performed to check whether LA and LC have the
 * prescribed structure.
 *
 */

#include <string.h>
#include "mex.h"
#include "matrix.h"

#define DEBUG(x) mexPrintf x; // mexEvalString("pause");
//#define DEBUG(x)

#define IsNonZero(d) ((d)!=0.0)
#define IsZero(d) ((d)==0.0)

// Linear index into compact 3-by-3 array h=[a11,a21,a31,a22,a32,a33].
// NB: i and j are one-based. Linear index is zero-based.
#define C3IX(i,j) (((i)-1)+((j)-1)*3-((j)-1+((j-1)/2)))

#define MAX(a,b) ((a) > (b) ? (a) : (b))

// Linear index for element (I,J) into full array of size M-by-N
// stored in column-major format. (I,J) are one-based. Linear index is
// zero-based.
#define LIX(m,n,i,j) (((i)-1) + (m)*((j)-1))

/*
 * Compute inverse of subdiagonal 3-by-3 block starting at element
 * A(i,i). A is assumed to have non-zero elements only in subdiagonal
 * 3-by-3 block.
 * i is zero-based index of first row/column of block.
 *
 * Compact result (diagonal + subdiagonal elements only) is stored in
 * buffer D.
 *
 * Return 0 if ok, -1 otherwise (division by zero).
 */
int invertAblock(const mwIndex *Ar,const mwIndex *Ac,const mxDouble *Av,
		mwIndex i,double *D)
{
	DEBUG(("\n------- invertAblock entry ------\n"));
	/*
	 * | D11 |     |     |   | A11 |     |     |   | 1 |   |   |
	 * | D21 | D22 |     | * | A21 | A22 |     | = |   | 1 |   |
	 * | D31 | D32 | D33 |   | A31 | A32 | A33 |   |   |   | 1 |
	 *
	 * | D11 A11 |   |         |   |         |   | 1 |
	 * | D21 A11 | + | D22 A21 | + |         | = | 0 |
	 * | D31 A11 |   | D32 A21 |   | D33 A31 |   | 0 |
	 *
	 * |   |   |         |   |         |   | 0 |
	 * |   | + | D22 A22 | + |         | = | 1 |
	 * |   |   | D32 A22 |   | D33 A32 |   | 0 |
	 *
	 * |   |   |   |   |         |   | 0 |
	 * |   | + |   | + |         | = | 0 |
	 * |   |   |   |   | D33 A33 |   | 1 |
	 *
	 * D33 = 1/A33
	 * D22 = 1/A22
	 * D32 = -D33 A32 / A22
	 * D11 = 1/A11
	 * D21 = -D22 A21 / A11
	 * D31 = (-D32 A21 - D33 A31) / A11
	 *
	 */

	// Column i+2
	mwIndex ix3=Ac[i+2];
	mwIndex row3=Ar[ix3];
	mxDouble a33=Av[ix3];

	DEBUG(("a33 = %g\n",a33));

	// row index should be equal to column index
	if (row3!=i+2 || IsZero(a33))
		return -1;

	// D33 = 1/A33
	D[C3IX(3,3)]=1/a33;

	DEBUG(("d33 = %g\n",D[C3IX(3,3)]));

	// Column i+1
	mwIndex ix2=Ac[i+1];
	mwIndex row2=Ar[ix2];
	mxDouble a22=Av[ix2];

	DEBUG(("a22 = %g\n",a22));

	// row index should be equal to column index
	if (row2!=i+1 || IsZero(a22))
		return -1;

	// D22 = 1/A22
	D[C3IX(2,2)]=1/a22;

	DEBUG(("d22 = %g\n",D[C3IX(2,2)]));

	// Advance to next row,value
	ix2++;

	// D32 = -D33 A32 / A22
	if (ix2<ix3 && Ar[ix2]==i+2) {
		mxDouble a32=Av[ix2];
		DEBUG(("a32 = %g\n",a32));
		D[C3IX(3,2)]=-D[C3IX(3,3)]*a32/a22;
	} else {
		// Next non-zero element is in next column => A32=0 => D32=0.
		DEBUG(("a32 = %g\n",0.0));
		D[C3IX(3,2)]=0;
	}
	DEBUG(("d32 = %g\n",D[C3IX(3,2)]));

	// Column i
	mwIndex ix1=Ac[i];
	mwIndex row1=Ar[ix1];
	mxDouble a11=Av[ix1];

	DEBUG(("a11 = %g\n",a11));

	// row index should be equal to column index
	if (row1!=i || IsZero(a11))
		return -1;

	// D11 = 1/A11
	D[C3IX(1,1)]=1/a11;

	DEBUG(("d11 = %g\n",D[C3IX(1,1)]));

	// Advance to next row, value
	ix1++;

	mxDouble a21=0,a31=0;
	if (ix1<ix2 && Ar[ix1]==i+1) {
		a21=Av[ix1];
		ix1++;
	}
	DEBUG(("a21 = %g\n",a21));

	// D21 = -D22 A21 / A11
	D[C3IX(2,1)]=-D[C3IX(2,2)]*a21/a11;

	DEBUG(("d21 = %g\n",D[C3IX(2,1)]));

	if (ix1<ix2 && Ar[ix1]==i+2) {
		a31=Av[ix1];
		ix1++;
	}
	DEBUG(("a31 = %g\n",a31));

	// D31 = (-D32 A21 - D33 A31) / A11
	D[C3IX(3,1)]=(-D[C3IX(3,2)]*a21 -D[C3IX(3,3)]*a31)/a11;

	DEBUG(("d31 = %g\n",D[C3IX(3,1)]));

	return 0;
}

/*
 * Compute inner products of the UA = inv(LA) block. UA is stored
 * compactly in D = [D11,D21,D31,D22,D32,D33]. Only compute the
 * diagonal elements, i.e., BUF(i)=sum(D(:,i).^2).
 */
void computeUAinnerProductsBlock1(const mxDouble *D,mxDouble *BUF)
{
	DEBUG(("\n------- computeUAinnerProductsBlock1 entry ------\n"));
	for (int j=1; j<=3; j++) {
		mxDouble sum=0;
		for (int i=j; i<=3; i++) {
			DEBUG(("(i,j,v)=(%ld,%ld,%g)\n",i,j,D[C3IX(i,j)]));
			sum += D[C3IX(i,j)] * D[C3IX(i,j)];
		}
		DEBUG(("BUF[%ld]=%g\n",j-1,sum));
		BUF[j-1]=sum;
	}
}

/*
 * Compute inner products of the UA = inv(LA) block. UA is stored
 * compactly in D = [D11,D21,D31,D22,D32,D33]. Compute all subdiagonal
 * elements, i.e., BUF(i,j)=D(:,i)'*D(:,j) for i>=j. BUF has the same
 * storage pattern as D.
 *
 */
void computeUAinnerProductsBlock3(const mxDouble *D,mxDouble *BUF)
{
	DEBUG(("\n------- computeUAinnerProductsBlock3 entry ------\n"));
	for (int j=1; j<=3; j++) {
		for (int i=j; i<=3; i++) {
			// Element BUF(i,j) is sum(D(:,i).*D(:,j)).
			// D(:,j) has non-zeros for D(i,j), i>=j.
			mxDouble sum=0;
			for (int k=MAX(i,j); k<=3; k++) {
				DEBUG(("(i,j,k,v1,v2)=(%ld,%ld,%ld,%g,%g)\n",
						i,j,k,D[C3IX(k,i)],
						D[C3IX(k,j)]));
				sum += D[C3IX(k,i)] * D[C3IX(k,j)];
			}
			DEBUG(("(i,j)=(%ld,%ld), BUF[%ld]=%g\n",i,j,C3IX(i,j),sum));
			BUF[C3IX(i,j)]=sum;
		}
	}
}

/*
 * Compute one row of one block column of y=LB * UA
 *
 * Bs is stop index vector, i.e., once we've reached these indices,
 * there are no more non-zeros in this column.
 * Br is vector with actual row numbers.
 * Bv contain non-zero values of B.
 * Bix is 3-vector with next index for each column to be considered.
 * Bix is updated on output if an element on this row is referred to.
 * i is current row number, zero-based.
 * j is first column of block, zero-based.
 * UA is buffer with lower triangular part of inv(LA).
 * Y is buffer for row.
 *
 */

void computeRowInBlock(const mwIndex *Br,const mwIndex *Bs,const mxDouble *Bv,
			mwIndex *Bix,mwIndex i,const mxDouble *UA,
			mxDouble *Y)
{
	DEBUG(("\n------- computeRowInBlock entry ------\n"));
	/*
	 * | y11 | = | A11 b11 | + | A21 b12 | + | A31 b13 |
	 * | y12 | = |         |   | A22 b12 | + | A32 b13 |
	 * | y13 | = |         |   |         |   | A33 b13 |
	 */

	// Set all Y values to zero.
	memset(Y,0,3*sizeof(*Y));

	DEBUG(("Bix=[ "));
	for (int i=1; i<=3; i++) {
		DEBUG(("%ld ",Bix[i-1]));
	}
	DEBUG(("]\n"));

	// Check whether we have non-zeros on this row in each column.
	for (int k=1; k<=3; k++) {
		// Bix[k-1] is index into Br and Bv.
		// If Bix[k-1]>=Bc[j+k-1+1], next non-zero is in next column.
		// If Br[Bix[k-1]]!=i, next non-zero is further down in this
		// column.
		DEBUG(("k=%ld, Bix=%ld, Br=%ld, i=%ld, Bs=%ld, ",k,Bix[k-1],
				Br[Bix[k-1]],i,Bs[k-1]));
		if (Br[Bix[k-1]]==i && Bix[k-1]<Bs[k-1]) {
			// We have a non-zero
			mxDouble b=Bv[Bix[k-1]];
			DEBUG(("b=%g\n",b));
			// Advance index
			Bix[k-1]++;
			// Update Y values
			for (int j=1; j<=k; j++) {
				Y[j-1] += UA[C3IX(k,j)] * b;
			}
		} else {
			DEBUG(("b=zero\n"));
		}
	}
}

/*
 * Compute UB = LC \ (LB * inv(LA)) for each block column. Compute and
 * update the requested dot products. The UB block is re-used between
 * block columns.
 *
 * LA is 3M-by-3M sparse, block diagonal with 3-by-3 lower triangular blocks.
 * LB is K-by-3M and sparse.
 * LC is K-by-K and full lower triangular.
 * n is the requested block size, n=1 or n=3.
 * UB is preallocated to hold K-by-3 elements.
 * C is sparse return array, preallocated to hold n non-zero elements per
 * column.
 */

int computeUBandDotProducts(const mxArray *LA,const mxArray *LB,
			const mxArray *LC,int n,mxDouble *UB,
			mxArray *C)
{
	const mwIndex *LAr=mxGetIr(LA); // Row number vector
	const mwIndex *LAc=mxGetJc(LA); // Column index vector
	const mxDouble *LAv=mxGetDoubles(LA); // Value vector

	const mwIndex *LBr=mxGetIr(LB); // Row number vector
	const mwIndex *LBc=mxGetJc(LB); // Column index vector
	const mxDouble *LBv=mxGetDoubles(LB); // Value vector

	mwSize K=mxGetM(LB);
	mwSize L=mxGetN(LB); // =3M

	const mxDouble *LCv=mxGetDoubles(LB); // Value vector

	for (mwIndex blk=0; blk<L; blk+=n) {
		DEBUG(("\n\nblk=%ld\n",blk));

		// Buffer for lower triangular UA=inv(LA) block.
		mxDouble UA[6];

		// Buffer for inner products. Will only use 3 for
		// diagonal-only computations.
		mxDouble IP[6];

		// Compute inv(LA) of one block.
		if (invertAblock(LAr,LAc,LAv,blk,UA)) {
			mexErrMsgIdAndTxt("DBAT:icpc_mex:badInvert",
					"Bad inverse of LA block.");
		}

		DEBUG(("\nUA printout\n"));

		// Print UA block
		for (int i=0; i<6; i++) {
			DEBUG(("UA[%ld]=%g\n",i,UA[i]));
		}

		// Initialize inner product buffer to inner products of UA.

		// Compute diagonal inner products of UA block.
		if (n==1) {
			computeUAinnerProductsBlock1(UA,IP);
		} else {
			computeUAinnerProductsBlock3(UA,IP);
		}

		// Print IP values
		for (int i=0; i<(n==1 ? 3 : 6); i++) {
			DEBUG(("IP[%ld]=%g\n",i,IP[i]));
		}

		// Compute UB for this block.

		for (mwIndex row=1; row<=K; row++) {
			// Indices into each column within block column.
			mwIndex LBix[3];

			// Initialize indices to top of each column.
			memcpy(LBix,LBc+blk,3* sizeof(*LBc));

			// Stop indices are the first in next column.
			const mwIndex *LBs=LBc+blk+1;

			for (int i=0; i<3; i++) {
				DEBUG(("LBs[%ld]=%ld\n",i,LBs[i]));
			}

			// Buffer for row of LBUA.
			mxDouble LBUA[3];

			DEBUG(("\nrow=%ld\n",row));
			computeRowInBlock(LBr,LBs,LBv,LBix,row,UA,LBUA);

			for (int i=0; i<3; i++) {
				DEBUG(("LBUA[%ld]=%g\n",i,LBUA[i]));
			}

			// Subtract L(row,1:row-1)*UB(1:row-1,:) from LBUA.
			for (int k=1; k<=3; k++) {
				for (mwIndex j=1; j<row; j++) {
					LBUA[k-1] -= LCv[LIX(K,K,row,j)] *
						UB[LIX(K,3,j,k)];
				}
				// Divide by Lii
				LBUA[k-1] /= LCv[LIX(K,K,row,row)];
				// Copy to UB
				UB[LIX(K,3,row,k)] = LBUA[k-1];
			}

			// Update inner products.
			if (n==1) {
				for (int k=0; k<3; k++) {
					IP[k] += LBUA[k]*LBUA[k];
				}
			} else {
				for (int j=1; j<=3; j++) {
					for (int i=j; i<=3; i++) {
						IP[C3IX(i,j)] += 
							LBUA[i-1]*LBUA[j-1];
					}
				}
			}
		}

	}
	return 0;
}

#define C plhs[0]

// The gateway function
void mexFunction(int nlhs, mxArray *plhs[],
		int nrhs, const mxArray *prhs[])
{
	// Check for proper number of arguments.
	if (nrhs!=4) {
		mexErrMsgIdAndTxt("DBAT:icpc_mex:nrhs","Four inputs required.");
	}
	if (nlhs>1) {
		mexErrMsgIdAndTxt("DBAT:icpc_mex:nlhs","One output required.");
	}

	const mxArray *LA=prhs[0];
	const mxArray *LB=prhs[1];
	const mxArray *LC=prhs[2];
	const mxArray *N=prhs[3];

	// Make sure N is real scalar, double.
	if (!mxIsDouble(N) || mxIsComplex(N) || mxGetNumberOfElements(N)!=1) {
		mexErrMsgIdAndTxt("DBAT:icpc_mex:notScalar",
				"N must be a real double scalar.");
	}

	// Make sure N is 1 or 3.
	double nVal=mxGetScalar(N);

	if (nVal!=3.0 && nVal!=1.0) {
		mexErrMsgIdAndTxt("DBAT:icpc_mex:badN",
				"N must be 1 or 3.");
	}

	// Store block size
	mwSize n=(mwSize)nVal;

	// Make sure A is double real
	if (!mxIsDouble(LA) || mxIsComplex(LA)) {
		mexErrMsgIdAndTxt("DBAT:icpc_mex:notDouble",
				"A must be type real double.");
	}

	// Make sure LA is square and of size 3M.
	mwSize LAm=mxGetM(LA);
	mwSize LAn=mxGetN(LA);

	if (LAm!=LAn || mxGetNumberOfDimensions(LA) != 2 || LAm % 3!=0) {
		mexErrMsgIdAndTxt("DBAT:icpc_mex:notSquare",
				"LA must be 2D square of size 3M-by-3M.");
	}

	// Make sure LA is full.
	if ( !mxIsSparse(LA) ) {
		mexErrMsgIdAndTxt("DBAT:icpc_mex:notSparse",
				"LA must be sparse.");
	}

	// Make sure LC is square.
	mwSize LCm=mxGetM(LC);
	mwSize LCn=mxGetN(LC);

	if (LCm!=LCn || mxGetNumberOfDimensions(LC) != 2) {
		mexErrMsgIdAndTxt("DBAT:icpc_mex:notSquare",
				"LC must be 2D square of size 3M-by-3M.");
	}

	// Make sure LC is full.
	if ( mxIsSparse(LC) ) {
		mexErrMsgIdAndTxt("DBAT:icpc_mex:notFull",
				"LC must be full.");
	}

	// Make sure LB is double real.
	if (!mxIsDouble(LB) || mxIsComplex(LB)) {
		mexErrMsgIdAndTxt("DBAT:icpc_mex:notDouble",
				"LB must be type real double.");
	}

	mwSize LBm=mxGetM(LB);
	mwSize LBn=mxGetN(LB);

	// Make sure LB has correct size.
	if (LAn!=LBn || LCm!=LBm || mxGetNumberOfDimensions(LB) != 2) {
		mexErrMsgIdAndTxt("DBAT:icpc_mex:notSquare",
				"Incorrect dimensions.");
	}

	// Pre-allocate sparse output array.
	mwSize nzmax=n*LAn;

	DEBUG(("nzmax = %ld\n",nzmax));

	C=mxCreateSparse(LAm,LAm,nzmax,false);
	mwIndex *Cr=mxGetIr(C); // Row number vector
	mwIndex *Cc=mxGetJc(C); // Column index vector
	mxDouble *Cv=mxGetDoubles(C); // Value vector

	if (C==NULL || Cr==NULL || Cc==NULL || Cv==NULL) {
		mexErrMsgIdAndTxt("DBAT:icpc_mex:memErr",
				"Memory allocation error for C.");
	}

	// Buffer for lower triangular UA=inv(LA) block.
	mxDouble UA[6];

	const mwIndex *LAr=mxGetIr(LA); // Row number vector
	const mwIndex *LAc=mxGetJc(LA); // Column index vector
	const mxDouble *LAv=mxGetDoubles(LA); // Value vector

	// Compute inv(LA) of one block.
	if (invertAblock(LAr,LAc,LAv,0,UA)) {
		mexErrMsgIdAndTxt("DBAT:icpc_mex:badInvert",
				"Bad inverse of LA block.");
	}

	DEBUG(("\nUA printout\n"));

	// Print UA block
	for (int i=0; i<6; i++) {
		DEBUG(("UA[%ld]=%g\n",i,UA[i]));
	}

	mwIndex *LBr=mxGetIr(LB); // Row number vector
	mwIndex *LBc=mxGetJc(LB); // Column index vector
	mxDouble *LBv=mxGetDoubles(LB); // Value vector

	// Compute first row of first block in LBUA.
	mwIndex blk=0;

	for (int i=0; i<3; i++) {
		DEBUG(("LBc[%ld]=%ld\n",i,LBc[i]));
	}

	// Stop indices are the first in next column.
	const mwIndex *LBs=LBc+blk+1;

	for (int i=0; i<3; i++) {
		DEBUG(("LBs[%ld]=%ld\n",i,LBs[i]));
	}

	// Indices into each column within block.
	mwIndex LBix[3];
	// Initialize indices to top of each column.
	memcpy(LBix,LBc+blk,3* sizeof(*LBc));

	// Output buffer.
	mxDouble Y[3];

	for (int row=1; row<=LBm; row++) {
		DEBUG(("\nrow=%ld\n",row));
		computeRowInBlock(LBr,LBs,LBv,LBix,row-1,UA,Y);

		for (int i=0; i<3; i++) {
			DEBUG(("Y[%ld]=%g\n",i,Y[i]));
		}
	}

#if 0
	// Buffer for inner products. Will only use 3 for
	// diagonal-only computations.
	mxDouble IP[6];

	// Compute diagonal inner products of UA block.
	computeUAinnerProductsBlock1(UA,IP);

	DEBUG(("\nIP printout\n"));

	// Print IP values
	for (int i=0; i<3; i++) {
		DEBUG(("IP[%ld]=%g\n",i,IP[i]));
	}

	// Compute all inner products of UA block.
	computeUAinnerProductsBlock3(UA,IP);

	DEBUG(("\nIP printout\n"));

	// Print IP values
	for (int i=0; i<6; i++) {
		DEBUG(("IP[%ld]=%g\n",i,IP[i]));
	}
#endif

	// Local buffer
	double *buf=mxCalloc(3*LBm,sizeof(*buf));

	if (buf==NULL) {
		mexErrMsgIdAndTxt("DBAT:icpc_mex:outOfMemory",
				"Out of memory for local buffer.");
	}

	if (buf) {
		mxFree(buf);
		buf=NULL;
	}
}
