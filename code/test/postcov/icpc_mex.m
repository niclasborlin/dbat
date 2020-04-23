function C=icpc_mex(LA,LB,LC,N)
%ICPC_MEX Inverse Cholesky Post Cov computation
%
%   C=ICPC_MEX(LA,LB,LC,N) computes a subset of the posterior
%   covariance matrix C of the object points from the Cholesky
%   factorisation L=[LA, 0; LB, LC] of the normal matrix. The matrix
%   LA is 3M-by-3M sparse, block diagonal 3-by-3 with lower triangular
%   blocks. The matrix LB is K-by-3M and sparse. The K-by-K matrix LC
%   is full lower triangular. All matrices are real.
%
%   The returned 3M-by-3M matrix C is symmetric block diagonal with
%   N-by-N blocks, where N is 1 or 3.

error('Mex file not found.');
