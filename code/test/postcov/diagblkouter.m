function C=diagblkouter(A,B,BT,n)
%DIAGBLKOUTER Perform an outer product, keeping only diagonal blocks.
%
%   C=DIAGBLKOUTER(A,B,N) computes the N-by-N diagonal blocks of the
%   product B'*A*B, where the array A is real double, full and
%   symmetric. The array B is real of conforming size. The number of
%   columns of B is an integer multiple of N. The return matrix C is
%   sparse with N-by-N diagonal blocks.
%
%   C=DIAGBLKOUTER(A,B,BT,N) computes the same result for sparse B,
%   where BT is assumed to be BT=B'.
%
%   No test is performed as to whether A is symmetric not that BT=B.'.
    
error('MEX file not found.');
