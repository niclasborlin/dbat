function b=isblkdiag(C,m,n)
%ISBLKDIAG True for a block diagonal matrix.
%
%   ISBLKDIAG(C,M,N) returns true if C is block diagonal with MxN blocks.
%
%   ISBLKDIAG(C,[M,N]) does the same.
%
%   ISBLKDIAG(C,N) assumes NxN blocks.
%
%See also: blkdiag, mkblkdiag.


narginchk(2,3);

if (nargin==2)
    if (length(m)==1)
        n=m;
    else
        [m,n]=deal(m(1),m(2));
    end
end

% Do some cheap test(s) first.
if (nnz(C)>min(ceil(size(C)./[m,n]))*m*n)
    % More non-zeros than can fit in full diagonal blocks.
    b=false;
    return;
end

% Find non-zero elements.
[i,j]=find(C);
% Calculate block indices.
bi=floor((i-1)/m)+1;
bj=floor((j-1)/n)+1;
% All elements on block diagonal?
b=all(bi==bj);
