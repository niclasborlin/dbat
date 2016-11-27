function A=mkblkdiag(C,varargin)
%MKBLKDIAG Make a matrix block diagonal.
%
%   A=MKBLKDIAG(C,N) return a block diagonal copy of C with NxN-sized blocks.
%
%   A=MKBLKDIAG(C,M,N) and A=MKBLKDIAG(C,[M,N]) does the same for
%   MxN-sized diagonal block.
%
%   Default is to return A with same storage mode as C.  Use
%   A=MKBLKDIAG(...,true) or A=MKBLKDIAG(...,'sparse') to enforce sparse
%   return and A=MKBLKDIAG(...,false) or A=MKBLKDIAG(...,'full') to
%   enforce full return.
%
%See also: blkdiag, mkblkdiag.


narginchk(2,inf);

useSparse=issparse(C);
if (~isempty(varargin))
    if (islogical(varargin{end}))
        useSparse=varargin{end};
        varargin(end)=[];
    elseif (ischar(varargin{end}))
        arg=[varargin{end},' '];
        switch (lower(arg(1)))
        case 's'
            useSparse=true;
        case 'f'
            useSparse=false;
        otherwise
            error('Bad string argument');
        end
        varargin(end)=[];
    end
end

if (length(varargin)==1)
    m=varargin{1};
    if (length(m)==1)
        n=m;
    else
        [m,n]=deal(m(1),m(2));
    end
else
    [m,n]=deal(varargin{1:2});
end

% Find non-zero elements of C.
[i,j,v]=find(C);
% Calculate block indices.
bi=floor((i-1)/m)+1;
bj=floor((j-1)/n)+1;
onBlkDiag=bi==bj;

% Create a new matrix containing the block diagonal elements.
A=sparse(i(onBlkDiag),j(onBlkDiag),v(onBlkDiag),size(C,1),size(C,2));
if (~useSparse)
    A=full(A);
end
