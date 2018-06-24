function p=blkcolperm(M,bixIO,bixEO,bixOP)
%BLKCOLPERM Block column permutation.
%
%   p=BLKCOLPERM(M,bixIO,bixEO,bixOP) returns a permutation vector p the
%   reorders the block columns of the sparse N-by-N matrbix M in
%   nondecreasing order of nonzero count. The bixNN matrices are index
%   matrices where each column is treated as a block. Only non-zero entries
%   in the bixNN matrices are considered. The number of non-zero indices in
%   bixIO, bixEO, and bixOP should sum to N.
%
%   The matrix M(p,p) is a reordering of M such that the most dense
%   block columns are moved to the right.
%
%See also: COLPERM.


[m,n]=size(M);

if m~=n, error('DBAT:blkcolperm:badSize','M matrix must be square'); end

% Remove any all-zero columns.
bixIO(:,all(bixIO==0,1))=[];
bixEO(:,all(bixEO==0,1))=[];
bixOP(:,all(bixOP==0,1))=[];

% Stack zero-row-padded IO matrices into one index matrix.

% Determine largest block.
maxRows=max([0,find(~all(bixIO==0,2),1,'last'),...
             find(~all(bixEO==0,2),1,'last'),find(~all(bixOP==0,2),1,'last')]);

% Create main index matrix.
ix=zeros(maxRows,size(bixIO,2)+size(bixEO,2)+size(bixOP,2));

% Insert bixIO.
[mIO,nIO]=size(bixIO);
ix(1:min(maxRows,mIO),1:nIO)=bixIO(1:min(maxRows,mIO),:);

% Insert bixEO.
[mEO,nEO]=size(bixEO);
ix(1:min(maxRows,mEO),nIO+(1:nEO))=bixEO(1:min(maxRows,mEO),:);

% Insert bixOP.
[mOP,nOP]=size(bixOP);
ix(1:min(maxRows,mOP),nIO+nEO+(1:nOP))=bixOP(1:min(maxRows,mOP),:);

% Verify that we have all indices.
bix=sort(ix(ix~=0));

if length(bix)~=m || any(bix'~=1:m)
    error('DBAT:blkcolperm:invalidArgument','Bad index count');
end

% Column count first, prepend zero to enable vectorized summation.
colCount=[0,full(sum(M~=0,1))];

% Compute the average column count within each block.
blockCount=sum(colCount(ix+1),1);
avgCount=blockCount./sum(ix~=0,1);

% Sort block column count nondecreasingly.
[~,q]=sort(avgCount);

% Extract column indices in correct block order.
p=ix(:,q);
% Unroll and remove zero indices.
p=reshape(p(p~=0),1,[]);

