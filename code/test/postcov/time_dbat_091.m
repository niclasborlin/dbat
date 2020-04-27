function [times,C]=time_dbat_091(s,e)
%Returned times are [chol,...,total];

% Code from bundle_cov.m

% We may need J'*J many times. Precalculate and prefactor.
JTJ=e.final.weighted.J'*e.final.weighted.J;

% IO blocks.
[i,j]=ind2sub(size(s.bundle.est.IO),s.bundle.serial.IO.src);
bixIO=full(sparse(i,j,s.bundle.serial.IO.dest,size(s.bundle.est.IO,1),size(s.bundle.est.IO,2)));
% EO blocks.
[i,j]=ind2sub(size(s.bundle.est.EO),s.bundle.serial.EO.src);
bixEO=full(sparse(i,j,s.bundle.serial.EO.dest,size(s.bundle.est.EO,1),size(s.bundle.est.EO,2)));
% OP blocks.
[i,j]=ind2sub(size(s.bundle.est.OP),s.bundle.serial.OP.src);
bixOP=full(sparse(i,j,s.bundle.serial.OP.dest,size(s.bundle.est.OP,1),size(s.bundle.est.OP,2)));

p=blkcolperm(JTJ,bixIO,bixEO,bixOP);
JTJ=JTJ(p,p);

startClock=now;

% Perform Cholesky on permuted J'*J.
[LT,fail]=chol(JTJ);

L=LT';

% Memory limit in elements.
memLimit=1e7;

warnState=[];
% Turn off some numerical warnings
warnState=warning('off','MATLAB:singularMatrix');
warnState(2)=warning('off','MATLAB:nearlySingularMatrix');

C=BlockDiagonalC(L,p,s.bundle.est.OP,s.bundle.deserial.OP.src,...
                 memLimit,'Computing OP covariances');

C=C(s.bundle.est.OP(:),s.bundle.est.OP(:));

stopClock=now;
times=(stopClock-startClock)*86400;

if ~isempty(warnState)
    % Restore warning state.
    warning(warnState);
end

function C=BlockDiagonalC(L,p,calc,xIx,bsElems,msg)
%L       - Cholesky factor of the permuted normal matrix A(p,p).
%p       - Corresponding permutation vector.
%calc    - logical M-by-N array indicating which data elements have been
%          estimated and whose covariances should be computed.
%xIx     - vector of length nnz(calc) with indices of the estimated
%          elements in the x vector.
%bsElems - how many elements of the inverse should at most be calculated
%          at each block iteration? bsElems=1 will use unblocked algorithm.
%msg     - waitbar message to present. No waitbar is shown if message is empty.
%C       - M-by-M block-diagonal covariance matrix of size
%          nnz(calc)-by-nnz(calc).

% Delayed progress dialog.
start=clock;
lapTime=start;
h=[];

% Pre-allocate sparse matrix with place for all diagonal blocks.
C=spalloc(numel(calc),numel(calc),sum(sum(calc,1).^2));

% Pack indices as the data
ix=zeros(size(calc));
ix(calc)=xIx;

% Construct corresponding indexing among the OP parameters.
ixInt=reshape(1:numel(calc),size(calc));

% Determine block column size such that computed part of inverse is
% approximately bsElems elements.
bsCols=floor(bsElems/size(L,1)/max(sum(ix~=0,1)));
bsCols=min(max(bsCols,1),size(ix,2));

% This function is not memory-limited. A block column size of 30
% seems to give reasonably good results.
% TODO: Refine bsCols value.
bsCols=30;

% Create inverse permutation.
invP=zeros(size(p));
invP(p)=1:length(p);

% Sort by original column.

% Current column (one from each block). Will be zero for fixed data.
cCol=max(ix,[],1);
% Corresponding original column.
oCol=zeros(size(cCol));
oCol(cCol>0)=invP(cCol(cCol>0));

% Get permutation to sort by increasing original column number. This will
% move fixed columns to the beginning, becoming part of the first blocks.
[~,colPerm]=sort(oCol);

% Loop over each column in bsCols blocks.
for j=1:bsCols:size(ix,2)
    % Columns in block.
    jCols=colPerm(j:min(j+bsCols-1,size(ix,2)));

    % Indices into J.
    jixBlock=ix(:,jCols);
    jix=jixBlock(calc(:,jCols));
    % Indices into data.
    eixBlock=ixInt(:,jCols);
    eix=eixBlock(calc(:,jCols));

    % Compute needed part of inverse and put it into a temporary matrix.
    Cblock=spalloc(size(C,1),size(C,1),length(eix)^2);
    Z2=invblock(L,p,jix,'sqrt');
    %Z3=invblock(L,p,jix,'sqrtsplit');
    %Z4=invblock(L,p,jix,'direct');
    %max(max(abs(Z2-Z3)))
    %max(max(abs(Z3-Z4)))
    Cblock(eix,eix)=Z2;
    % Make it block-diagonal...
    Cblock=mkblkdiag(Cblock,size(ix,1));
    % ...and put it into the right part of C.
    C(eix,eix)=Cblock(eix,eix);
end
%etime(clock,start)

function x=invblock(L,p,ix,method)
%INVBLOCK Compute one block of normal matrix inverse.
%
%   X=INVBLOCK(L,p,IX) computes the block M(IX,IX) of the inverse M of the
%   normal matrix N. L should have been computed as chol(N(p,p))', where p
%   is a permutation vector.


if isempty(ix)
    x=[];
    return;
end

if isnan(L(1,1))
    x=nan(nnz(ix));
end

% Create inverse permutation.
invP=zeros(size(p));
invP(p)=1:length(p);

switch method
  case 'verify'
    % Reconstruct N and solve via full inverse.
    LLT=L*L';
    N=LLT(invP,invP);
    invN=inv(N);
    x=invN(ix,ix);
  case 'sqrt'
    % Let E contain the wanted columns of I. Then our wanted block is
    % B=E'*(A\E)=E'*(inv(A)*E). With A=L*L',
    % B=E'*(inv(L*L')*E)=E'*inv(L')*inv(L)*E. With W=L\E=inv(L)*E, B=W'*W.

    % Permute the identity matrix.
    Ip=speye(size(L,2));
    Ip=Ip(p,:);

    % Run computation on multiple rows to enable profiling.
    v1=Ip(:,ix);
    v2=L\v1;
    x=v2'*v2;
  case 'sqrtsplit'
    % Let E contain the wanted columns of I. Then our wanted block is
    % B=E'*(A\E)=E'*(inv(A)*E). With A=L*L',
    % B=E'*(inv(L*L')*E)=E'*inv(L')*inv(L)*E. With W=L\E=inv(L)*E, B=W'*W.

    % Do a split of L into [L1,  0;
    %                       L21, L2]
    % L1 will generally be block-diagonal - use sparse storage.
    % L21 will generally be sparse - use sparse storage.
    % L2 will generally be dense - use full storage.
    
    % Permute the identity matrix.
    Ip=speye(size(L,2));
    Ip=Ip(p,:);

    % Find jump in row density.
    colDens=full(sum(L~=0,2))/size(L,1);
    [~,i]=max(diff(colDens));
    ix1=1:i;
    ix2=i+1:size(L,1);
    
    % Extract (1,1) block and (2,1) block as sparse matrices.
    L1=L(ix1,ix1);
    L21=L(ix2,ix1);
    % Extract (2,2) block as dense matrix.
    L2=full(L(ix2,ix2));

    % Permute the identity matrix.
    Ip=speye(size(L,2));
    Ip=Ip(p,:);

    % Run computation on multiple rows to enable profiling.
    % [L1    0  ] [ a ]   [ c ]
    % [L21   L2 ] [ b ] = [ d ]
    % a=L1\c
    % b=L2\(d-L21*a)
    cd=Ip(:,ix);
    c=cd(ix1,:);
    d=cd(ix2,:);
    a=L1\c;
    L21a=L21*a;
    dmL21a=d-L21a;
    b=L2\dmL21a;
    x=a'*a+b'*b;
  case 'direct'
    Ip=speye(size(L,2));
    Ip=Ip(p,:);
    A=L'\(L\Ip(:,ix));
    x=A(invP(ix),:);
  otherwise
    error('Bad method');
end

