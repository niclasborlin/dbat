function varargout=bundle_cov(s,e,varargin)
%BUNDLE_COV Compute covariances from bundle result.
%
%   BUNDLE_COV(S,E,CC) computes and returns a component covariance matrix
%   specified by the string CC. CC should be one of 
%       'CIO' - covariance of estimated internal parameters,
%       'CEO' - covariance of estimated external parameters, or
%       'COP' - covariance of estimated object points.
%   The returned matrices will be sparse,
%   block-diagonal matrices - 6x6 block for CEO, 3x3 for COP, and
%   camera-model-dependent for CIO. The blocks are zero-padded for fixed
%   elements, i.e. elements not estimated by BUNDLE. The blocking
%   corresponds to that any correlations with other cameras/images/points
%   are ignored.
%
%   To include camera-camera, image-image, etc., covariance matrices, use
%   CC='CIOF', 'CEOF', 'COPF', where e.g. CEOF will be a full
%   (6*N)-by-(6*N) matrix, where N is the number of images.
%
%   The special covariance matrix CC='CXX' returns the CXX matrix that
%   corresponds directly to the estimated vector. CXX is not zero-padded.
%
%   Warning! Especially CXX and COPF may require a lot of memory!
%
%   [CA,CB,...]=BUNDLE_COV(S,E,CCA,CCB,...) returns multiple covariance
%   matrices CA, CB, etc.
%
%See also: BUNDLE.

% $Id$

if isempty(varargin)
    return;
end

for i=1:length(varargin)
    varargin{i}=lower(varargin{i});
    switch lower(varargin{i})
      case {'cxx','cio','ceo','cop','ciof','ceof','copf'}
        % OK, do nothing.
      otherwise
        error('DBAT:bundle_cov:badInput',...
              ['Bad covariance string ''',varargin{i},'''']);
    end
end

% Create indices into the vector of unknowns.
[ixIO,ixEO,ixOP]=indvec([nnz(s.cIO),nnz(s.cEO),nnz(s.cOP)]);

% We may need J'*J many times. Precalculate and prefactor.
JTJ=e.s0^2*e.J'*e.J;
    
% Use block column count reordering to reduce fill-in in Cholesky factor.
    
% IO blocks.
bixIO=double(s.cIO);
bixIO(s.cIO)=ixIO;
% EO blocks.
bixEO=double(s.cEO);
bixEO(s.cEO)=ixEO;
% OP blocks.
bixOP=double(s.cOP);
bixOP(s.cOP)=ixOP;
    
p=blkcolperm(JTJ,bixIO,bixEO,bixOP);

% Perform Cholesky on permuted J'*J.
R=chol(JTJ(p,p));

% Construct inverse permutation.
invP=zeros(size(p));
invP(p)=1:length(p);
    
for i=1:length(varargin)
    switch varargin{i}
      case 'cxx' % Raw, whole covariance matrix.

        C=invblock(R,p,1:size(R,1),'direct');
            
      case 'ciof' % Whole CIO covariance matrix.
            
        % Pre-allocate matrix with place for nnz(s.cIO)^2 elements.
        C=spalloc(numel(s.IO),numel(s.IO),nnz(s.cIO)^2);
        
        % Compute needed part of inverse and put it into the right
        % part of C.
        C(s.cIO(:),s.cIO(:))=invblock(R,p,ixIO,'split');
        
      case 'ceof' % Whole CEO covariance matrix.
        
        start=clock;
        % Pre-allocate matrix with place for nnz(s.cEO)^2 elements.
        C=spalloc(numel(s.EO),numel(s.EO),nnz(s.cEO)^2);
        
        % Compute needed part of inverse and put it into the right
        % part of C.
        C(s.cEO(:),s.cEO(:))=invblock(R,p,ixEO,'split');

        % Remove axis indicator psuedo-elements.
        keep=rem(1:size(C,1),7)~=0;
        C=C(keep,keep);
        
        %etime(clock,start)
        
      case 'copf' % Whole COP covariance matrix.
        
        start=clock;
        % Pre-allocate matrix with place for nnz(s.cOP)^2 elements.
        C=spalloc(numel(s.OP),numel(s.OP),nnz(s.cOP)^2);
        
        % Compute needed part of inverse and put it into the right
        % part of C.
        C(s.cOP(:),s.cOP(:))=invblock(R,p,ixOP,'split');
        
        %etime(clock,start)
        
      case 'cio' % Block-diagonal CIO
        
        C=BlockDiagonalC(R,p,s.cIO,ixIO,1,'Computing IO covariances');
        
      case 'ceo' % Block-diagonal CEO
        
        C=BlockDiagonalC(R,p,s.cEO,ixEO,1,'Computing EO covariances');

        % Remove axis indicator psuedo-elements.
        keep=rem(1:size(C,1),7)~=0;
        C=C(keep,keep);
        
      case 'cop' % Block-diagonal COP
        
        C=BlockDiagonalC(R,p,s.cOP,ixOP,1,'Computing OP covariances');
        
    end
    varargout{i}=C;
end


function C=BlockDiagonalC(R,p,calc,xIx,bsElems,msg)
%R       - Cholesky factor of the permuted normal matrix A(p,p).
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
bsCols=floor(bsElems/size(R,1)/max(sum(ix~=0,1)));
bsCols=min(max(bsCols,1),size(ix,2));

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
[dummy,colPerm]=sort(oCol);

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
    Cblock(eix,eix)=invblock(R,p,jix,'split');
    % Make it block-diagonal...
    Cblock=mkblkdiag(Cblock,size(ix,1));
    % ...and put it into the right part of C.
    C(eix,eix)=Cblock(eix,eix);
    
    if isempty(h) && etime(clock,start)>1 && j+bsCols-1<size(ix,2)
        % Only create dialog if execution takes more than 1s and this
        % iteration is not the last.
        h=waitbar(min(j+bsCols-1,size(ix,2))/size(ix,2),msg);
        lapTime=clock;
    elseif etime(clock,lapTime)>1
        % Update dialog every 1 s.
        if ishandle(h) % Guard against window close.
            waitbar(min(j+bsCols-1,size(ix,2))/size(ix,2),h);
        end
        lapTime=clock;
    end
end
if ishandle(h), close(h), end
%etime(clock,start)
