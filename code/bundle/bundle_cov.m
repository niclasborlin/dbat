function varargout=bundle_cov(s,e,varargin)
%BUNDLE_COV Compute covariances from bundle result.
%
%   BUNDLE_COV(S,E,CC) computes and returns a component covariance matrix
%   specified by the string CC from the BUNDLE result files S and E. CC
%   should be one of
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
%   E=BUNDLE_COV(S,E,'prepare') performs some initial factorisations
%   that will speed up later covariance computations.
%
%See also: BUNDLE.


if isempty(varargin)
    if nargout>0
        error('DBAT:BUNDLE_COV:badInput','Too few input arguments');
    end
    return;
end

doPrepare=false;

for i=1:length(varargin)
    varargin{i}=lower(varargin{i});
    switch lower(varargin{i})
      case {'cxx','cio','ceo','cop','ciof','ceof','copf'}
        % OK, do nothing.
      case 'prepare'
        doPrepare=true;
      otherwise
        error('DBAT:bundle_cov:badInput',...
              ['Bad covariance string ''',varargin{i},'''']);
    end
end

if doPrepare && ~isempty(e.final.factorized)
    % Preparation already done, return.
    varargout{1}=e;
    return;
end

if isempty(e.final.factorized) || doPrepare
    % We may need J'*J many times. Precalculate and prefactor.
    JTJ=e.final.weighted.J'*e.final.weighted.J;
    
    % Use block column count reordering to reduce fill-in in Cholesky factor.
    
    % IO blocks.
    [i,j]=ind2sub(size(s.bundle.est.IO),s.bundle.serial.IO.src);
    bixIO=full(sparse(i,j,s.bundle.serial.IO.dest,size(s.bundle.est.IO,1),size(s.bundle.est.IO,2)));
    % EO blocks.
    [i,j]=ind2sub(size(s.bundle.est.EO),s.bundle.serial.EO.src);
    bixEO=full(sparse(i,j,s.bundle.serial.EO.dest,size(s.bundle.est.EO,1),size(s.bundle.est.EO,2)));
    % OP blocks.
    [i,j]=ind2sub(size(s.bundle.est.OP),s.bundle.serial.OP.src);
    bixOP=full(sparse(i,j,s.bundle.serial.OP.dest,size(s.bundle.est.OP,1),size(s.bundle.est.OP,2)));

    %p=blkcolperm(JTJ,bixIO,bixEO,bixOP);
    % Permute OP first, then EO, then IO.
    p=[bixOP(:);bixEO(:);bixIO(:)];
    p=p(p~=0);

    % Perform Cholesky on permuted J'*J.
    [LT,fail]=chol(JTJ(p,p));

    if fail==0
        L=LT';
        % Number of OP parameters
        nOP=nnz(bixOP);
        % Extract blocks of L = [ A, 0; B, C].
        % Diagonal OP block of L
        LA=LT(1:nOP,1:nOP)';
        % Diagonal non-OP block of L
        LC=full(LT(nOP+1:end,nOP+1:end))';
        % Subdiagonal block
        LB=full(LT(1:nOP,nOP+1:end))';
    else
        warning(['Posterior covariance matrix was not positive definite. ' ...
                 'Results will be inaccurate.'])
        n=size(JTJ,1);
        L=sparse(1:n,1:n,nan,n,n);
        LA=[];
        LB=[];
        LC=[];
    end
    Lblocks=struct('LA',LA,'LB',LB,'LC',LC);
    e.final.factorized=struct('p',p,'L',L,'fail',fail,'Lblocks',Lblocks);
    ok=~fail;
    if doPrepare
        varargout{1}=e;
        return;
    end
else
    p=e.final.factorized.p;
    L=e.final.factorized.L;
    Lblocks=e.final.factorized.Lblocks;
    fail=e.final.factorized.fail;
end

% Memory limit in elements.
memLimit=1e7;

warnState=[];
if fail
    % Turn off some numerical warnings
    warnState=warning('off','MATLAB:singularMatrix');
    warnState(2)=warning('off','MATLAB:nearlySingularMatrix');
end

for i=1:length(varargin)
    switch varargin{i}
      case 'cxx' % Raw, whole covariance matrix.

        if fail
            C=nan(size(L));
        else
            C=invblock(L,p,1:size(L,1),'direct');
        end
        
      case 'ciof' % Whole CIO covariance matrix.
            
        % Pre-allocate matrix with place for nnz(s.estIO)^2 elements.
        C=spalloc(numel(s.IO.val),numel(s.IO.val),nnz(s.bundle.est.IO)^2);
        
        if fail
            C(s.bundle.deserial.IO.dest,s.bundle.deserial.IO.dest)=nan;
        else
            % Compute needed part of inverse and put it into the right
            % part of C.
            C(s.bundle.deserial.IO.dest,s.bundle.deserial.IO.dest)=...
                invblock(L,p,s.bundle.deserial.IO.src,'sqrt');
        end
        
      case 'ceof' % Whole CEO covariance matrix.
        
        start=clock; %#ok<NASGU>
        % Pre-allocate matrix with place for nnz(s.estEO)^2 elements.
        C=spalloc(numel(s.EO.val),numel(s.EO.val),nnz(s.bundle.est.EO)^2);
        
        if fail
            C(s.bundle.deserial.EO.dest,s.bundle.deserial.EO.dest)=nan;
        else
            % Compute needed part of inverse and put it into the right
            % part of C.
            C(s.bundle.deserial.EO.dest,s.bundle.deserial.EO.dest)=...
                invblock(L,p,s.bundle.deserial.EO.src,'sqrt');
        end
        
        %etime(clock,start)
        
      case 'copf' % Whole COP covariance matrix.
        
        start=clock; %#ok<NASGU>
        % Pre-allocate matrix with place for nnz(s.estOP)^2 elements.
        C=spalloc(numel(s.OP.val),numel(s.OP.val),nnz(s.bundle.est.OP)^2);
        
        if fail
            C(s.bundle.deserial.OP.dest,s.bundle.deserial.OP.dest)=nan;
        else
            % Compute needed part of inverse and put it into the right
            % part of C.
            C(s.bundle.deserial.OP.dest,s.bundle.deserial.OP.dest)=...
                invblock(L,p,s.bundle.deserial.OP.src,'sqrt');
        end
        
        %etime(clock,start)
        
      case 'cio' % Block-diagonal CIO
        
        C=BlockDiagonalC(L,p,s.bundle.est.IO,s.bundle.deserial.IO.src,...
                         memLimit,'Computing IO covariances');
        
      case 'ceo' % Block-diagonal CEO
        
        if ~isempty(s.post.cov.CEO)
            C=s.post.cov.CEO;
        else
            C=BlockDiagonalC(L,p,s.bundle.est.EO,s.bundle.deserial.EO.src,...
                             memLimit,'Computing EO covariances');
        end
      case 'cop' % Block-diagonal COP
        if ~isempty(s.post.cov.COP)
            C=s.post.cov.COP;
        else
            tic
            C=BlockDiagonalC(L,p,s.bundle.est.OP,s.bundle.deserial.OP.src,...
                             memLimit,'Computing OP covariances');
            toc
            tic
            CC=VectorizedCOP(Lblocks,s.bundle.est.OP);
            toc
            
            cm2=full(diag(C,-2));
            ccm2=full(diag(CC,-2));
            cm1=full(diag(C,-1));
            ccm1=full(diag(CC,-1));
            c0=full(diag(C,0));
            cc0=full(diag(CC,0));
            c1=full(diag(C,1));
            cc1=full(diag(CC,1));
            c2=full(diag(C,2));
            cc2=full(diag(CC,2));

            dm2=abs(cm2-ccm2)./abs(cm2);
            dm1=abs(cm1-ccm1)./abs(cm1);
            d0=abs(c0-cc0)./abs(c0);
            d1=abs(c1-cc1)./abs(c1);
            d2=abs(c2-cc2)./abs(c2);
            max([max(dm2),max(dm1),max(d0),max(d1),max(d2)])
        end
    end
    varargout{i}=e.s0^2*C; %#ok<*AGROW>
end

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

function CC=VectorizedCOP(Lblocks,calc)

% L = [LA, 0; LB, LC]. LA is square diagonal block corresponding to
% the OPs. LC is square diagonal block corresponding to the non-OPs.
% LB is rectangular subdiagonal block. LA is sparse block-diagonal
% with lower triangular 3-by-3 blocks. LB is dense. LC is dense lower
% triangular.
LA=Lblocks.LA;
LB=Lblocks.LB;
LC=Lblocks.LC;

% So, L*L'=J'*J.
% We want inv(J'*J) = inv(L*L') = inv(L')*inv(L) = U'*U.

% With U = [UA, 0; UB, UC] blocked as L, we are only interested in
% UA and UB.

LCs=sparse(LC);

% Invert LA. Actually faster than LA\speye(size(LA)).
UA=inv(LA);
UB=-(LC\(LB*UA));

% Extract staggered columns of U with stride 3
ua1=UA(:,1:3:end);
ua2=UA(:,2:3:end);
ua3=UA(:,3:3:end);
ub1=UB(:,1:3:end);
ub2=UB(:,2:3:end);
ub3=UB(:,3:3:end);

% Each diagonal 3-by-3 block of
% U'*U = [u11 u12 u13
%         u12 u22 u23
%         u13 u23 u33]

% Compute diagonal elements of U'*U. ud contains 3 elements per
% block [u11 u22 u33]
ud0=(sum(UA.^2,1)+sum(UB.^2,1))';

% Compute the (1,2), (1,3), (2,3) elements of each block. Each
% vector will contain one element per block.
u12=sum(ua1.*ua2,1)+sum(ub1.*ub2,1);
u13=sum(ua1.*ua3,1)+sum(ub1.*ub3,1);
u23=sum(ua2.*ua3,1)+sum(ub2.*ub3,1);

% Create superdiagonal vectors
udm1=reshape([u12;u23;zeros(size(u12))],[],1);
udm2=reshape([u13;zeros(2,length(u13))],[],1);
ud1=reshape([zeros(size(u12));u12;u23],[],1);
ud2=reshape([zeros(2,length(u13));u13],[],1);

% Preallocate place for all diagonals, including those
% corresponding to unestimated OPs.

ud=zeros(numel(calc),5);
if ~isempty(udm2)
    ud(calc(:),1)=udm2;
    ud(calc(:),2)=udm1;
    ud(calc(:),3)=ud0;
    ud(calc(:),4)=ud1;
    ud(calc(:),5)=ud2;
end

CC=spdiags(ud,-2:2,size(ud,1),size(ud,1));
