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
    % Time the preparation.
    stopWatch=cputime;
    
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
    prepTime=cputime-stopWatch;
    e.final.factorized=struct('p',p,'L',L,'Lblocks',Lblocks,'fail',fail,...
                              'prepTime',prepTime);
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
            C=VectorizedCOP(Lblocks,s.bundle.est.OP);
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

function C=VectorizedCOP(Lblocks,calc)

% L = [LA, 0; LB, LC]. LA is square diagonal block corresponding to
% the OPs. LC is square diagonal block corresponding to the non-OPs.
% LB is rectangular subdiagonal block. LA is sparse block-diagonal
% with lower triangular 3-by-3 blocks. LB is dense. LC is dense lower
% triangular.
LA=Lblocks.LA;
LB=Lblocks.LB;
LC=Lblocks.LC;

% Early return if prior error signalled by empty LA.
if isempty(LA)
    nans=zeros(numel(calc),1);
    nans(calc)=nan;
    C=spdiags(nans,0,numel(calc),numel(calc));
    return;
end

% So, L*L'=J'*J.
% We want inv(J'*J) = inv(L*L') = inv(L')*inv(L) = K'*K.

% With K = [KA, 0; KB, KC] blocked as L, we are only interested in
% P=KA and Q=KB.

% Invert LA. Actually faster than LA\speye(size(LA)).
P=inv(LA);

% Compute LB*inv(LA).
LBP=LB*P; %#ok<MINV>

% Indices of elements that have been estimated. The index calcIx(i)
% indicates what element in COP the column P(:,i) corresponds to.
calcIx=find(calc(:));
% Indices for elements in each row
ixIs1=rem(calcIx-1,3)==0;
calcIx1=find(calc(1,:));
ixIs2=rem(calcIx-1,3)==1;
calcIx2=find(calc(2,:));
ixIs3=rem(calcIx-1,3)==2;
calcIx3=find(calc(3,:));

% Compute diagonal elements of P'*P. c0 contains 3 elements per block
% [c11 c22 c33]. Introduce zeros for covariances that are not to be
% computed.
if all(calc(:))
    c0=full(sum(P.^2,1));
else
    c0=zeros(1,numel(calc));
    c0(calcIx)=full(sum(P.^2,1));
end
    
if false && onlyDiag
    c12=[];
    c13=[];
    c23=[];
else
    % Compute off-diagonal elements
    if all(calc(:))
        % Extract staggered columns of P with stride 3
        p1=P(:,1:3:end);
        p2=P(:,2:3:end);
        p3=P(:,3:3:end);
    else
        % Extract staggered columns of P with stride 3. Expand with zeros for
        % elements that have not been estimated.
        p1=sparse(size(P,1),size(calc,2));
        p2=p1;
        p3=p1;
        p1(:,calcIx1)=P(:,ixIs1);
        p2(:,calcIx2)=P(:,ixIs2);
        p3(:,calcIx3)=P(:,ixIs3);
    end        

    % Compute the (1,2), (1,3), (2,3) elements of each block. Each
    % vector will contain one element per block.
    c12=full(sum(p1.*p2,1));
    c13=full(sum(p1.*p3,1));
    c23=full(sum(p2.*p3,1));
end

% Computing Q directly may cause out of memory. Accept blocks of about
% 256MB.

blockSize=256*1024^2;
blockCols=floor(min(round(blockSize/8/size(LB,1)),size(LB,2))/3)*3;

% Index of first column in calcIx to use.
baseIx=1;

while baseIx<=length(calcIx)
    % Corresponding base column in calc
    baseCalcCol=floor((calcIx(baseIx)-1)/3)+1;
    % Look blockCols forward.
    lastIx=min(baseIx+blockCols-1,length(calcIx));
    % Corresponding column in calc
    lastCalcCol=floor((calcIx(lastIx)-1)/3)+1;
    % Adjust to include all elements of last column.
    lastIx=find(calcIx<=lastCalcCol*3,1,'last');
    
    % Columns in LBP in this block.
    ix=baseIx:lastIx;
    % Corresponding columns in c0
    cix=calcIx(ix);
    % Corresponding columns in cols
    calcCols=baseCalcCol:lastCalcCol;

    LBPblk=LBP(:,ix);
    Q=-LC\LBPblk;

    % Update c0 elements with diagonal elements of Q'*Q.
    c0(cix)=c0(cix)+sum(Q.^2,1);

    if true || ~onlyDiag
        % Compute off-diagonal elements
        if all(calc(:))
            % Extract staggered columns of Q with stride 3
            q1=Q(:,1:3:end);
            q2=Q(:,2:3:end);
            q3=Q(:,3:3:end);
        else
            % Extract staggered columns of Q with stride 3. Expand with zeros for
            % elements that have not been estimated.
            q1=zeros(size(Q,1),length(calcCols));
            q2=q1;
            q3=q1;
            q1(:,calc(1,calcCols))=Q(:,rem(cix-1,3)==0);
            q2(:,calc(2,calcCols))=Q(:,rem(cix-1,3)==1);
            q3(:,calc(3,calcCols))=Q(:,rem(cix-1,3)==2);
        end        

        % Compute the (1,2), (1,3), (2,3) elements of each block. Each
        % vector will contain one element per block.
        c12(calcCols)=c12(calcCols)+sum(q1.*q2,1);
        c13(calcCols)=c13(calcCols)+sum(q1.*q3,1);
        c23(calcCols)=c23(calcCols)+sum(q2.*q3,1);
    end

    baseIx=lastIx+1;
end

if false && onlyDiag
    C=spdiags(ud0(:),0,length(ud0),length(ud0));
else
    % Create subdiagonal vectors...
    cdm1=reshape([c12;c23;zeros(size(c12))],[],1);
    cdm2=reshape([c13;zeros(2,length(c13))],[],1);
    % ...and superdiagonal vectors
    cd1=reshape([zeros(size(c12));c12;c23],[],1);
    cd2=reshape([zeros(2,length(c13));c13],[],1);

    % Preallocate place for all diagonals, including those
    % corresponding to unestimated OPs.

    cd=zeros(numel(calc),5);
    cd(:,1)=cdm2;
    cd(:,2)=cdm1;
    cd(:,3)=c0;
    cd(:,4)=cd1;
    cd(:,5)=cd2;

    C=spdiags(cd,-2:2,size(cd,1),size(cd,1));
end
