function [s,ok,iters,s0,E,varargout]=bundle(s,varargin)
%BUNDLE Run bundle adjustment iterations on a camera network.
%
%   [S,OK,N]=BUNDLE(S), where S is a struct returned by PROB2DBATSTRUCT,
%   runs the damped bundle adjustment on the camera network in the structure
%   S. The parameter values in S are used as initial values. The cIO, cEO,
%   cOP fields of S are used to indicate which parameters are free. OK is
%   returned as true if the bundle converged within the allowed number of
%   iterations. N gives the number of iterations. On return, the parameter
%   values in S are updated if the bundle converged.
%
%   ...=BUNDLE(S,...,K), where K is an integer, sets the maximum number of
%   iterations to K (default: 20).
%
%   ...=BUNDLE(S,...,DAMP), where DAMP is a string, specifies which damping
%   to use: 'none' or 'GM' (classic bundle with no damping), 'GNA'
%   (Gauss-Newton with Armijo linesearch, default), 'LM' (original
%   Levenberg-Marquardt) , 'LMP' (Levenberg-Marquardt with Powell dogleg).
%
%   ...=BUNDLE(S,...,TRACE), where TRACE is a string='trace' specifies
%   that the bundle should print an trace during the iterations.
%
%   ...=BUNDLE(S,...,SINGULAR_TEST), where SINGULAR_TEST is the
%   string='singular_test' specifies that the bundle should stop
%   immediately if a 'Matrix is singular' or 'Matrix is almost singular'
%   warning is issued on the normal matrix.
%
%   ...=BUNDLE(S,...,CHI), where CHI is a logical scalar, specifies if
%   chirality veto damping should be used (default: false). Chirality
%   veto damping is ignored for the undamped bundle.
%
%   [S,OK,N,S0]=... returns the sigma0 for the last iteration.
%
%   [S,OK,N,S0,E]=... returns damping-dependent trace data in the struct E:
%       E.damping - string with the name of the damping scheme.
%       E.trace   - NOBS-by-(N+1) array with successive parameter estimates.
%       E.res     - (N+1)-vector with the residual norm at every iteration.
%   Furthermore,
%       for GNA damping: E.alpha  - (N+1)-vector with used steplength.
%       for LM damping:  E.lambda - (N+1)-vector with used lambda values.
%       for LMP damping: E.delta  - (N+1)-vector with used trust region sizes,
%                        E.rho    - (N+1)-vector with gain ratios.
%
%   [S,OK,N,S0,E,CXX]=BUNDLE(S,...,'CXX') computes and returns the
%   covariance matrix CXX for all estimated parameters, scaled by sigma0^2.
%
%   [S,OK,N,S0,E,CA]=BUNDLE(S,...,CCA) computes and returns a selected
%   covariance matrix CA. CCA should be one of 'CIO' (covariance of internal
%   parameters), 'CEO' (external parameters), 'COP' (object points), or
%   'CXX' (all parameters). CIO, CEO, COP will be zero-padded for fixed
%   elements. CXX corresponds directly to the estimated vector and is not
%   zero-padded. As default, the returned component covariance matrices CIO,
%   CEO, and COP are sparse matrices with the block-diagonal part of CXX
%   only, i.e. the inter- and intra-IO/EO/OP parts are not returned. Use
%   CCA='CIOF','CEOF', or 'COPF' obtain full component covariance matrices.
%
%   Warning! Especially CXX and COPF may require a lot of memory!
%
%   [S,OK,N,S0,E,CA,CB,...]=BUNDLE(S,...,CCA,CCB,...) returns multiple
%   covariance matrices CA, CB, ...
%
%   References: Börlin, Grussenmeyer (2013), "Bundle Adjustment With and
%       Without Damping". Photogrammetric Record 28(144), pp. 396-415. DOI
%       10.1111/phor.12037.
%
%   See also PROB2DBATSTRUCT, BROWN_EULER_CAM.

% $Id$

maxIter=20;
damping='gna';
veto=false;
singular_test=false;
trace=false;
covMatrices={};

while ~isempty(varargin)
    if isnumeric(varargin{1}) && isscalar(varargin{1})
        % N
        maxIter=varargin{1};
        varargin(1)=[];
    elseif ischar(varargin{1})
        % DAMP, TRACE, or CXX
        switch lower(varargin{1})
          case {'none','gm','gna','lm','lmp'}
            % OK
            damping=varargin{1};
            varargin(1)=[];
          case 'trace'
            trace=true;
            varargin(1)=[];
          case 'singular_test'
            singular_test=true;
            varargin(1)=[];
          case {'cxx','cio','ceo','cop','ciof','ceof','copf'}
            covMatrices{end+1}=lower(varargin{1});
            varargin(1)=[];
          otherwise
            error('DBAT:bundle:badInput','Unknown damping');
        end
    elseif islogical(varargin{1})
        veto=varargin{1};
        varargin(1)=[];
    else
        error('DBAT:bundle:badInput','Unknown parameter');
    end
end

% Create indices into the vector of unknowns.
[ixIO,ixEO,ixOP]=indvec([nnz(s.cIO),nnz(s.cEO),nnz(s.cOP)]);

% Number of unknowns.
n=max(ixOP);

% Set up vector of initial values.
x0=nan(n,1);
x0(ixIO)=s.IO(s.cIO);
x0(ixEO)=s.EO(s.cEO);
x0(ixOP)=s.OP(s.cOP);

% Residual function.
resFun=@brown_euler_cam;

if veto
    vetoFun=@chirality;
else
    vetoFun='';
end

% Convergence tolerance.
convTol=1e-3;

% Set up cell array of extra parameters.
params={s};

% For all optimization methods below, the final estimate is returned in x.
% The final residual vector and jacobian are returned as f and J. Successive
% estimates of x are returned as columns of X. Furthermore, a status code (0
% - ok, -1 - too many iterations) and the number of required iterations are
% returned.

E=struct('damping','');

switch lower(damping)
  case {'none','gm'}
    % Gauss-Markov with no damping.
    
    % Call Gauss-Markov optimization routine.
    stopWatch=cputime;
    [x,code,iters,f,J,X,res]=gauss_markov(resFun,x0,maxIter,convTol,trace, ...
                                      singular_test,params);
    time=cputime-stopWatch;
    E.damping='gm';
  case 'gna'
    % Gauss-Newton with Armijo linesearch.

    % Armijo parameter.
    mu=0.1;
    % Shortest allowed step length.
    alphaMin=1e-3;
    
    % Call Gauss-Newton-Armijo optimization routine. The vector alpha is
    % returned with the step lengths used at each iteration.
    stopWatch=cputime;
    [x,code,iters,f,J,X,res,alpha]=gauss_newton_armijo(resFun,vetoFun,x0, ...
                                                   mu,alphaMin,maxIter, ...
                                                   convTol,trace, ...
                                                   singular_test,params);
    time=cputime-stopWatch;
    E.alpha=alpha;
    E.damping='gna';
  case 'lm'
    % Original Levenberg-Marquardt "lambda"-version.

    % Estimate initial lambda value.
    [f0,J0]=resfun(x0,params{:});
    lambda0=1e-10*trace(J0'*J0)/n;

    % Set value below which all lambda values are considered zero.
    lambdaMin=lambda0;

    % Call Levenberg-Marquardt optimization routine. Supply f0, J0 to avoid
    % recomputing them. The vector lambdas is returned with the lambda
    % values used at each iteration.
    stopWatch=cputime;
    [x,code,iters,f,J,X,lambda]=levenberg_marquardt(resFun,vetoFun,x0, ...
                                                    f0,J0,maxIter, ...
                                                    convTol,lambda0, ...
                                                    lambdaMin,trace,params);
    time=cputime-stopWatch;
    E.lambda=lambda;
    E.damping='lm';
  case 'lmp'
    % Levenberg-Marquardt-Powell trust-region, "delta"-version with
    % Powell dogleg.
    
    % Limits on the gain ratio.
    rhoBad=0.25; % Below this is bad.
    rhoGood=0.75; % Above this is good.
    
    % Initial delta value.
    delta0=norm(x0);
    
    % Call Levenberg-Marquardt-Powell optimization routine. The vector deltas
    % and rhos are returned with the used delta and computed rho values for
    % each iteration.
    stopWatch=cputime;
    [x,code,iters,f,J,X,delta,rho]=levenberg_marquardt_powell(resFun, ...
                                                      vetoFun,x0,delta0, ...
                                                      maxIter,convTol, ...
                                                      rhoBad,rhoGood, ...
                                                      trace,params);
    time=cputime-stopWatch;
    E.delta=delta;
    E.rho=rho;
    E.damping='lmp';
  otherwise
    error('DBAT:bundle:internal','Unknown damping');
end
E.res=res;
E.trace=X;
E.time=time;

% Handle returned values.
ok=code==0;

% Update s if optimization converged.
if ok
    s.IO(s.cIO)=x(ixIO);
    s.EO(s.cEO)=x(ixEO);
    s.OP(s.cOP)=x(ixOP);
end

% Sigma0 is sqrt(f'*f/(m-n)) in mm, convert to pixels.
s0=sqrt(f'*f/(length(f)-length(x)))*mean(s.IO(end-1:end));

% Compute covariance matrices.
if ~isempty(covMatrices)
    % We may need J'*J many times. Precalculate and prefactor.
    JTJ=s0^2*J'*J;
    
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

    % Inverse permutation.
    invP=zeros(size(p));
    invP(p)=1:length(p);
    
    for i=1:length(covMatrices)
        switch covMatrices{i}
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

          case 'cop' % Block-diagonal COP
            
            C=BlockDiagonalC(R,p,s.cOP,ixOP,1,'Computing OP covariances');

        end
        varargout{i}=C;
    end
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
