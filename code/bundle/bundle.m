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

E=struct('damping',damping);

switch lower(damping)
  case {'none','gm'}
    % Gauss-Markov with no damping.
    
    % Call Gauss-Markov optimization routine.
    stopWatch=cputime;
    [x,code,iters,f,J,X]=gauss_markov(resFun,x0,maxIter,convTol,trace, ...
                                      singular_test,params);
    time=cputime-stopWatch;
  case 'gna'
    % Gauss-Newton with Armijo linesearch.

    % Armijo parameter.
    mu=0.1;
    % Shortest allowed step length.
    alphaMin=1e-3;
    
    % Call Gauss-Newton-Armijo optimization routine. The vector alpha is
    % returned with the step lengths used at each iteration.
    stopWatch=cputime;
    [x,code,iters,f,J,X,alpha]=gauss_newton_armijo(resFun,vetoFun,x0, ...
                                                   mu,alphaMin,maxIter, ...
                                                   convTol,trace, ...
                                                   singular_test,params);
    time=cputime-stopWatch;
    E.alpha=alpha;
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
  otherwise
    error('DBAT:bundle:internal','Unknown damping');
end
E.trace=X;
E.time=time;

% Handle returned values.
ok=code==0;

% Update s if optimization converged.
if ok
    s.IO(s.cIO)=x0(ixIO);
    s.EO(s.cEO)=x0(ixEO);
    s.OP(s.cOP)=x0(ixOP);
end

% s0=sqrt(f'*f/(m-n)) in mm, convert to pixels.
s0=sqrt(f'*f/(length(f)-length(x)))*mean(s.IO(end-1:end));

covMatrices={'ceof'};
% Compute covariance matrices.
if ~isempty(covMatrices)
    % We may need J'*J many times, precalculate and prefactor.
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
    
    % Inverse permutation.
    invP=zeros(size(p));
    invP(p)=1:length(p);
    
    R=chol(JTJ(p,p));
    RT=R';
    % Do a split of R into [R1, R12;
    %                        0   R2]
    % R1 will generally be block-diagonal - use sparse storage.
    % R12, R2 will generally be dense - use full storage.
    
    % Find jump in column density.
    colDens=full(sum(R~=0,1))/size(R,1);
    [dummy,i]=max(diff(colDens));
    ix1=1:i;
    ix2=i+1:size(R,1);
    R1=R(ix1,ix1);
    R12s=R(ix1,ix2);
    R12=full(R12s);
    R2s=R(ix2,ix2);
    R2=full(R2s);
    
    % Permute the identity to match permutation of JTJ.
    Ip=speye(size(JTJ,2));
    Ip=Ip(p,:);
    for i=1:length(covMatrices)
        switch covMatrices{i}
          case 'cxx' % Raw, whole covariance matrix.

            C=zeros(size(JTJ));
            C(p,:)=R\(R'\Ip);
            
          case 'ciof' % Whole CIO covariance matrix.
            
            % Pre-allocate matrix with place for nnz(s.cIO)^2 elements.
            C=spalloc(numel(s.IO),numel(s.IO),nnz(s.cIO)^2);
            
            % Compute needed part of inverse.
            cc=zeros(size(JTJ,1),length(ixIO));
            cc(p,:)=R\(R'\Ip(:,ixIO));
            
            % We get the full columns of the inverse. Extract only the
            % interesting block and put it into the right part of C.
            C(s.cIO(:),s.cIO(:))=cc(ixIO,:);
            
          case 'ceof' % Whole CEO covariance matrix.
            
            start=clock;
            % Pre-allocate matrix with place for nnz(s.cEO)^2 elements.
            C=spalloc(numel(s.EO),numel(s.EO),nnz(s.cEO)^2);
            
            % Compute needed part of inverse.
            cc=zeros(size(JTJ,1),length(ixEO));
            cc(p,:)=R\(R'\Ip(:,ixEO));
            
            % We get the full columns. Extract only the interesting block
            % and put it into the right part of C.
            C(s.cEO(:),s.cEO(:))=cc(ixEO,:);

            etime(clock,start)
            
          case 'copf' % Whole COP covariance matrix.
            
            % Pre-allocate matrix with place for nnz(s.cOP)^2 elements.
            C=spalloc(numel(s.OP),numel(s.OP),nnz(s.cOP)^2);
            
            % Compute needed part of inverse.
            cc=zeros(size(JTJ,1),length(ixOP));
            cc(p,:)=R\(R'\Ip(:,ixOP));
            
            % We get the full columns. Extract only the interesting block
            % and put it into the right part of C.
            C(s.cOP(:),s.cOP(:))=cc(ixOP,:);
            
          case 'cio' % Block-diagonal CIO

            % Delayed progress dialog.
            start=clock;
            h=[];
            
            % Pre-allocate matrix with place for the diagonal blocks.
            C=spalloc(numel(s.IO),numel(s.IO),sum(sum(s.cIO,1).^2));
                     
            % Pack indices as the data
            ix=zeros(size(s.IO));
            ix(s.cIO)=ixIO;
            
            % Construct corresponding indexing among the IO parameters.
            ixInt=reshape(1:numel(s.IO),size(s.IO));
            
            % Loop over each IO column.
            for j=1:size(ix,2)
                % Indices into J.
                jix=ix(s.cIO(:,j),j);
                % Indices into s.IO
                eix=ixInt(s.cIO(:,j),j);

                % Compute needed part of inverse.
                cc=zeros(size(JTJ,1),length(jix));
                cc(p,:)=R\(R'\Ip(:,jix));
                % We get the full columns. Extract only the interesting block
                % and put it into the right part of C.
                C(eix,eix)=cc(jix,:);
                
                if (isempty(h) && etime(clock,start)>1)
                    % Create dialog.
                    h=waitbar(j/size(ix,2),'Computing IO covariances');
                elseif rem(j,floor(sqrt(size(ix,2))))==0
                    % Update dialog.
                    if ishandle(h) % Guard against window close.
                        waitbar(j/size(ix,2),h);
                    end
                end
            end
            if ishandle(h), close(h), end
            
          case 'ceo' % Block-diagonal CEO
            
            % Delayed progress dialog.
            start=clock;
            lapTime=start;
            h=[];
            
            % Pre-allocate matrix with place for the diagonal blocks.
            C=spalloc(numel(s.EO),numel(s.EO),sum(sum(s.cEO,1).^2));
                     
            % Pack indices as the data
            ix=zeros(size(s.EO));
            ix(s.cEO)=ixEO;
            
            % Construct corresponding indexing among the EO parameters.
            ixInt=reshape(1:numel(s.EO),size(s.EO));

            % Determine block column size such that computed part of inverse is
            % approximately 10M-elements.
            bsElems=10*1024^2;
            bsCols=floor(bsElems/size(JTJ,1)/max(sum(ix~=0,1)));
            bsCols=min(max(bsCols,1),size(ix,2));
            
            % Loop over each EO column.
            for j=1:bsCols:size(ix,2)
                % Columns in block.
                jCols=j:min(j+bsCols-1,size(ix,2));
                
                % Indices into J.
                jixBlock=ix(:,jCols);
                jix=jixBlock(s.cEO(:,jCols));
                % Indices into s.EO
                eixBlock=ixInt(:,jCols);
                eix=eixBlock(s.cEO(:,jCols));

                % Compute needed part of inverse.
                cc=zeros(size(JTJ,1),length(jix));
                cc(p,:)=R\(R'\Ip(:,jix));
                % We get the full columns of the inverse. Extract only
                % the interesting block...
                Cblock=spalloc(size(C,1),size(C,1),length(eix)^2);
                Cblock(eix,eix)=cc(jix,:);
                % ...make it block-diagonal...
                Cblock=mkblkdiag(Cblock,size(ix,1));
                % ...and put it into the right part of C.
                C(eix,eix)=Cblock(eix,eix);
                
                if (isempty(h) && etime(clock,start)>1) && ...
                        jCols(end)~=size(ix,2)
                    % Only create dialog if execution takes more than 1s
                    % and this iteration is not the last.
                    h=waitbar(jCols(end)/size(ix,2),'Computing EO covariances');
                    lapTime=clock;
                elseif etime(clock,lapTime)>1
                    % Update dialog.
                    if ishandle(h) % Guard against window close.
                        waitbar(jCols(end)/size(ix,2),h);
                    end
                    lapTime=clock;
                end
            end
            if ishandle(h), close(h), end
            etime(clock,start)
            %            max(max(abs(mkblkdiag(varargout{i-1},size(ix,1))-C)))
            
          case 'cop' % Block-diagonal COP
            
            % Delayed progress dialog.
            start=clock;
            lapTime=start;
            h=[];

            % Pre-allocate matrix with place for the diagonal blocks.
            C=spalloc(numel(s.OP),numel(s.OP),sum(sum(s.cOP,1).^2));
                     
            % Pack indices as the data
            ix=zeros(size(s.OP));
            ix(s.cOP)=ixOP;
            
            % Construct corresponding indexing among the OP parameters.
            ixInt=reshape(1:numel(s.OP),size(s.OP));
            
            % Determine block column size such that computed part of inverse is
            % approximately 10M-elements.
            bsElems=50*1024^2;
            bsCols=floor(bsElems/size(JTJ,1)/max(sum(ix~=0,1)))
            bsCols=min(max(bsCols,1),size(ix,2));
            
            % Loop over each OP column.
            for j=1:bsCols:size(ix,2)
                % Columns in block.
                jCols=j:min(j+bsCols-1,size(ix,2));

                % Indices into J.
                jixBlock=ix(:,jCols);
                jix=jixBlock(s.cOP(:,jCols));
                % Indices into s.OP
                eixBlock=ixInt(:,jCols);
                eix=eixBlock(s.cOP(:,jCols));

                % Compute needed part of inverse.
                cc=zeros(size(JTJ,1),length(jix));
                % Original command:
                %cc(p,:)=R\(R'\Ip(:,jix));
                
                % Split computation into several steps:
                tmp1=Ip(:,jix);
                tmp2=RT\tmp1;
                cc0pt2=R2\tmp2(ix2,:);
                tmp3s=R12s*cc0pt2;
                tmp4=tmp2(ix1,:)-tmp3s;
                cc0pt1=R1\tmp4;
                
                cc(p,:)=[cc0pt1;cc0pt2];
                
                % Compare actual computation with slow result.
                % max(max(abs([cc0pt1;cc0pt2]-R\(R'\Ip(:,jix)))))
                
                % We get the full columns of the inverse. Extract only
                % the interesting block...
                Cblock=spalloc(size(C,1),size(C,1),length(eix)^2);
                Cblock(eix,eix)=cc(jix,:);
                % ...make it block-diagonal...
                Cblock=mkblkdiag(Cblock,size(ix,1));
                % ...and put it into the right part of C.
                C(eix,eix)=Cblock(eix,eix);
                
                if (isempty(h) && etime(clock,start)>1) && ...
                        jCols(end)~=size(ix,2)
                    % Only create dialog if execution takes more than 1s
                    % and this iteration is not the last.
                    h=waitbar(jCols(end)/size(ix,2),'Computing OP covariances');
                    lapTime=clock;
                elseif etime(clock,lapTime)>1
                    % Update dialog.
                    if ishandle(h) % Guard against window close.
                        waitbar(jCols(end)/size(ix,2),h);
                    end
                    lapTime=clock;
                end
            end
            if ishandle(h), close(h), end
            etime(clock,start)
        end
        varargout{i}=C;
    end
end
