function [s,ok,iters,s0,E]=bundle(s,varargin)
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
%   ...=BUNDLE(S,...,'trace') specifies that the bundle should
%   print an trace during the iterations.
%
%   ...=BUNDLE(S,...,'singulartest') specifies that the bundle should
%   stop immediately if a 'Matrix is singular' or 'Matrix is almost
%   singular' warning is issued on the normal matrix.
%
%   ...=BUNDLE(S,...,CHI), where CHI is a logical scalar, specifies if
%   chirality veto damping should be used (default: false). Chirality
%   veto damping is ignored for the undamped bundle.
%
%   ...=BUNDLE(S,...,'pmdof') uses Photomodeler degrees-of-freedom
%   computation.
%
%   ...=BUNDLE(S,...,'dofverb') outputs how the degrees of freedom
%   are calculated.
%
%   [S,OK,N,S0]=... returns the sigma0 (in pixel units) for the last iteration.
%
%   [S,OK,N,S0,E]=... returns a struct E with information about the iterations:
%       E.trace   - NOBS-by-(N+1) array with successive parameter estimates.
%       E.res     - (N+1)-vector with the residual norm at every iteration.
%       E.damping - struct with damping-specific information, including
%           E.damping.name - the name of the damping scheme used.
%
%   Use BUNDLE_COV to compute covariances, etc., of the result or
%   BUNDLE_RESULT_FILE to generate a result file.
%
%   References:
%     Börlin, Grussenmeyer (2013), "Bundle Adjustment With and Without
%       Damping". Photogrammetric Record 28(144), pp. 396-415. DOI
%       10.1111/phor.12037.
%
%See also: GAUSS_MARKOV, GAUSS_NEWTON_ARMIJO, LEVENBERG_MARQUARDT,
%   LEVENBERG_MARQUARDT_POWELL, BROWN_EULER_CAM, BUNDLE_COV,
%   PROB2DBATSTRUCT, BUNDLE_RESULT_FILE


maxIter=20;
damping='gna';
veto=false;
singularTest=false;
trace=false;
dofVerb=false;
pmDof=false;

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
          case 'singulartest'
            singularTest=true;
            varargin(1)=[];
          case 'pmdof'
            pmDof=true;
            varargin(1)=[];
          case 'dofverb'
            dofVerb=true;
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

% Verify a consistent estimation/use observations arrays. If we
% shouldn't estimate a parameter, we should not use the prior
% observation of it during the bundle.
if any(s.useIOobs(~s.estIO))
    error('DBAT:bundle:badInput','IO estimate/use prior obs mismatch');
end
if any(s.useEOobs(~s.estEO))
    error('DBAT:bundle:badInput','EO estimate/use prior obs mismatch');
end
if any(s.useOPobs(~s.estOP))
    error('DBAT:bundle:badInput','OP estimate/use prior obs mismatch');
end

% Create indices into the vector of unknowns. n is the total number of unknowns.
[ixIO,ixEO,ixOP,n]=indvec([nnz(s.estIO),nnz(s.estEO),nnz(s.estOP)]);

% Set up vector of initial values.
x0=nan(n,1);
x0(ixIO)=s.IO(s.estIO);
x0(ixEO)=s.EO(s.estEO);
x0(ixOP)=s.OP(s.estOP);

% Residual function.
resFun=@(x)brown_euler_cam(x,s);

if veto
    vetoFun=@chirality;
else
    vetoFun='';
end

% Covariance estimates for observations.

% Mark points. Standard deviation is given in pixel units, the
% residuals are in mm, so scale the variance.
ptCols=s.colPos(s.vis);
varMark=(s.markStd(:,ptCols)./s.IO(end-1:end,s.ptCams(ptCols))).^2;
% Variance of prior IO observations.
varIO=s.prior.IOstd(s.useIOobs).^2;
% Prior EO observations.
varEO=s.prior.EOstd(s.useEOobs).^2;
% Prior OP observations.
varOP=s.prior.OPstd(s.useOPobs).^2;

% Create indices into the residual vector. nObs is the total number
% of observations.
[resIxMarkPt,resIxIO,resIxEO,resIxOP,nObs]=indvec(...
    [numel(varMark),numel(varIO),numel(varEO),numel(varOP)]);

% Weight matrix.
varAll=[varMark(:);varIO(:);varEO(:);varOP(:)];
Cobs=spdiags(varAll,0,nObs,nObs);
% Use the inverse as the weight matrix.
W=inv(Cobs);

% Convergence tolerance.
convTol=1e-6;

% For all optimization methods below, the final estimate is returned
% in x.  The final weighted and unweighted residual vectors and
% Jacobians are returned in the struct final. Successive estimates of x and
% norm(r) are returned as columns of X and elements of res,
% respectively. Furthermore, a status code (0 - ok, -1 - too many
% iterations) and the number of required iterations are returned.

% Set up return struct with bundle setup.

% TODO: Verify consistent return sizes of trace array, residual
% estimates, and damping et al. vectors given
% 1) x* = x0
% 2) x* found before maxIter
% 3) x* found at maxIter
% 4) x* not found at maxIter

% Version string.
[v,d]=dbatversion;
E=struct('maxIter',maxIter,'convTol',convTol,'singularTest',singularTest,...
         'chirality',veto,'dateStamp',datestr(now),...
         'version',sprintf('%s (%s)',v,d));

switch lower(damping)
  case {'none','gm'}
    % Gauss-Markov with no damping.
    
    % Call Gauss-Markov optimization routine.
    stopWatch=cputime;
    [x,code,iters,final,X,res]=gauss_markov(resFun,x0,W,maxIter, ...
                                          convTol,trace, singularTest);
    time=cputime-stopWatch;
    E.damping=struct('name','gm');
  case 'gna'
    % Gauss-Newton with Armijo linesearch.

    % Armijo parameter.
    mu=0.1;
    % Shortest allowed step length.
    alphaMin=1e-9;
    
    % Call Gauss-Newton-Armijo optimization routine. The vector alpha is
    % returned with the step lengths used at each iteration.
    stopWatch=cputime;
    [x,code,iters,final,X,res,alpha]=gauss_newton_armijo(resFun, ...
                                                      vetoFun,x0,W, ...
                                                      maxIter, ...
                                                      convTol,trace, ...
                                                      singularTest, ...
                                                      mu,alphaMin);  
    time=cputime-stopWatch;
    E.damping=struct('name','gna','alpha',alpha,'mu',mu,'alphaMin',alphaMin);
  case 'lm'
    % Original Levenberg-Marquardt "lambda"-version.

    % Initial lambda value. A negative value tells levenberg_marquardt to
    % scale lambda0 by trace(J0'*J0)/size(J0,2).
    lambda0=-1e-10;

    % Use lambda0 as lower damping cutoff value.
    lambdaMin=lambda0;

    % Call Levenberg-Marquardt optimization routine. The vector lambdas is
    % returned with the lambda values used at each iteration.
    stopWatch=cputime;
    [x,code,iters,final,X,res,lambda]=levenberg_marquardt(resFun,vetoFun,x0, ...
                                                      W,maxIter, ...
                                                      convTol,trace, ...
                                                      lambda0,lambdaMin);
    time=cputime-stopWatch;
    E.damping=struct('name','lm','lambda',lambda,'lambda0',lambda(1), ...
                     'lambdaMin',lambda(1));
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
    [x,code,iters,final,X,res,delta,rho,step]=levenberg_marquardt_powell(...
        resFun,vetoFun,x0,W,maxIter,convTol,trace, delta0,rhoBad,rhoGood);
    time=cputime-stopWatch;
    E.damping=struct('name','lmp','delta',delta,'rho',rho,'delta0',delta0,...
                     'rhoBad',rhoBad,'rhoGood',rhoGood,'step',step);
  otherwise
    error('DBAT:bundle:internal','Unknown damping');
end

% Store iteration results.
E.res=res;
E.trace=X;
E.time=time;
E.code=code;
E.usedIters=iters;

% Store final weighted residual and Jacobian for later covariance calculations.
E.final=final;

% Handle returned values.
ok=code==0;

% Update s if optimization converged.
if ok
    s.IO(s.estIO)=x(ixIO);
    s.EO(s.estEO)=x(ixEO);
    s.OP(s.estOP)=x(ixOP);
end

% Always update the residuals.
s.residuals.markPt(:)=final.unweighted.r(resIxMarkPt);
% Mark pt residuals are in mm, scale to pixels.
s.residuals.markPt=s.residuals.markPt.*s.IO(end-1:end,s.ptCams(ptCols));

s.residuals.IO(s.useIOobs)=final.unweighted.r(resIxIO);
s.residuals.EO(s.useEOobs)=final.unweighted.r(resIxEO);
s.residuals.OP(s.useOPobs)=final.unweighted.r(resIxOP);

% Sigma0 is sqrt(r'*r/(m-n)), where m is the number of
% observations, and n is the number of unknowns.

% Extra observations without a specific reisdual (fixed control points
% that have been measured, fixed camera stations that have been
% used).
if pmDof
    p=nnz(s.estOP(:,any(s.vis,2))==0)+nnz(s.estEO(1:6,any(s.vis,1))==0);
else
    p=0;
end

r=E.final.weighted.r;
lenR=length(r);
lenX=length(x);
dof=(lenR+p-lenX);
if dofVerb
    fprintf('%s: dof=%d+%d-%d=%d.\n',mfilename,lenR,p,lenX,dof);
end
s0=sqrt((r'*r)/dof);
sigmas=s0*s.prior.sigmas;

E.numObs=lenR;
E.numParams=lenX;
E.redundancy=dof;
E.s0=s0;
E.sigmas=sigmas;
