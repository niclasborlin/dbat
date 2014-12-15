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
%   ...=BUNDLE(S,...,TRACE), where TRACE is a string='trace' specifies
%   that the bundle should print an trace during the iterations.
%
%   ...=BUNDLE(S,...,SINGULARTEST), where SINGULARTEST is the
%   string='singulartest' specifies that the bundle should stop
%   immediately if a 'Matrix is singular' or 'Matrix is almost singular'
%   warning is issued on the normal matrix.
%
%   ...=BUNDLE(S,...,CHI), where CHI is a logical scalar, specifies if
%   chirality veto damping should be used (default: false). Chirality
%   veto damping is ignored for the undamped bundle.
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

% $Id$

maxIter=20;
damping='gna';
veto=false;
singularTest=false;
trace=false;

while ~isempty(varargin)
    if isnumeric(varargin{1}) && isscalar(varargin{1})
        % N
        maxIter=varargin{1};
        varargin(1)=[];
    elseif ischar(varargin{1})
        % DAMP, TRACE, or CXX
        switch lower(varargin{1})
          case {'none','gm','gna','lm','lmp','gh','lwb'}
            % OK
            damping=varargin{1};
            varargin(1)=[];
          case 'trace'
            trace=true;
            varargin(1)=[];
          case 'singulartest'
            singularTest=true;
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

% Create indices into the vector of unknowns. n is the number of unknowns.
[ixIO,ixEO,ixOP,n]=indvec([nnz(s.cIO),nnz(s.cEO),nnz(s.cOP)]);

% Set up vector of initial values.
x0=nan(n,1);
x0(ixIO)=s.IO(s.cIO);
x0(ixEO)=s.EO(s.cEO);
x0(ixOP)=s.OP(s.cOP);

% Residual and constraint functions.
resFun=@brown_euler_cam;
conFun=@brown_euler_cam_c;

if veto
    vetoFun=@chirality;
else
    vetoFun='';
end

% Convergence tolerance.
convTol=1e-3;
% Constraint tolerance.
constrTol=1e-8;

% Set up cell array of extra parameters.
params={s};

% For all optimization methods below, the final estimate is returned in x.
% The final residual vector and Jacobian are returned as r and J. Successive
% estimates of x and norm(r) are returned as columns of X and elements of
% res, respectively. Furthermore, a status code (0 - ok, -1 - too many
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

% Vector of Lagrange multipliers. Empty for unconstrained problems.
l=[];

switch lower(damping)
  case {'none','gm'}
    % Gauss-Markov with no damping.
    
    % Call Gauss-Markov optimization routine.
    stopWatch=cputime;
    [x,code,iters,r,J,X,res]=gauss_markov(resFun,x0,maxIter,convTol,trace, ...
                                      singularTest,params);
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
    [x,code,iters,r,J,X,res,alpha]=gauss_newton_armijo(resFun,vetoFun,x0, ...
                                                      maxIter,convTol,trace, ...
                                                      singularTest, ...
                                                      mu,alphaMin,params);  
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
    [x,code,iters,r,J,X,res,lambda]=levenberg_marquardt(resFun,vetoFun,x0, ...
                                                      maxIter,convTol,trace, ...
                                                      lambda0,lambdaMin,params);
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
    [x,code,iters,r,J,X,res,delta,rho,step]=levenberg_marquardt_powell(...
        resFun,vetoFun,x0,maxIter,convTol,trace, delta0,rhoBad,rhoGood,params);
    time=cputime-stopWatch;
    E.damping=struct('name','lmp','delta',delta,'rho',rho,'delta0',delta0,...
                     'rhoBad',rhoBad,'rhoGood',rhoGood,'step',step);
  case 'gh'
    % Gauss-Helmert constrained.

    % Armijo parameter.
    mu=-inf; % Always accept.
    nu0=0.1;
    % Shortest allowed step length.
    alphaMin=1e-9;
    
    % Call Gauss-Newton-Armijo optimization routine. The vector alpha is
    % returned with the step lengths used at each iteration.
    stopWatch=cputime;
    [x,code,iters,l,X,alpha,C,L,nus,r,J,A,res]=sqpsq(resFun,conFun,vetoFun,x0, ...
                                         maxIter,convTol,constrTol,trace, ...
                                         singularTest,mu,nu0,alphaMin,params);  
    time=cputime-stopWatch;
    E.damping=struct('name','gh','alpha',alpha,'mu',mu,'alphaMin',alphaMin);
  case 'lwb'
    % Lindstrom-Wedin-Borlin constrained.

    % Armijo parameter.
    mu=0.1;
    nu0=0.1;
    % Shortest allowed step length.
    alphaMin=1e-9;
    
    % Call Gauss-Newton-Armijo optimization routine. The vector alpha is
    % returned with the step lengths used at each iteration.
    stopWatch=cputime;
    [x,code,iters,l,X,alpha,C,L,nus,r,J,A,res]=sqpsq(resFun,conFun,vetoFun,x0, ...
                                         maxIter,convTol,constrTol,trace, ...
                                         singularTest,mu,nu0,alphaMin,params);  
    time=cputime-stopWatch;
    E.damping=struct('name','lwb','alpha',alpha,'mu',mu,'alphaMin',alphaMin);
  otherwise
    error('DBAT:bundle:internal','Unknown damping');
end

% Store iteration results.
E.res=res;
E.trace=X;
E.time=time;
E.code=code;
E.usedIters=iters;

% Store final residual and Jacobian for later covariance calculations.
E.J=J;
E.r=r;

% Handle returned values.
ok=code==0;

% Update s if optimization converged.
if ok
    s.IO(s.cIO)=x(ixIO);
    s.EO(s.cEO)=x(ixEO);
    s.OP(s.cOP)=x(ixOP);
end

% Sigma0 is sqrt(r'*r/(m-n)) in mm, convert to pixels.
s0mm=sqrt(r'*r/(length(r)-length(x)+length(l)));
s0px=s0mm*mean(s.IO(end-1:end));
s0=s0px;

E.s0mm=s0mm;
E.s0px=s0px;
