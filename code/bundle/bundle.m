function [s,X]=bundle(s,varargin)
%BUNDLE Run bundle adjustment iterations on a camera network.
%
%   [S,OK]=BUNDLE(S), where S is a struct returned by PROB2DBATSTRUCT, runs
%   the damped bundle adjustment on the camera network in the structure
%   S. The parameter values in S are used as initial values. The cIO, cEO,
%   cOP fields of S are used to indicate which parameters are free. OK is
%   returned as true if the bundle converged within the allowed number of
%   iterations. On return, the parameter values in S are updated if the
%   bundle converged.
%
%   ...=BUNDLE(S,...,N), where N is an integer, sets the maximum number of
%   iterations to N (default: 20).
%
%   ...=BUNDLE(S,...,DAMP), where DAMP is a string, specifies which damping
%   to use: 'none' or 'GM' (classic bundle with no damping), 'GNA'
%   (Gauss-Newton with Armijo linesearch, default), 'LM' (original
%   Levenberg-Marquardt) , 'LMP' (Levenberg-Marquardt with Powell dogleg).
%
%   ...=BUNDLE(S,...,CHI), where CHI is a logical scalar, specifies if
%   chirality veto damping should be used (default: false). Chirality
%   veto damping is ignored for the undamped bundle.
%
%   [S,OK,S0]=... returns the sigma0 for the last iteration.
%
%   [S,OK,S0,X]=... returns successive parameter estimates as columns in X.
%
%   [S,OK,S0,X,CXX]=... returns the covariance matrix CXX of the final X,
%   scaled by sigma0.
%
%References: Börlin, Grussenmeyer (2013), "Bundle Adjustment with and
%   without Damping", Photogrammetric Record 28(144), pp. XX-YY. DOI
%   10.1111/phor.12037.

% $Id$

maxIter=20;
damp='gna';
veto=false;

while ~isempty(varargin)
    if isnumeric(varargin{1}) && isscalar(varargin{1})
        % N
        maxIter=varargin{1};
        varargin(1)=[];
    elseif ischar(varargin{1})
        % DAMP
        switch lower(varargin{1})
          case {'none','gm','gna','lm','lmp'}
            % OK
            damp=varargin{1};
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
ixIO=find(s.cIO);
ixEO=nnz(ixIO)+find(s.cEO);
ixOP=nnz(ixIO)+nnz(ixEO)+find(s.cOP);

% Number of unknowns.
n=nnz(ixIO)+nnz(ixEO)+nnz(ixOP);

% Set up vector of initial values.
x0=nan(n,1);
x0(ixIO)=s.IO(s.cIO);
x0(ixEO)=s.EO(s.cEO);
x0(ixOP)=s.OP(s.cOP);

% Residual function.
resFun=@pm_eulerbundle1p_f;

if veto
    vetoFun=@pm_chirality;
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

switch lower(damp)
  case {'none','gm'}
    % Gauss-Markov with no damping.
    
    % Call Gauss-Markov optimization routine.
    stopWatch=cputime;
    [x,code,iters,X,f,J]=gauss_markov(resFun,x0,convTol,maxIter,params);
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
    [x,code,iters,X,f,J,alpha]=gauss_newton_armijo(resFun,vetoFun,x0, ...
                                                   alphaMin,convTol, ...
                                                   maxIter,params);
    time=cputime-stopWatch;
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
    [x,code,iters,X,f,J,lambdas]=levenberg_marquardt(resFun,vetoFun,x0, ...
                                                     f0,J0,lambda0, ...
                                                     lambdaMin,params);
    time=cputime-stopWatch;
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
    [x,code,iters,X,f,J,deltas,rhos]=levenberg_marquardt_powell(resFun, ...
                                                      vetoFun,x0,delta0, ...
                                                      rhoBad,rhoGood,params);
    time=cputime-stopWatch;
  otherwise
    error('DBAT:bundle:internal','Unknown damping');
end

% Handle returned values...