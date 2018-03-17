function [Q,dQ,dQn]=scale2(U,k,varargin)
%SCALE2 2D isotropic scaling for the DBAT projection model.
%
%   Q=SCALE2(U,K) scales the 2D points in the 2-by-N array U by the
%   scalar K.
%
%   [Q,dQ]=... also returns a struct dQ with the analytical Jacobians
%   with respect to U and K in the fields dU and dK, respectively. For
%   more details, see DBAT_BUNDLE_FUNCTIONS.
%
%SEE ALSO: DBAT_BUNDLE_FUNCTIONS

% Treat selftest call separately.
if nargin>=1 && ischar(U), Q=selftest(nargin>1 && k); return; end

% Otherwise, verify number of parameters.
narginchk(2,4);

Q=[]; %#ok<NASGU>
dQ=[];
dQn=[];

if nargout>1
    % Construct empty Jacobian struct.
    dQ=struct('dU',[],...
              'dK',[]);
    dQn=dQ;
end

% What Jacobians to compute?
cU=nargout>1 && (length(varargin)<1 || varargin{1});
cK=nargout>1 && (length(varargin)<2 || varargin{2});

%% Test parameters
[m,~]=size(U);
if m~=2 || ~isscalar(k)
    error([mfilename,': bad size']);
end

%% Actual function code
Q=k*U;

if nargout>2
    %% Numerical Jacobian

    % FMT is function handle to repackage vector argument to what
    % the function expects.
    if cU
        fmt=@(U)reshape(U,2,[]);
        fun=@(U)feval(mfilename,fmt(U),k);
        dQn.dU=jacapprox(fun,U);
    end
    if cK
        fun=@(k)feval(mfilename,U,k);
        dQn.dK=jacapprox(fun,k);
    end
end

if nargout>1
    %% Analytical Jacobian
    if cU
        dQ.dU=k*speye(numel(U));
    end
    if cK
        dQ.dK=vec(U);
    end    
end


function fail=selftest(verbose)

% Set up test data.
n=2;
m=5;
k=rand+1;
U=rand(n,m);

fail=full_self_test(mfilename,{U,k},1e-8,1e-8,verbose);
