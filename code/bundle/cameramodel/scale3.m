function [Q,dQ,dQn]=scale3(P,k,varargin)
%SCALE3 3D isotropic scaling for the DBAT projection model.
%
%   Q=SCALE3(P,K) scales the 3D points in the 3-by-N array P by the
%   scalar K.
%
%   [Q,dQ]=... also returns a struct dQ with the analytical Jacobians
%   with respect to P and K in the fields dP and dK, respectively. For
%   more details, see DBAT_BUNDLE_FUNCTIONS.
%
%SEE ALSO: DBAT_BUNDLE_FUNCTIONS

% Treat selftest call separately.
if nargin>=1 && ischar(P), Q=selftest(nargin>1 && k); return; end

% Otherwise, verify number of parameters.
narginchk(2,4);

Q=[]; %#ok<NASGU>
dQ=[];
dQn=[];

if nargout>1
    % Construct empty Jacobian struct.
    dQ=struct('dP',[],...
              'dK',[]);
    dQn=dQ;
end

% What Jacobians to compute?
cP=nargout>1 && (length(varargin)<1 || varargin{1});
cK=nargout>1 && (length(varargin)<2 || varargin{2});

%% Test parameters
[m,~]=size(P);
if m~=3 || ~isscalar(k)
    error([mfilename,': bad size']);
end

%% Actual function code
Q=k*P;

if nargout>2
    %% Numerical Jacobian

    % FMT is function handle to repackage vector argument to what
    % the function expects.
    if cP
        fmt=@(P)reshape(P,3,[]);
        fun=@(P)feval(mfilename,fmt(P),k);
        dQn.dP=jacapprox(fun,P);
    end
    if cK
        fun=@(k)feval(mfilename,P,k);
        dQn.dK=jacapprox(fun,k);
    end
end

if nargout>1
    %% Analytical Jacobian
    if cP
        dQ.dP=k*speye(numel(P));
    end
    if cK
        dQ.dK=vec(P);
    end    
end


function fail=selftest(verbose)

% Set up test data.
n=3;
m=5;
k=rand+1;
P=rand(n,m);

fail=full_self_test(mfilename,{P,k},1e-8,1e-8,verbose);
