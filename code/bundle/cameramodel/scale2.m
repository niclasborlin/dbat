function [Q,dQ,dQn]=scale2(k,P,varargin)
%SCALE2 2D scaling for the DBAT projection model
%
%   Q=SCALE2(K,P) scales the 2D points in the 2-by-N array P by the
%   scalar K.
%
%   [Q,dQ]=... also returns a struct dQ with the analytical Jacobians
%   with respect to P and K in the fields dP and dK. For more details,
%   see DBAT_BUNDLE_FUNCTIONS.
%
%SEE ALSO: DBAT_BUNDLE_FUNCTIONS

% Treat selftest call separately.
if nargin>=1 && ischar(k), Q=selftest(nargin>1 && P); return; end

% Otherwise, verify number of parameters.
narginchk(1,4);

Q=[];
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
[m,n]=size(P);
if m~=2 || ~isscalar(k)
    error([mfilename,': bad size']);
end

%% Actual function code
Q=k*P;

if nargout>2
    %% Numerical Jacobian

    % FMT is function handle to repackage vector argument to what
    % the function expects.
    if cP
        fmt=@(P)reshape(P,2,[]);
        fun=@(P)feval(mfilename,k,fmt(P));
        dQn.dP=jacapprox(fun,P);
    end
    if cK
        fun=@(k)feval(mfilename,k,P);
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
n=2;
m=5;
k=rand+1;
P=rand(n,m);

fail=full_self_test(mfilename,{k,P},1e-8,1e-8,verbose);
