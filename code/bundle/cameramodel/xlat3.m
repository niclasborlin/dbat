function [Q,dQ,dQn]=xlat3(P,c,varargin)
%XLAT3 3D translation for the DBAT projection model.
%
%   Q=XLAT3(P,C) translates the 3D points in the 3-by-N array P by
%   the 3-vector C.
%
%   [Q,dQ]=... also returns a struct dQ with the analytical Jacobians
%   with respect to P and C in the fields dP and dC, respectively. For
%   more details, see DBAT_BUNDLE_FUNCTIONS.
%
%SEE ALSO: DBAT_BUNDLE_FUNCTIONS

% Treat selftest call separately.
if nargin>=1 && ischar(P), Q=selftest(nargin>1 && c); return; end

% Otherwise, verify number of parameters.
narginchk(2,4);

Q=[]; %#ok<NASGU>
dQ=[];
dQn=[];

if nargout>1
    % Construct empty Jacobian struct.
    dQ=struct('dP',[],...
              'dC',[]);
    dQn=dQ;
end

% What Jacobians to compute?
cP=nargout>1 && (length(varargin)<1 || varargin{1});
cC=nargout>1 && (length(varargin)<2 || varargin{2});

%% Test parameters
[m,n]=size(P);
if m~=3 || any(size(c)~=[m,1])
    error([mfilename,': bad size']);
end

%% Actual function code
Q=P+repmat(c,1,n);

if nargout>2
    %% Numerical Jacobian

    % FMT is function handle to repackage vector argument to what
    % the function expects.
    if cP
        fmt=@(P)reshape(P,3,[]);
        fun=@(P)feval(mfilename,fmt(P),c);
        dQn.dP=jacapprox(fun,P);
    end
    if cC
        fun=@(c)feval(mfilename,P,c);
        dQn.dC=jacapprox(fun,c);
    end
end

if nargout>1
    %% Analytical Jacobian
    if cP
        dQ.dP=speye(numel(P));
    end
    if cC
        dQ.dC=repmat(speye(m),n,1);
    end    
end


function fail=selftest(verbose)

% Set up test data.
n=3;
m=5;
c=rand(n,1);
P=rand(n,m);

fail=full_self_test(mfilename,{P,c},1e-8,1e-8,verbose);
