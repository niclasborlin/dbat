function [Q,dQ,dQn]=xlat2(U,c,varargin)
%XLAT2 2D translation for the DBAT projection model.
%
%   Q=XLAT2(U,C) translates the 2D points in the 2-by-N array U by
%   the 2-vector C.
%
%   [Q,dQ]=... also returns a struct dQ with the analytical Jacobians
%   with respect to U and C in the fields dU and dC, respectively. For
%   more details, see DBAT_BUNDLE_FUNCTIONS.
%
%SEE ALSO: DBAT_BUNDLE_FUNCTIONS

% Treat selftest call separately.
if nargin>=1 && ischar(U), Q=selftest(nargin>1 && c); return; end

% Otherwise, verify number of parameters.
narginchk(2,4);

Q=[]; %#ok<NASGU>
dQ=[];
dQn=[];

if nargout>1
    % Construct empty Jacobian struct.
    dQ=struct('dU',[],...
              'dC',[]);
    dQn=dQ;
end

% What Jacobians to compute?
cU=nargout>1 && (length(varargin)<1 || varargin{1});
cC=nargout>1 && (length(varargin)<2 || varargin{2});

%% Test parameters
[m,n]=size(U);
if m~=2 || any(size(c)~=[m,1])
    error([mfilename,': bad size']);
end

%% Actual function code
Q=U+repmat(c,1,n);

if nargout>2
    %% Numerical Jacobian

    % FMT is function handle to repackage vector argument to what
    % the function expects.
    if cU
        fmt=@(U)reshape(U,2,[]);
        fun=@(U)feval(mfilename,fmt(U),c);
        dQn.dU=jacapprox(fun,U);
    end
    if cC
        fun=@(c)feval(mfilename,U,c);
        dQn.dC=jacapprox(fun,c);
    end
end

if nargout>1
    %% Analytical Jacobian
    if cU
        dQ.dU=speye(numel(U));
    end
    if cC
        dQ.dC=repmat(speye(m),n,1);
    end    
end


function fail=selftest(verbose)

% Set up test data.
n=2;
m=5;
c=rand(n,1);
U=rand(n,m);

fail=full_self_test(mfilename,{U,c},1e-8,1e-8,verbose);
