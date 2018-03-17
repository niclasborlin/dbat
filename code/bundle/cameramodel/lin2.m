function [Q,dQ,dQn]=lin2(U,M,varargin)
%LIN2 2D linear transform for the DBAT projection model.
%
%   Q=LIN2(U,M) applies the linear transform M to the 2D points in the
%   2-by-N array U, i.e. computes M*U.
%
%   [Q,dQ]=... also returns a struct dQ with the analytical Jacobians
%   with respect to M and U in the fields dU and dM, respectively. For
%   more details, see DBAT_BUNDLE_FUNCTIONS.
%
%SEE ALSO: DBAT_BUNDLE_FUNCTIONS

% Treat selftest call separately.
if nargin>=1 && ischar(U), Q=selftest(nargin>1 && M); return; end

% Otherwise, verify number of parameters.
narginchk(2,4);

Q=[]; %#ok<NASGU>
dQ=[];
dQn=[];

if nargout>1
    % Construct empty Jacobian struct.
    dQ=struct('dU',[],...
              'dM',[]);
    dQn=dQ;
end

% What Jacobians to compute?
cU=nargout>1 && (length(varargin)<1 || varargin{1});
cM=nargout>1 && (length(varargin)<2 || varargin{2});

%% Test parameters
[m,n]=size(U);
if m~=2 || any(size(M)~=m)
    error([mfilename,': bad size']);
end

%% Actual function code
Q=M*U;

if nargout>2
    %% Numerical Jacobian

    % FMT is function handle to repackage vector argument to what
    % the function expects.
    if cU
        fmt=@(U)reshape(U,2,[]);
        fun=@(U)feval(mfilename,fmt(U),M);
        dQn.dU=jacapprox(fun,U);
    end
    if cM
        fmt=@(M)reshape(M,2,[]);
        fun=@(M)feval(mfilename,U,fmt(M));
        dQn.dM=jacapprox(fun,M);
    end
end

if nargout>1
    %% Analytical Jacobian
    if cU
        dQ.dU=kron(speye(n),M);
    end
    if cM
        dQ.dM=kron(U',speye(m));
    end    
end


function fail=selftest(verbose)

% Set up test data.
m=2;
n=5;
M=rand(m);
U=rand(m,n);

fail=full_self_test(mfilename,{U,M},1e-8,1e-8,verbose);
