function [Q,dQ,dQn]=affscale2(U,k,b1,b2,varargin)
%AFFSCALE2 2D affine scaling for the DBAT projection model.
%
%   Q=AFFSCALE2(U,K,B1,B2) scales the 2D points in the 2-by-N array U
%   by the scalar K and applies an affine transformation. The affine
%   transformation matrix is [1+B1 B2;0 1].
%
%   [Q,dQ]=... also returns a struct dQ with the analytical Jacobians
%   with respect to U, K, B1, and B2 in the fields dU, dK, dB1, and
%   dB2, respectively. For more details, see DBAT_BUNDLE_FUNCTIONS.
%
%SEE ALSO: AFFINE2, DBAT_BUNDLE_FUNCTIONS

% Treat selftest call separately.
if nargin>=1 && ischar(U), Q=selftest(nargin>1 && k); return; end

% Otherwise, verify number of parameters.
narginchk(4,8);

Q=[]; %#ok<NASGU>
dQ=[];
dQn=[];

if nargout>1
    % Construct empty Jacobian struct.
    dQ=struct('dU',[],...
              'dB1',[],...
              'dB2',[],...
              'dK',[]);
    dQn=dQ;
end

% What Jacobians to compute?
cU=nargout>1 && (length(varargin)<1 || varargin{1});
cK=nargout>1 && (length(varargin)<2 || varargin{2});
cB1=nargout>1 && (length(varargin)<3 || varargin{3});
cB2=nargout>1 && (length(varargin)<4 || varargin{4});

%% Test parameters
[m,~]=size(U);
if m~=2 || ~isscalar(k) || ~isscalar(b1) || ~isscalar(b2)
    error([mfilename,': bad size']);
end

%% Actual function code
[A,dA]=affine2(U,b1,b2);
Q=k*A;

if nargout>2
    %% Numerical Jacobian

    % FMT is function handle to repackage vector argument to what
    % the function expects.
    if cU
        fmt=@(U)reshape(U,2,[]);
        fun=@(U)feval(mfilename,fmt(U),k,b1,b2);
        dQn.dU=jacapprox(fun,U);
    end
    if cK
        fun=@(k)feval(mfilename,U,k,b1,b2);
        dQn.dK=jacapprox(fun,k);
    end
    if cB1
        fun=@(b1)feval(mfilename,U,k,b1,b2);
        dQn.dB1=jacapprox(fun,b1);
    end
    if cB2
        fun=@(b2)feval(mfilename,U,k,b1,b2);
        dQn.dB2=jacapprox(fun,b2);
    end
end

if nargout>1
    %% Analytical Jacobian
    if cU
        dQ.dU=k*dA.dU;
    end
    if cK
        dQ.dK=A(:);
    end    
    if cB1
        dQ.dB1=k*dA.dB1;
    end
    if cB2
        dQ.dB2=k*dA.dB2;
    end
end


function fail=selftest(verbose)

% Set up test data.
n=2;
m=5;
k=rand+1;
b1=0.1+rand;
b2=0.1+rand;
U=rand(n,m);

fail=full_self_test(mfilename,{U,k,b1,b2},1e-8,1e-8,verbose);
