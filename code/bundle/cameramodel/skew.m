function [Q,dQ,dQn]=skew(U,k,varargin)
%SKEW 2D anisotropic scaling for the DBAT projection model.
%
%   Q=SKEW(U,K) applies a skew to the 2D points in the 2-by-N array
%   U. The actual transformation is Q=[1,k;0,1]*U.
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
[m,n]=size(U);
if m~=2 || ~isscalar(k)
    error([mfilename,': bad size']);
end

%% Actual function code
Q=[1,k;0,1]*U;

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
        % Each block is A = [1, K; 0, 1].

        % Values for Jacobian in memory-order.
        vv=repmat([1,k,1]',1,n);
        % Row indices.
        ii=repmat(0:2:2*n-1,3,1)+repmat([1,1,2]',1,n);
        % Column indices.
        jj=repmat(0:2:2*n-1,3,1)+repmat([1,2,2]',1,n);
        dQ.dU=sparse(ii,jj,vv,2*n,2*n);
    end
    if cK
        dQ.dK=zeros(2*n,1);
        dQ.dK(1:2:end)=U(2,:)';
    end    
end


function fail=selftest(verbose)

% Set up test data.
n=2;
m=5;
k=rand+0.01;
U=rand(n,m);

fail=full_self_test(mfilename,{U,k},1e-8,1e-8,verbose);
