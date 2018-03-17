function [v,dv,dvn]=power_vec(x,nn,varargin)
%POWER_VEC Vector of powers function used by DBAT lens distortion functions.
%
%   V=POWER_VEC(X,N) returns the vector of powers of the scalar X,
%   V=[X, X^2, ..., X^N]'. If X is an M-vector, V will be an N-by-M
%   array.
%
%   [V,dV]=... also returns a struct dV with the analytical Jacobian
%   with respect to X in the field dX. For more details, see
%   DBAT_BUNDLE_FUNCTIONS.
%
%SEE ALSO: DBAT_BUNDLE_FUNCTIONS

% Treat selftest call separately.
if nargin>=1 && ischar(x), v=selftest(nargin>1 && varargin{1}); return; end

% Otherwise, verify number of parameters.
% Actually, in production we only want 3 parameters, but the full
% testing we may get 4.
narginchk(1,4);

v=[]; %#ok<NASGU>
dv=[];
dvn=[];

if nargout>1
    % Construct empty Jacobian struct.
    dv=struct('dX',[]);
    dvn=dv;
end

% What Jacobians to compute?
cX=nargout>1 && (length(varargin)<1 || varargin{1});

%% Test parameters
[m,n]=size(x);
if m~=1
    error([mfilename,': bad size']);
end

%% Actual function code
v=repmat(x,nn,1).^repmat((1:nn)',1,n);

if nargout>2
    %% Numerical Jacobian

    % FMT is function handle to repackage vector argument to what
    % the function expects.
    if cX
        fmt=@(x)reshape(x,1,[]);
        fun=@(x)feval(mfilename,fmt(x),nn);
        dvn.dX=jacapprox(fun,x);
    end
end

if nargout>1
    %% Analytical Jacobian
    if cX
        % Every block is [1;2*x;...;n*x^(n-1)]
        
        % Values in memory order
        vv=repmat((1:nn)',1,n).*repmat(x,nn,1).^repmat((0:nn-1)',1,n);
        % Row indices
        ii=1:n*nn;
        % Column indices
        jj=repmat(1:n,nn,1);
        dv.dX=sparse(ii,jj,vv,n*nn,n);
    end
end


function fail=selftest(verbose)

% Set up test data.
m=5;
x=rand(1,m);
n=4;

fail=full_self_test(mfilename,{x,n},1e-8,1e-8,verbose);

