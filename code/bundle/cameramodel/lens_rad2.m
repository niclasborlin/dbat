function [r,dr,drn]=lens_rad2(U,varargin)
%LENS_RAD2 Radius squared function used by DBAT lens distortion functions.
%
%   R=LENS_RAD2(U) computes the radius squared for each column in
%   the 2-by-N array U.
%
%   [R,dR]=... also returns a struct dR with the analytical Jacobian
%   with respect to U in the field dU. For more details, see
%   DBAT_BUNDLE_FUNCTIONS.
%
%SEE ALSO: DBAT_BUNDLE_FUNCTIONS

% Treat selftest call separately.
if nargin>=1 && ischar(U), r=selftest(nargin>1 && varargin{1}); return; end

% Otherwise, verify number of parameters.
narginchk(1,2);

r=[]; %#ok<NASGU>
dr=[];
drn=[];

if nargout>1
    % Construct empty Jacobian struct.
    dr=struct('dU',[]);
    drn=dr;
end

% What Jacobians to compute?
cU=nargout>1 && (length(varargin)<1 || varargin{1});

%% Test parameters
[m,n]=size(U);
if m~=2
    error([mfilename,': bad size']);
end

%% Actual function code
r=sum(U.^2,1)';

if nargout>2
    %% Numerical Jacobian

    % FMT is function handle to repackage vector argument to what
    % the function expects.
    if cU
        fmt=@(U)reshape(U,2,[]);
        fun=@(U)feval(mfilename,fmt(U));
        drn.dU=jacapprox(fun,U);
    end
end

if nargout>1
    %% Analytical Jacobian
    if cU
        % Each block is u'
        
        % Values in memory order.
        vv=2*U;
        % Row indices
        ii=repmat(1:n,m,1);
        % Column indices
        jj=1:m*n;
        
        dr.dU=sparse(ii,jj,vv,n,m*n);
    end
end


function fail=selftest(verbose)

% Set up test data.
n=2;
m=5;
U=rand(n,m);

fail=full_self_test(mfilename,{U},1e-8,1e-8,verbose);

