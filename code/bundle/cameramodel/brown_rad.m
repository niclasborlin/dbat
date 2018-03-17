function [v,dv,dvn]=brown_rad(u,K,varargin)
%BROWN_RAD Brown radial distortion.
%
%   V=BROWN_RAD(U,K) returns the radial distortion of Brown (1971) for
%   each 2D point in the 2-by-N array U. The vector K contain the
%   radial coefficients.
%
%   [V,dV]=... also returns a struct dV with the analytical Jacobians
%   with respect to U and K in the field dU and dK, respectively. For
%   more details, see DBAT_BUNDLE_FUNCTIONS.
%
%   References: Brown (1971), "Close-range camera calibration".
%       Photogrammetric Engineering, 37(8): 855-866.
%
%SEE ALSO: BROWN_TANG, BROWN_DIST_ABS, DBAT_BUNDLE_FUNCTIONS.

% Treat selftest call separately.
if nargin>=1 && ischar(u), v=selftest(nargin>1 && K); return; end

% Otherwise, verify number of parameters.
narginchk(1,4);

v=[]; %#ok<NASGU>
dv=[];
dvn=[];

if nargout>1
    % Construct empty Jacobian struct.
    dv=struct('dK',[],...
              'dU',[]);
    dvn=dv;
end

% What Jacobians to compute?
cU=nargout>1 && (length(varargin)<1 || varargin{1});
cK=nargout>1 && (length(varargin)<2 || varargin{2});

%% Test parameters
[um,un]=size(u);
[km,kn]=size(K);
if um~=2 || kn~=1
    error([mfilename,': bad size']);
end

%% Actual function code
if nargout<2
    % Only need the function values.
    v=u.*repmat(rad_scale(u,K)',2,1);
else
    % Need Jacobians too.
    [rs,drs]=rad_scale(u,K,cU,cK);
    v=u.*repmat(rs',2,1);
end

if nargout>2
    %% Numerical Jacobian

    % FMT is function handle to repackage vector argument to what
    % the function expects.
    if cK
        fun=@(K)feval(mfilename,u,K);
        dvn.dK=jacapprox(fun,K);
    end
    if cU
        fmt=@(u)reshape(u,2,[]);
        fun=@(u)feval(mfilename,fmt(u),K);
        dvn.dU=jacapprox(fun,u);
    end
end

if nargout>1
    %% Analytical Jacobian
    if cK
        % Each row block is u*drs.dC
        
        % Expand u and K
        U=repmat(u(:),1,km);
        ii=repmat(1:un,2,1);
        dK=drs.dC(ii(:),:);
        dv.dK=U.*dK;
    end
    if cU
        % Each diagonal block is rs*eye(2)+u*drs.dU
        
        % Build values in memory order.
        vv=u*drs.dU;
        vv(1,1:2:end)=vv(1,1:2:end)+rs';
        vv(2,2:2:end)=vv(2,2:2:end)+rs';
        % Row indices.
        ii=repmat([1,2,1,2]',1,un)+repmat(2*(0:un-1),4,1);
        % Column indices.
        jj=repmat([1,1,2,2]',1,un)+repmat(2*(0:un-1),4,1);
        dv.dU=sparse(ii,jj,vv,2*un,2*un);
    end
end


function fail=selftest(verbose)

% Set up test data.
K=rand(4,1);
m=3;
u=rand(2,m);

fail=full_self_test(mfilename,{u,K},1e-8,1e-8,verbose);

