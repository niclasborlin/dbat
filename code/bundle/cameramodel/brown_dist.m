function [v,dv,dvn]=brown_dist(u,K,P,varargin)
%BROWN_DIST Brown radial distortion with absolute coordinates.
%
%   V=BROWN_DIST(U,K,P) applies lens distortion according to Brown
%   (1971) to the 2D points in the 2-by-N array U. The vector K
%   contain the radial coefficients. The vector P contain the
%   tangential coefficients.
%
%   [V,dV]=... also returns a struct dV with the analytical Jacobians
%   with respect to U, K, and P in the field dU, dK, and dP,
%   respectively. For more details, see DBAT_BUNDLE_FUNCTIONS.
%
%   References: Brown (1971), "Close-range camera calibration".
%       Photogrammetric Engineering, 37(8): 855-866.
%
%SEE ALSO: BROWN_TANG, BROWN_RAD, DBAT_BUNDLE_FUNCTIONS.

% Treat selftest call separately.
if nargin>=1 && ischar(u), v=selftest(nargin>1 && K); return; end

% Otherwise, verify number of parameters.
narginchk(1,6);

v=[]; %#ok<NASGU>
dv=[];
dvn=[];

if nargout>1
    % Construct empty Jacobian struct.
    dv=struct('dK',[],...
              'dP',[],...
              'dU',[]);
    dvn=dv;
end

% What Jacobians to compute?
cU=nargout>1 && (length(varargin)<1 || varargin{1});
cK=nargout>1 && (length(varargin)<2 || varargin{2});
cP=nargout>1 && (length(varargin)<3 || varargin{3});

%% Test parameters
[um,un]=size(u);
[~,kn]=size(K);
[~,pn]=size(P);
if um~=2 || (kn~=1 && ~isempty(K)) || (pn~=1 && ~isempty(P))
    error([mfilename,': bad size']);
end

%% Actual function code
if nargout<2
    % Only need the function values.
    v=u+brown_rad(u,K)+brown_tang(u,P);
else
    % Need Jacobians too.
    [br,dbr]=brown_rad(u,K,cU,cK);
    [bt,dbt]=brown_tang(u,P,cU,cP);
    v=u+br+bt;
end

if nargout>2
    %% Numerical Jacobian

    % FMT is function handle to repackage vector argument to what
    % the function expects.
    if cK
        fun=@(K)feval(mfilename,u,K,P);
        dvn.dK=jacapprox(fun,K);
    end
    if cP
        fun=@(P)feval(mfilename,u,K,P);
        dvn.dP=jacapprox(fun,P);
    end
    if cU
        fmt=@(u)reshape(u,2,[]);
        fun=@(u)feval(mfilename,fmt(u),K,P);
        dvn.dU=jacapprox(fun,u);
    end
end

if nargout>1
    %% Analytical Jacobian
    if cK
        dv.dK=dbr.dK;
    end
    if cP
        dv.dP=dbt.dP;
    end
    if cU
        dv.dU=speye(2*un)+dbr.dU+dbt.dU;
    end
end


function fail=selftest(verbose)

m=3;
u=rand(2,m);

fail=false;

% Set up test data.
for kl=0:4
    K=rand(kl,1);
    for pl=[0,2:5]
        P=rand(pl,1);
        fail=fail | full_self_test(mfilename,{u,K,P},1e-8,1e-8,verbose);
    end 
end
