function [v,dv,dvn]=brown_dist_rel(u,u0,K,P,varargin)
%BROWN_DIST_REL Brown radial distortion with relative coordinates.
%
%   V=BROWN_DIST_REL(U,U0,K,P) returns the lens distortion of Brown
%   (1971) for each 2D point in the 2-by-N array U. The vector K
%   contain the radial coefficients. The vector P contain the
%   tangential coefficients. The distortion should be computed with
%   respect to the 2-by-1 principal point U0.
%
%   [V,dV]=... also returns a struct dV with the analytical Jacobians
%   with respect to U, U0, K, and P in the field dU, dU0, dK, and dP,
%   respectively. For more details, see DBAT_BUNDLE_FUNCTIONS.
%
%   References: Brown (1971), "Close-range camera calibration".
%       Photogrammetric Engineering, 37(8): 855-866.
%
%SEE ALSO: BROWN_DIST_ABS, DBAT_BUNDLE_FUNCTIONS.

% Treat selftest call separately.
if nargin>=1 && ischar(u), v=selftest(nargin>1 && u0); return; end

% Otherwise, verify number of parameters.
narginchk(1,8);

v=[];
dv=[];
dvn=[];

if nargout>1
    % Construct empty Jacobian struct.
    dv=struct('dK',[],...
              'dP',[],...
              'dU0',[],...
              'dU',[]);
    dvn=dv;
end

% What Jacobians to compute?
cU=nargout>1 && (length(varargin)<1 || varargin{1});
cU0=nargout>1 && (length(varargin)<2 || varargin{2});
cK=nargout>1 && (length(varargin)<3 || varargin{3});
cP=nargout>1 && (length(varargin)<4 || varargin{4});

%% Test parameters
[um,un]=size(u);
[km,kn]=size(K);
[pm,pn]=size(P);
if um~=2 || (kn~=1 && ~isempty(K)) || (pn~=1 && ~isempty(P)) || ...
       any(size(u0)~=[2,1])
    error([mfilename,': bad size']);
end

%% Actual function code
if nargout<2
    % Only need the function values.
    v=xlat2(brown_dist_abs(xlat2(u,-u0),K,P),u0);
else
    % Need Jacobians too.
    [bda,dbda]=brown_dist_abs(xlat2(u,-u0),K,P,cU | cU0,cK,cP);
    v=xlat2(bda,u0);
end

if nargout>2
    %% Numerical Jacobian

    % FMT is function handle to repackage vector argument to what
    % the function expects.
    if cK
        fun=@(K)feval(mfilename,u,u0,K,P);
        dvn.dK=jacapprox(fun,K);
    end
    if cP
        fun=@(P)feval(mfilename,u,u0,K,P);
        dvn.dP=jacapprox(fun,P);
    end
    if cU0
        fun=@(u0)feval(mfilename,u,u0,K,P);
        dvn.dU0=jacapprox(fun,u0);
    end
    if cU
        fmt=@(u)reshape(u,2,[]);
        fun=@(u)feval(mfilename,fmt(u),u0,K,P);
        dvn.dU=jacapprox(fun,u);
    end
end

if nargout>1
    %% Analytical Jacobian
    if cK
        dv.dK=dbda.dK;
    end
    if cP
        dv.dP=dbda.dP;
    end
    if cU
        dv.dU=dbda.dU;
    end
    if cU0
        [ii,jj,vv]=find(dbda.dU);
        % Set to 1-vv on the main diagonal.
        vv=(ii==jj)-vv;
        % In most cases, the dU matrix will have full diagonal blocks.
        if length(vv)==un*4
            dU0=[vv(rem(jj,2)==1),vv(rem(jj,2)==0)];
        else
            dU0=zeros(2*un,2);
            if nnz(vv)>0
                % All-zero vv happens e.g. when K and P are empty.
                kk=sub2ind(size(dU0),ii,rem(jj-1,2)+1);
                dU0(kk)=vv;
            end
        end
        dv.dU0=dU0;
    end
end


function fail=selftest(verbose)

m=3;
u=rand(2,m);
u0=rand(2,1);

fail=false;

% Set up test data.
for kl=0:4
    K=rand(kl,1);
    for pl=[0,2:5]
        P=rand(pl,1);
        fail=fail | full_self_test(mfilename,{u,u0,K,P},1e-8,1e-8,verbose);
    end 
end
