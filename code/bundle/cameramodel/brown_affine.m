function [v,dv,dvn]=brown_affine(u,sz,u0,K,P,b1,b2,cs,varargin)
%BROWN_AFFINE Brown radial distortion with absolute coordinates.
%
%   V=BROWN_AFFINE(U,SZ,U0,K,P,B1,B2,CS) applies an affine
%   transformation and lens distortion to 2D measured image points in
%   the 2-by-N array U. The pixel size in mm is given by the scalar
%   SZ. The principal point in mm is the 2-by-1 vector U0. The Brown
%   radial and tangential lens distortion coefficients are given by
%   the vectors K and P, respectively. The affine transformation is
%   specified by the scalars B1 and B2. The integer CS determines
%   where the affine transformation is introduced:
%       0 - no affine transformation
%           - scale by sz
%           - subtract u0
%           - apply lens distortion, K, P
%       1 - before lens distortion:
%           - scale by sz
%           - subtract u0
%           - apply affine b1, b2
%           - apply lens distortion, K, P
%       2 - after lens distortion
%           - scale by sz
%           - subtract u0
%           - apply lens distortion, K, P
%           - apply affine b1, b2
%       3 - anisotropic scaling before lens distortion, skew after
%           lens distortion
%           - scale by sz
%           - apply anisotropic scaling by b1
%           - subtract u0
%           - apply lens distortion, K, P
%           - apply skew by b2
%
%   No affine transformation correponds to b1=b2=0 and any value of cs.
%
%   [V,dV]=... also returns a struct dV with the analytical Jacobians
%   with respect to U, SZ, U0, K, P, B1, and B2 in the fields dU, dSZ,
%   dU0, dK, dP, dB1, and dB2, respectively. For more details, see
%   DBAT_BUNDLE_FUNCTIONS.
%
%   References: Brown (1971), "Close-range camera calibration".
%       Photogrammetric Engineering, 37(8): 855-866.
%
%SEE ALSO: BROWN_DIST, AFFINE2, ANISCALE2B, SKEW, DBAT_BUNDLE_FUNCTIONS.

% Treat selftest call separately.
if nargin>=1 && ischar(u), v=selftest(nargin>1 && sz); return; end

% Otherwise, verify number of parameters.
narginchk(8,15);

v=[]; %#ok<NASGU>
dv=[];
dvn=[];

if nargout>1
    % Construct empty Jacobian struct.
    dv=struct('dU',[],...
              'dSZ',[],...
              'dU0',[],...
              'dK',[],...
              'dP',[],...
              'dB1',[],...
              'dB2',[]);
    dvn=dv;
end

% What Jacobians to compute?
cU =nargout>1 && (length(varargin)<1 || varargin{1});
cSZ=nargout>1 && (length(varargin)<2 || varargin{2});
cU0=nargout>1 && (length(varargin)<3 || varargin{3});
cK =nargout>1 && (length(varargin)<4 || varargin{4});
cP =nargout>1 && (length(varargin)<5 || varargin{5});
cB1=nargout>1 && (length(varargin)<6 || varargin{6});
cB2=nargout>1 && (length(varargin)<7 || varargin{7});

%% Test parameters
[um,~]=size(u);
[~,kn]=size(K);
[~,pn]=size(P);
if um~=2 || ~isscalar(sz) || any(size(u0)~=[2,1]) || ...
       (kn~=1 && ~isempty(K)) || (pn~=1 && ~isempty(P)) || ~isscalar(b1) ...
       || ~isscalar(b2) || ~isscalar(cs)
    error([mfilename,': bad size']);
end

%% Actual function code
if nargout<2
    % Only need the function values.
    switch cs
      case 0
        v=brown_dist(xlat2(scale2(u,sz),-u0),K,P);
      case 1
        v=brown_dist(affine2(xlat2(scale2(u,sz),-u0),b1,b2),K,P);
      case 2
        v=affine2(brown_dist(xlat2(scale2(u,sz),-u0),K,P),b1,b2);
      case 3
        v=skew(brown_dist(xlat2(aniscale2b(scale2(u,sz),b1),-u0),K,P),b2);
      otherwise
        error([mfilename,': bad case']);
    end        
else
    % Need Jacobians too.
    switch cs
      case 0
        cB1=false;
        cB2=false;
        [s,dS]=scale2(u,sz,cU,cSZ);
        [x,dX]=xlat2(s,-u0,cU | cSZ,cU0);
        [l,dL]=brown_dist(x,K,P,cU | cSZ | cU0,cK,cP);
        v=l;
      case 1
        [s,dS]=scale2(u,sz,cU,cSZ);
        [x,dX]=xlat2(s,-u0,cU | cSZ,cU0);
        [a,dA]=affine2(x,b1,b2,cU | cSZ | cU0,cB1,cB2);
        [l,dL]=brown_dist(a,K,P,cU | cSZ | cU0 | cB1 | cB2,cK,cP);
        v=l;
      case 2
        [s,dS]=scale2(u,sz,cU,cSZ);
        [x,dX]=xlat2(s,-u0,cU | cSZ,cU0);
        [l,dL]=brown_dist(x,K,P,cU | cSZ | cU0,cK,cP);
        [a,dA]=affine2(l,b1,b2,cU | cSZ | cU0 | cK | cP,cB1,cB2);
        v=a;
      case 3
        [s,dS]=scale2(u,sz,cU,cSZ);
        [as,dAS]=aniscale2b(s,b1,cU | cSZ,cB1);
        [x,dX]=xlat2(as,-u0,cU | cSZ | cB1,cU0);
        [l,dL]=brown_dist(x,K,P,cU | cSZ | cB1 | cU0,cK,cP);
        [sk,dSK]=skew(l,b2,cU | cSZ | cB1 | cU0 | cK | cP,cB2);
        v=sk;
      otherwise
        error([mfilename,': bad case']);
    end        
end

if nargout>2
    %% Numerical Jacobian

    % FMT is function handle to repackage vector argument to what
    % the function expects.
    if cU
        fmt=@(u)reshape(u,2,[]);
        fun=@(u)feval(mfilename,fmt(u),sz,u0,K,P,b1,b2,cs);
        dvn.dU=jacapprox(fun,u);
    end
    if cSZ
        fun=@(sz)feval(mfilename,u,sz,u0,K,P,b1,b2,cs);
        dvn.dSZ=jacapprox(fun,sz);
    end
    if cU0
        fun=@(u0)feval(mfilename,u,sz,u0,K,P,b1,b2,cs);
        dvn.dU0=jacapprox(fun,u0);
    end
    if cK
        fun=@(K)feval(mfilename,u,sz,u0,K,P,b1,b2,cs);
        dvn.dK=jacapprox(fun,K);
    end
    if cP
        fun=@(P)feval(mfilename,u,sz,u0,K,P,b1,b2,cs);
        dvn.dP=jacapprox(fun,P);
    end
    if cB1
        fun=@(b1)feval(mfilename,u,sz,u0,K,P,b1,b2,cs);
        dvn.dB1=jacapprox(fun,b1);
    end
    if cB2
        fun=@(b2)feval(mfilename,u,sz,u0,K,P,b1,b2,cs);
        dvn.dB2=jacapprox(fun,b2);
    end
end

if nargout>1
    %% Analytical Jacobian
    switch cs
      case 0
        if cU
            dv.dU=dL.dU*dX.dU*dS.dU;
        end
        if cSZ
            dv.dSZ=dL.dU*dX.dU*dS.dK;
        end
        if cU0
            dv.dU0=-dL.dU*dX.dC;
        end
        if cK
            dv.dK=dL.dK;
        end
        if cP
            dv.dP=dL.dP;
        end
      case 1
        if cU
            dv.dU=dL.dU*dA.dU*dX.dU*dS.dU;
        end
        if cSZ
            dv.dSZ=dL.dU*dA.dU*dX.dU*dS.dK;
        end
        if cU0
            dv.dU0=-dL.dU*dA.dU*dX.dC;
        end
        if cB1
            dv.dB1=dL.dU*dA.dB1;
        end
        if cB2
            dv.dB2=dL.dU*dA.dB2;
        end
        if cK
            dv.dK=dL.dK;
        end
        if cP
            dv.dP=dL.dP;
        end
      case 2
        if cU
            dv.dU=dA.dU*dL.dU*dX.dU*dS.dU;
        end
        if cSZ
            dv.dSZ=dA.dU*dL.dU*dX.dU*dS.dK;
        end
        if cU0
            dv.dU0=-dA.dU*dL.dU*dX.dC;
        end
        if cK
            dv.dK=dA.dU*dL.dK;
        end
        if cP
            dv.dP=dA.dU*dL.dP;
        end
        if cB1
            dv.dB1=dA.dB1;
        end
        if cB2
            dv.dB2=dA.dB2;
        end
      case 3
        if cU
            dv.dU=dSK.dU*dL.dU*dX.dU*dAS.dU*dS.dU;
        end
        if cSZ
            dv.dSZ=dSK.dU*dL.dU*dX.dU*dAS.dU*dS.dK;
        end
        if cB1
            dv.dB1=dSK.dU*dL.dU*dX.dU*dAS.dK;
        end
        if cU0
            dv.dU0=-dSK.dU*dL.dU*dX.dC;
        end
        if cK
            dv.dK=dSK.dU*dL.dK;
        end
        if cP
            dv.dP=dSK.dU*dL.dP;
        end
        if cB2
            dv.dB2=dSK.dK;
        end
      otherwise
        error([mfilename,': bad case']);
    end        
end


function fail=selftest(verbose)

m=5;
u=rand(2,m);
K=rand(4,1);
P=rand(3,1);
b1=rand;
b2=rand;
sz=rand/10;
u0=rand(2,1);

fail=false;

fail=fail | full_self_test(mfilename,{u,sz,u0,K,P,0,0,0},1e-8,1e-8,verbose,5);

for cs=1:3
    fail=fail | full_self_test(mfilename,{u,sz,u0,K,P,b1,b2,cs},1e-8,1e-8,verbose,7);
end
