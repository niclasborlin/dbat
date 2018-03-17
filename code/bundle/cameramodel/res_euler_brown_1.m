function [v,dv,dvn]=res_euler_brown_1(Q,q0,ang,f,u,sz,u0,K,P,b,varargin)
%RES_EULER_BROWN_1 Residual function 1 for DBAT.
%
%   V=RES_EULER_BROWN_1(Q,Q0,A,F,U,SZ,U0,K,P,B) computes the residual
%   between the projected object points in the 3-by-N array Q and the
%   measured image coordinates in the 2-by-N array U. The object
%   points are projected using the pinhole projection model with the
%   camera center in Q0 (3-by-1). The rotation from the world to the
%   camera is given by the Euler x-y-z angles in A (3-by-1). The
%   camera has a focal length F (scalar), principal point U0 (2-by-1),
%   and pixel size SZ (scalar). The residual is formed between the
%   projected object points and the measured image points after the
%   image points have been converted to mm and corrected for lens
%   distortion. The vectors K and P contain the radial and tangential
%   distortion coefficients, respectively. Furthermore, an affine
%   transformation with parameters B (2-by-1) is applied before
%   correction for lens distortion.
%
%   [V,dV]=... also returns a struct dV with the analytical Jacobians
%   with respect to Q, Q0, A, F, U, SZ, U0, K, P, and B in the fields
%   dQ, dQ0, dA, dF, dU, dSZ, dU0, dK, dP, and dB, respectively. For
%   more details, see DBAT_BUNDLE_FUNCTIONS.
%
%   References: Brown (1971), "Close-range camera calibration".
%       Photogrammetric Engineering, 37(8): 855-866.
%
%SEE ALSO: BROWN_DIST, EULERPINHOLE2, AFFINE2, DBAT_BUNDLE_FUNCTIONS.

% Treat selftest call separately.
if nargin>=1 && ischar(Q), v=selftest(nargin>1 && q0); return; end

% Otherwise, verify number of parameters.
narginchk(10,20);

v=[]; %#ok<NASGU>
dv=[];
dvn=[];

if nargout>1
    % Construct empty Jacobian struct.
    dv=struct('dQ',[],...
              'dA',[],...
              'dQ0',[],...
              'dF',[],...
              'dU',[],...
              'dSZ',[],...
              'dU0',[],...
              'dK',[],...
              'dP',[],...
              'dB',[]);
    dvn=dv;
end

% What Jacobians to compute?
cQ =nargout>1 && (length(varargin)<1 || varargin{1});
cQ0=nargout>1 && (length(varargin)<2 || varargin{2});
cA =nargout>1 && (length(varargin)<3 || varargin{3});
cF =nargout>1 && (length(varargin)<4 || varargin{4});
cU =nargout>1 && (length(varargin)<5 || varargin{5});
cSZ=nargout>1 && (length(varargin)<6 || varargin{6});
cU0=nargout>1 && (length(varargin)<7 || varargin{7});
cK =nargout>1 && (length(varargin)<8 || varargin{8});
cP =nargout>1 && (length(varargin)<9 || varargin{9});
cB =nargout>1 && (length(varargin)<10 || varargin{10});

%% Test parameters
[qm,qn]=size(Q);
[um,un]=size(u);
[~,kn]=size(K);
[~,pn]=size(P);
if qm~=3 || any(size(ang)~=[3,1]) || any(size(q0)~=[3,1]) || ~ ...
       isscalar(f) || um~=2 || ~isscalar(sz) || any(size(u0)~=[2,1]) || ...
       (kn~=1 && ~isempty(K)) || (pn~=1 && ~isempty(P)) || any(size(b)~=[2,1])
    error([mfilename,': bad size']);
end

if qn~=un
    error([mfilename,': Q and U should have the same number of columns']);
end

%% Actual function code
if nargout<2
    % Only need the function values.
    lhs=eulerpinhole2(Q,q0,ang,-f);
    rhs=brown_dist(affine2(xlat2(aniscale2(scale2(u,sz),[1;-1]),-u0),b),-K,-P);
    v=lhs-rhs;
else
    % Need Jacobians too.
    [lhs,dlhs]=eulerpinhole2(Q,q0,ang,-f,cQ,cQ0,cA,cF);
    [s,dS]=scale2(u,sz,cU,cSZ);
    [as,dAS]=aniscale2(s,[1;-1],cU | cSZ,false);
    [x,dX]=xlat2(as,-u0,cU | cSZ,cU0);
    [a,dA]=affine2(x,b,cU | cSZ | cU0,cB);
    [l,dL]=brown_dist(a,-K,-P,cU | cSZ | cU0 | cB,cK,cP);
    v=lhs-l;
end

if nargout>2
    %% Numerical Jacobian

    % FMT is function handle to repackage vector argument to what
    % the function expects.
    if cQ
        fmt=@(Q)reshape(Q,3,[]);
        fun=@(Q)feval(mfilename,fmt(Q),q0,ang,f,u,sz,u0,K,P,b);
        dvn.dQ=jacapprox(fun,Q);
    end
    if cA
        fun=@(ang)feval(mfilename,Q,q0,ang,f,u,sz,u0,K,P,b);
        dvn.dA=jacapprox(fun,ang);
    end
    if cQ0
        fun=@(q0)feval(mfilename,Q,q0,ang,f,u,sz,u0,K,P,b);
        dvn.dQ0=jacapprox(fun,q0);
    end
    if cF
        fun=@(f)feval(mfilename,Q,q0,ang,f,u,sz,u0,K,P,b);
        dvn.dF=jacapprox(fun,f);
    end
    if cU
        fmt=@(u)reshape(u,2,[]);
        fun=@(u)feval(mfilename,Q,q0,ang,f,fmt(u),sz,u0,K,P,b);
        dvn.dU=jacapprox(fun,u);
    end
    if cSZ
        fun=@(sz)feval(mfilename,Q,q0,ang,f,u,sz,u0,K,P,b);
        dvn.dSZ=jacapprox(fun,sz);
    end
    if cU0
        fun=@(u0)feval(mfilename,Q,q0,ang,f,u,sz,u0,K,P,b);
        dvn.dU0=jacapprox(fun,u0);
    end
    if cK
        fun=@(K)feval(mfilename,Q,q0,ang,f,u,sz,u0,K,P,b);
        dvn.dK=jacapprox(fun,K);
    end
    if cP
        fun=@(P)feval(mfilename,Q,q0,ang,f,u,sz,u0,K,P,b);
        dvn.dP=jacapprox(fun,P);
    end
    if cB
        fun=@(b)feval(mfilename,Q,q0,ang,f,u,sz,u0,K,P,b);
        dvn.dB=jacapprox(fun,b);
    end
end

if nargout>1
    %% Analytical Jacobian
    if cQ
        dv.dQ=dlhs.dP;
    end
    if cA
        dv.dA=dlhs.dA;
    end
    if cQ0
        dv.dQ0=dlhs.dP0;
    end
    if cF
        dv.dF=-dlhs.dF;
    end
    if cU
        dv.dU=-dL.dU*dA.dU*dX.dU*dAS.dU*dS.dU;
    end
    if cSZ
        dv.dSZ=-dL.dU*dA.dU*dX.dU*dAS.dU*dS.dK;
    end
    if cU0
        dv.dU0=dL.dU*dA.dU*dX.dC;
    end
    if cK
        dv.dK=dL.dK;
    end
    if cP
        dv.dP=dL.dP;
    end
    if cB
        dv.dB=-dL.dU*dA.dB;
    end
end


function fail=selftest(verbose)

m=5;
Q=3+rand(3,m);
ang=rand(3,1)*pi/6;
q0=rand(3,1);
f=1+rand;
u=rand(2,m);
K=rand(4,1);
P=rand(3,1);
sz=rand/10;
u0=rand(2,1);
b=rand(2,1);

fail=full_self_test(mfilename,{Q,q0,ang,f,u,sz,u0,K,P,b},1e-8,1e-8,verbose);
