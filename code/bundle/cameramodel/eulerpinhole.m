function [Q,dQ,dQn]=eulerpinhole(f,k,P,p0,varargin)
%EULERPINHOLE 3D-to-2D pinhole camera transformation with Euler angles.
%
%   Q=EULERPINHOLE(F,K,P,P0) performs a 3D-to-2D transformation of the
%   3D points in the 3-by-N array P according to the pinhole camera
%   model. The scalar F contains the focal length. The 3-vector P0
%   contains the camera center. The 3-vector K contain the
%   omega-phi-kappa Euler angles. The returned 2-by-N array with
%   image coordinates in the same unit as F.
%
%   [Q,dQ]=... also returns a struct dQ with the analytical Jacobians
%   with respect to F, K, P, and P0 in the fields dF, dK, dP, and dP0.
%   For more details, see DBAT_BUNDLE_FUNCTIONS.
%
%SEE ALSO: WORLD2CAM, ROTMATEULER, DBAT_BUNDLE_FUNCTIONS

% Treat selftest call separately.
if nargin>=1 && ischar(f), Q=selftest(nargin>1 && k); return; end

% Otherwise, verify number of parameters.
narginchk(4,8);

Q=[];
dQ=[];
dQn=[];

if nargout>1
    % Construct empty Jacobian struct.
    dQ=struct('dK',[],...
              'dF',[],...
              'dP0',[],...
              'dP',[]);
    dQn=dQ;
end

% What Jacobians to compute?
cF=nargout>1 && (length(varargin)<1 || varargin{1});
cK=nargout>1 && (length(varargin)<2 || varargin{2});
cP=nargout>1 && (length(varargin)<3 || varargin{3});
cP0=nargout>1 && (length(varargin)<4 || varargin{4});

%% Test parameters
[m,n]=size(P);
if m~=3 || any(size(p0)~=[3,1]) || any(size(k)~=[3,1]) || ~isscalar(f)
    error([mfilename,': bad size']);
end

%% Actual function code
if nargout<2
    Q=f*pinhole(world2cam(eulerrotmat(k,123,false)',P,p0));
else
    if cK
        [M,dM]=eulerrotmat(k,123,false);

        % We want M', so apply corresponding permutation to dM.
        i=[1,4,7,2,5,8,3,6,9];
        MT=M';
        dMT=dM;
        dMT.dK=dMT.dK(i,:);
    else
        MT=eulerrotmat(k,123,false)';
    end
    
    [W2C,dW2C]=world2cam(MT,P,p0,cK,cP,cP0);
    [PH,dPH]=pinhole(W2C,cK | cP | cP0);
    Q=f*PH;
end

if nargout>2
    %% Numerical Jacobian

    % FMT is function handle to repackage vector argument to what
    % the function expects.
    if cF
        fun=@(f)feval(mfilename,f,k,P,p0);
        dQn.dF=jacapprox(fun,f);
    end
    if cK
        fun=@(k)feval(mfilename,f,k,P,p0);
        dQn.dK=jacapprox(fun,k);
    end
    if cP
        fmt=@(P)reshape(P,3,[]);
        fun=@(P)feval(mfilename,f,k,fmt(P),p0);
        dQn.dP=jacapprox(fun,P);
    end
    if cP0
        fun=@(p0)feval(mfilename,f,k,P,p0);
        dQn.dP0=jacapprox(fun,p0);
    end
end

if nargout>1
    %% Analytical Jacobian
    if cF
        dQ.dF=vec(PH);
    end
    if cK
        dQ.dK=f*dPH.dP*dW2C.dM*dMT.dK;
    end    
    if cP
        dQ.dP=f*dPH.dP*dW2C.dP;
    end
    if cP0
        dQ.dP0=f*dPH.dP*dW2C.dP0;
    end
end


function fail=selftest(verbose)

% Set up test data.
m=3;
n=5;

f=1+rand;
k=rand(3,1)*pi/6;
p0=rand(m,1);
P=rand(m,n)+3;

fail=full_self_test(mfilename,{f,k,P,p0},1e-8,1e-8,verbose);
