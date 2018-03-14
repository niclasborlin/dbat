function [Q,dQ,dQn]=eulerpinhole(P,ang,p0,f,varargin)
%EULERPINHOLE 3D-to-2D pinhole camera transformation with Euler angles.
%
%   Q=EULERPINHOLE(P,A,P0,F) performs a 3D-to-2D transformation of the
%   3D points in the 3-by-N array P according to the pinhole camera
%   model. The position of the camera center is given by the 3-vector
%   P0. The world-to-camera rotation is specified by the 3-vector A
%   with the omega-phi-kappa Euler angles. The focal length of the
%   camera is given by the scalar F. The returned 2-by-N array with
%   image coordinates in the same unit as F.
%
%   [Q,dQ]=... also returns a struct dQ with the analytical Jacobians
%   with respect to P, A, P0, and F in the fields dP, dA, dP0, and dF,
%   respectively. For more details, see DBAT_BUNDLE_FUNCTIONS.
%
%SEE ALSO: WORLD2CAM, ROTMATEULER, DBAT_BUNDLE_FUNCTIONS

% Treat selftest call separately.
if nargin>=1 && ischar(P), Q=selftest(nargin>1 && ang); return; end

% Otherwise, verify number of parameters.
narginchk(4,8);

Q=[];
dQ=[];
dQn=[];

if nargout>1
    % Construct empty Jacobian struct.
    dQ=struct('dA',[],...
              'dF',[],...
              'dP0',[],...
              'dP',[]);
    dQn=dQ;
end

% What Jacobians to compute?
cP=nargout>1 && (length(varargin)<1 || varargin{1});
cA=nargout>1 && (length(varargin)<2 || varargin{2});
cP0=nargout>1 && (length(varargin)<3 || varargin{3});
cF=nargout>1 && (length(varargin)<4 || varargin{4});

%% Test parameters
[m,n]=size(P);
if m~=3 || any(size(p0)~=[3,1]) || any(size(ang)~=[3,1]) || ~isscalar(f)
    error([mfilename,': bad size']);
end

%% Actual function code
if nargout<2
    Q=f*pinhole(world2cam(P,eulerrotmat(ang,123,false)',p0));
else
    if cA
        [M,dM]=eulerrotmat(ang,123,false);

        % We want M', so apply corresponding permutation to dM.
        i=[1,4,7,2,5,8,3,6,9];
        MT=M';
        dMT=dM;
        dMT.dA=dMT.dA(i,:);
    else
        MT=eulerrotmat(ang,123,false)';
    end
    
    [W2C,dW2C]=world2cam(P,MT,p0,cP,cA,cP0);
    [PH,dPH]=pinhole(W2C,cA | cP | cP0);
    Q=f*PH;
end

if nargout>2
    %% Numerical Jacobian

    % FMT is function handle to repackage vector argument to what
    % the function expects.
    if cF
        fun=@(f)feval(mfilename,P,ang,p0,f);
        dQn.dF=jacapprox(fun,f);
    end
    if cA
        fun=@(ang)feval(mfilename,P,ang,p0,f);
        dQn.dA=jacapprox(fun,ang);
    end
    if cP
        fmt=@(P)reshape(P,3,[]);
        fun=@(P)feval(mfilename,fmt(P),ang,p0,f);
        dQn.dP=jacapprox(fun,P);
    end
    if cP0
        fun=@(p0)feval(mfilename,P,ang,p0,f);
        dQn.dP0=jacapprox(fun,p0);
    end
end

if nargout>1
    %% Analytical Jacobian
    if cF
        dQ.dF=vec(PH);
    end
    if cA
        dQ.dA=f*dPH.dP*dW2C.dM*dMT.dA;
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
ang=rand(3,1)*pi/6;
p0=rand(m,1);
P=rand(m,n)+3;

fail=full_self_test(mfilename,{P,ang,p0,f},1e-8,1e-8,verbose);
