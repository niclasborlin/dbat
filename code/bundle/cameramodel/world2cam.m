function [Q,dQ,dQn]=world2cam(P,p0,M,varargin)
%WORLD2CAM Transform 3D points from the world to the camera coordinate system.
%
%   Q=WORLD2CAM(P,P0,M), where P is 3-by-N with 3D points in world
%   coordinates, M is a 3-by-3 rotation matrix, and P0 is 3-by-1 with
%   the camera center, transforms the points to the camera coordinate
%   system. The actual transformation is M*(P-P0).
%
%   [Q,dQ]=... also returns a struct dQ with the analytical Jacobians
%   with respect to P, P0, and M in the fields dP, dP0, and dM,
%   respectively. For more details, see DBAT_BUNDLE_FUNCTIONS.
%
%SEE ALSO: XLAT3, LIN3, DBAT_BUNDLE_FUNCTIONS

% Treat selftest call separately.
if nargin>=1 && ischar(P), Q=selftest(nargin>1 && p0); return; end

% Otherwise, verify number of parameters.
narginchk(3,6);

Q=[]; %#ok<NASGU>
dQ=[];
dQn=[];

if nargout>1
    % Construct empty Jacobian struct.
    dQ=struct('dP',[],...
              'dP0',[],...
              'dM',[]);
    dQn=dQ;
end

% What Jacobians to compute?
cP=nargout>1 && (length(varargin)<1 || varargin{1});
cP0=nargout>1 && (length(varargin)<2 || varargin{2});
cM=nargout>1 && (length(varargin)<3 || varargin{3});

%% Test parameters
[m,n]=size(P);
if m~=3 || any(size(p0)~=[3,1]) || any(size(M)~=3)
    error([mfilename,': bad size']);
end

%% Actual function code
if nargout<2
    Q=lin3(xlat3(P,-p0),M);
else
    X=xlat3(P,-p0);
    [Q,dL]=lin3(X,M,cP | cP0,cM);
end

if nargout>2
    %% Numerical Jacobian

    % FMT is function handle to repackage vector argument to what
    % the function expects.
    if cM
        fmt=@(M)reshape(M,3,3);
        fun=@(M)feval(mfilename,P,p0,fmt(M));
        dQn.dM=jacapprox(fun,M);
    end
    if cP
        fmt=@(P)reshape(P,3,[]);
        fun=@(P)feval(mfilename,fmt(P),p0,M);
        dQn.dP=jacapprox(fun,P);
    end
    if cP0
        fun=@(p0)feval(mfilename,P,p0,M);
        dQn.dP0=jacapprox(fun,p0);
    end
end

if nargout>1
    %% Analytical Jacobian
    if cM
        dQ.dM=dL.dM;
    end    
    if cP
        dQ.dP=dL.dP;
    end
    if cP0
        dQ.dP0=repmat(-M,n,1);
    end
end


function fail=selftest(verbose)

% Set up test data.
m=3;
n=5;

A=rand(3);
[M,~]=qr(A);
if det(M)<0
    M=-M;
end
p0=rand(m,1);
P=rand(m,n);

fail=full_self_test(mfilename,{P,p0,M},1e-8,1e-8,verbose);
