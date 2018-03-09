function [M,dM,dMn]=eulerrotmat(k,seq,fixed)
%EULERROTMAT 3D Euler rotation matrix.
%
%   M=EULERROTMAT(K,SEQ,FIXED) computes the 3D rotation matrix that
%   correspond to the Euler angles in the 3-vector K for the axis
%   sequence given by the 3-digit integer SEQ. Each digit can be 1
%   (X), 2 (Y), or 3 (Z), and determines the axes about which the
%   elementary rotations takes place. For instance, SEQ=123
%   corresponds to x-y-z rotations, SEQ=313 corresponds to z-x-z
%   rotations. If FIXED is TRUE, the rotations are performed in a
%   fixed coordinate system. If FIXED is FALSE, the rotations are
%   performed in a moving coordinate system.
%
%   [M,dM]=... also returns a struct dM with the analytical Jacobian
%   with respect to K in the field dK. For more details, see
%   DBAT_BUNDLE_FUNCTIONS.
%
%SEE ALSO: DBAT_BUNDLE_FUNCTIONS

% Treat selftest call separately.
if nargin>=1 && ischar(k), M=selftest(nargin>1 && seq); return; end

% Otherwise, verify number of parameters.
narginchk(3,3);

M=[];
dM=[];
dMn=[];

if nargout>1
    % Construct empty Jacobian struct.
    dM=struct('dK',[]);
    dMn=dM;
end

%% Test parameters
if length(k)~=3
    error([mfilename,': bad size']);
end

%% Actual function code
M=P+repmat(c,1,n);

if nargout>2
    %% Numerical Jacobian

    % FMT is function handle to repackage vector argument to what
    % the function expects.
    if cP
        fmt=@(P)reshape(P,2,[]);
        fun=@(P)feval(mfilename,fmt(P),c);
        dMn.dP=jacapprox(fun,P);
    end
    if cC
        fun=@(c)feval(mfilename,P,c);
        dMn.dC=jacapprox(fun,c);
    end
end

if nargout>1
    %% Analytical Jacobian
    if cP
        dM.dP=speye(numel(P));
    end
    if cC
        dM.dC=repmat(speye(m),n,1);
    end    
end

% Elementary rotations about each axis.

function R=R1(alpha)

R=[1,0,0;0,cos(alpha),-sin(alpha);0,sin(alpha),cos(alpha)];


function R=R2(alpha)

R=[cos(alpha),0,sin(alpha);0,1,0;-sin(alpha),0,cos(alpha)];


function R=R3(alpha)

R=[cos(alpha),-sin(alpha),0;sin(alpha),cos(alpha),0;0,0,1];


function fail=selftest(verbose)

% Set up test data.
n=2;
m=5;
c=rand(n,1);
P=rand(n,m);

fail=full_self_test(mfilename,{P,c},1e-8,1e-8,verbose);
