function [Q,dQ,dQn]=similar3(k,P,p0,varargin)
%SIMILAR3 3D similarity transform for the DBAT projection model
%
%   Q=SIMILAR3(K,P,P0) applies a translation and rotation on the 3D
%   points in the 3-by-N array P. The 3-vector P0 holds the camera
%   center. The vector K contains the rotation angles for the rotation
%   matrix M. The actual transformation is M(K)*(P-P0).
%
%   [Q,dQ]=... also returns a struct dQ with the analytical Jacobians
%   with respect to K, P, and P0 in the fields dK, dP, and dP0. For
%   more details, see DBAT_BUNDLE_FUNCTIONS.
%
%SEE ALSO: XLAT3, ROTMATEULER, DBAT_BUNDLE_FUNCTIONS

% Treat selftest call separately.
if nargin>=1 && ischar(K), Q=selftest(nargin>1 && P); return; end

% Otherwise, verify number of parameters.
narginchk(3,6);

Q=[];
dQ=[];
dQn=[];

if nargout>1
    % Construct empty Jacobian struct.
    dQ=struct('dP',[],...
              'dP0',[],...
              'dK',[]);
    dQn=dQ;
end

% What Jacobians to compute?
cK=nargout>1 && (length(varargin)<1 || varargin{1});
cP=nargout>1 && (length(varargin)<2 || varargin{2});
cP0=nargout>1 && (length(varargin)<3 || varargin{3});

%% Test parameters
[m,n]=size(P);
if m~=3 || any(size(p0)~=[3,1]) || any(size(k)~=[3,1])
    error([mfilename,': bad size']);
end

%% Actual function code
if nargout<2
    M=eulerrotmat(k,123,false)';
    Q=lin3(M,xlat3(P,p0));
else
    [M,dM]=eulerrotmat(k,123,false)';
    Q=
end

if nargout>2
    %% Numerical Jacobian

    % FMT is function handle to repackage vector argument to what
    % the function expects.
    if cP
        fmt=@(P)reshape(P,3,[]);
        fun=@(P)feval(mfilename,M,fmt(P));
        dQn.dP=jacapprox(fun,P);
    end
    if cM
        fmt=@(M)reshape(M,3,[]);
        fun=@(M)feval(mfilename,fmt(M),P);
        dQn.dM=jacapprox(fun,M);
    end
end

if nargout>1
    %% Analytical Jacobian
    if cP
        dQ.dP=kron(speye(n),M);
    end
    if cM
        dQ.dM=kron(P',speye(m));
    end    
end


function fail=selftest(verbose)

% Set up test data.
m=3;
n=5;
M=rand(m);
P=rand(m,n);

fail=full_self_test(mfilename,{M,P},1e-8,1e-8,verbose);
