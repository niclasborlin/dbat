function [Q,dQ,dQn]=pinhole(P,varargin)
%PINHOLE Pinhole projection for the DBAT projection model.
%
%   Q=PINHOLE(P) applies the pinhole projection on the 3D points in
%   the 3-by-N array P, i.e. computes P(1:2)/P(3) for each column.
%
%   [Q,dQ]=... also returns a struct dQ with the analytical Jacobian
%   with respect to P in the field dP. For more details, see
%   DBAT_BUNDLE_FUNCTIONS.
%
%SEE ALSO: DBAT_BUNDLE_FUNCTIONS

% Treat selftest call separately.
if nargin>=1 && ischar(P), Q=selftest(nargin>1 && varargin{1}); return; end

% Otherwise, verify number of parameters.
narginchk(1,2);

Q=[]; %#ok<NASGU>
dQ=[];
dQn=[];

if nargout>1
    % Construct empty Jacobian struct.
    dQ=struct('dP',[]);
    dQn=dQ;
end

% What Jacobians to compute?
cP=nargout>1 && (length(varargin)<1 || varargin{1});

%% Test parameters
[m,n]=size(P);
if m~=3
    error([mfilename,': bad size']);
end

%% Actual function code
Q=P(1:2,:)./P([3,3],:);

if nargout>2
    %% Numerical Jacobian

    % FMT is function handle to repackage vector argument to what
    % the function expects.
    if cP
        fmt=@(P)reshape(P,3,[]);
        fun=@(P)feval(mfilename,fmt(P));
        dQn.dP=jacapprox(fun,P);
    end
end

if nargout>1
    %% Analytical Jacobian
    if cP
        % Each block is 1/P(3)*[1, 0, P(1)/P(3)
        %                       0, 1, P(2)/P(3)]
        P12=P(1:2,:);
        P33i=1./P([3,3],:);
        % Values for Jacobian in memory-order.
        vv=[P33i;-P12.*P33i.^2];
        % Row indices.
        ii=repmat(0:2:2*n-1,4,1)+repmat([1,2,1,2]',1,n);
        % Column indices.
        jj=repmat(0:3:3*n-1,4,1)+repmat([1,2,3,3]',1,n);
        dQ.dP=sparse(ii,jj,vv,2*n,3*n);
    end
end


function fail=selftest(verbose)

% Set up test data.
m=3;
n=5;
P=rand(m,n)+1;

fail=full_self_test(mfilename,{P},1e-8,1e-8,verbose);
