function [Q,dQ,dQn]=affine2(b1,b2,P,varargin)
%AFFINE2 2D affine transform for the DBAT projection model
%
%   Q=AFFINE2(B1,B2,P) applies the linear transform AFFINE2MAT(B1,B2)
%   to the 2D points in the 2-by-N array P, i.e. computes
%   AFFINE2MAT(B1,B2)*P.
%
%   [Q,dQ]=... also returns a struct dQ with the analytical Jacobians
%   with respect to B1, B2, and P in the fields dB1, dB2, and dP. For
%   more details, see DBAT_BUNDLE_FUNCTIONS.
%
%SEE ALSO: DBAT_BUNDLE_FUNCTIONS

% Treat selftest call separately.
if nargin>=1 && ischar(b1), Q=selftest(nargin>1 && b2); return; end

% Otherwise, verify number of parameters.
narginchk(3,6);

Q=[];
dQ=[];
dQn=[];

if nargout>1
    % Construct empty Jacobian struct.
    dQ=struct('dP',[],...
              'dB1',[],...
              'dB2',[]);
    dQn=dQ;
end

% What Jacobians to compute?
cB1=nargout>1 && (length(varargin)<1 || varargin{1});
cB2=nargout>1 && (length(varargin)<2 || varargin{2});
cP=nargout>1 && (length(varargin)<3 || varargin{3});

%% Test parameters
[m,n]=size(P);
if m~=2 || ~isscalar(b1) || ~isscalar(b2)
    error([mfilename,': bad size']);
end

%% Actual function code
A=affine2mat(b1,b2);
Q=A*P;

if nargout>2
    %% Numerical Jacobian

    % FMT is function handle to repackage vector argument to what
    % the function expects.
    if cP
        fmt=@(P)reshape(P,2,[]);
        fun=@(P)feval(mfilename,b1,b2,fmt(P));
        dQn.dP=jacapprox(fun,P);
    end
    if cB1
        fun=@(b1)feval(mfilename,b1,b2,P);
        dQn.dB1=jacapprox(fun,b1);
    end
    if cB2
        fun=@(b2)feval(mfilename,b1,b2,P);
        dQn.dB2=jacapprox(fun,b2);
    end
end

if nargout>1
    %% Analytical Jacobian
    if cP
        % Each block is A = [1+B1, B2; 0 1].

        % Values for Jacobian in memory-order.
        vv=repmat(A([1,3,4])',1,n);
        % Row indices.
        ii=repmat(0:2:2*n-1,3,1)+repmat([1,1,2]',1,n);
        % Column indices.
        jj=repmat(0:2:2*n-1,3,1)+repmat([1,2,2]',1,n);
        dQ.dP=sparse(ii,jj,vv,2*n,2*n);
    end
    if cB1
        dQ.dB1=zeros(2*n,1);
        dQ.dB1(1:2:end)=P(1,:)';
    end    
    if cB2
        dQ.dB2=zeros(2*n,1);
        dQ.dB2(1:2:end)=P(2,:)';
    end    
end


function fail=selftest(verbose)

% Set up test data.
m=2;
n=5;
b1=rand+0.1;
b2=rand+0.1;
P=rand(m,n);

fail=full_self_test(mfilename,{b1,b2,P},1e-8,1e-8,verbose);
