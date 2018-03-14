function [Q,dQ,dQn]=affine2(U,b1,b2,varargin)
%AFFINE2 2D affine transform for the DBAT projection model.
%
%   Q=AFFINE2(U,B1,B2) applies the linear transform AFFINE2MAT(B1,B2)
%   to the 2D points in the 2-by-N array U, i.e. computes
%   AFFINE2MAT(B1,B2)*U.
%
%   [Q,dQ]=... also returns a struct dQ with the analytical Jacobians
%   with respect to U, B1, and B2 in the fields dU, dB1, and dB2,
%   respectively. For more details, see DBAT_BUNDLE_FUNCTIONS.
%
%SEE ALSO: DBAT_BUNDLE_FUNCTIONS

% Treat selftest call separately.
if nargin>=1 && ischar(U), Q=selftest(nargin>1 && b1); return; end

% Otherwise, verify number of parameters.
narginchk(3,6);

Q=[];
dQ=[];
dQn=[];

if nargout>1
    % Construct empty Jacobian struct.
    dQ=struct('dU',[],...
              'dB1',[],...
              'dB2',[]);
    dQn=dQ;
end

% What Jacobians to compute?
cU=nargout>1 && (length(varargin)<1 || varargin{1});
cB1=nargout>1 && (length(varargin)<2 || varargin{2});
cB2=nargout>1 && (length(varargin)<3 || varargin{3});

%% Test parameters
[m,n]=size(U);
if m~=2 || ~isscalar(b1) || ~isscalar(b2)
    error([mfilename,': bad size']);
end

%% Actual function code
A=affine2mat(b1,b2);
Q=A*U;

if nargout>2
    %% Numerical Jacobian

    % FMT is function handle to repackage vector argument to what
    % the function expects.
    if cU
        fmt=@(U)reshape(U,2,[]);
        fun=@(U)feval(mfilename,fmt(U),b1,b2);
        dQn.dU=jacapprox(fun,U);
    end
    if cB1
        fun=@(b1)feval(mfilename,U,b1,b2);
        dQn.dB1=jacapprox(fun,b1);
    end
    if cB2
        fun=@(b2)feval(mfilename,U,b1,b2);
        dQn.dB2=jacapprox(fun,b2);
    end
end

if nargout>1
    %% Analytical Jacobian
    if cU
        % Each block is A = [1+B1, B2; 0 1].

        % Values for Jacobian in memory-order.
        vv=repmat(A([1,3,4])',1,n);
        % Row indices.
        ii=repmat(0:2:2*n-1,3,1)+repmat([1,1,2]',1,n);
        % Column indices.
        jj=repmat(0:2:2*n-1,3,1)+repmat([1,2,2]',1,n);
        dQ.dU=sparse(ii,jj,vv,2*n,2*n);
    end
    if cB1
        dQ.dB1=zeros(2*n,1);
        dQ.dB1(1:2:end)=U(1,:)';
    end    
    if cB2
        dQ.dB2=zeros(2*n,1);
        dQ.dB2(1:2:end)=U(2,:)';
    end    
end


function fail=selftest(verbose)

% Set up test data.
m=2;
n=5;
b1=rand+0.1;
b2=rand+0.1;
U=rand(m,n);

fail=full_self_test(mfilename,{U,b1,b2},1e-8,1e-8,verbose);
