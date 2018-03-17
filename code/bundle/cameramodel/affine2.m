function [Q,dQ,dQn]=affine2(U,b,varargin)
%AFFINE2 2D affine transform for the DBAT projection model
%
%   Q=AFFINE2(U,B) applies the linear transform AFFINE2MAT(B) to the
%   2D points in the 2-by-N array U, i.e. computes AFFINE2MAT(B)*U.
%
%   [Q,dQ]=... also returns a struct dQ with the analytical Jacobians
%   with respect to B and U in the fields dB and dU. For more details,
%   see DBAT_BUNDLE_FUNCTIONS.
%
%SEE ALSO: AFFINE2MAT, DBAT_BUNDLE_FUNCTIONS

% Treat selftest call separately.
if nargin>=1 && ischar(U), Q=selftest(nargin>1 && b); return; end

% Otherwise, verify number of parameters.
narginchk(2,4);

Q=[]; %#ok<NASGU>
dQ=[];
dQn=[];

if nargout>1
    % Construct empty Jacobian struct.
    dQ=struct('dU',[],...
              'dB',[]);
    dQn=dQ;
end

% What Jacobians to compute?
cP=nargout>1 && (length(varargin)<1 || varargin{1});
cB=nargout>1 && (length(varargin)<2 || varargin{2});

%% Test parameters
[m,n]=size(U);
if m~=2 || any(size(b)~=[2,1])
    error([mfilename,': bad size']);
end

%% Actual function code
A=affine2mat(b);
Q=A*U;

if nargout>2
    %% Numerical Jacobian

    % FMT is function handle to repackage vector argument to what
    % the function expects.
    if cP
        fmt=@(U)reshape(U,2,[]);
        fun=@(U)feval(mfilename,fmt(U),b);
        dQn.dU=jacapprox(fun,U);
    end
    if cB
        fun=@(b)feval(mfilename,U,b);
        dQn.dB=jacapprox(fun,b);
    end
end

if nargout>1
    %% Analytical Jacobian
    if cP
        % Each block is A = [1+B(1), B(2); 0 1].

        % Values for Jacobian in memory-order.
        vv=repmat(A([1,3,4])',1,n);
        % Row indices.
        ii=repmat(0:2:2*n-1,3,1)+repmat([1,1,2]',1,n);
        % Column indices.
        jj=repmat(0:2:2*n-1,3,1)+repmat([1,2,2]',1,n);
        dQ.dU=sparse(ii,jj,vv,2*n,2*n);
    end
    if cB
        dQ.dB=zeros(2*n,2);
        dQ.dB(1:2:end,:)=U';
    end    
end


function fail=selftest(verbose)

% Set up test data.
m=2;
n=5;
b=rand(m,1);
U=rand(m,n);

fail=full_self_test(mfilename,{U,b},1e-8,1e-8,verbose);
