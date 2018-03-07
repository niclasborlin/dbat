function [Q,dQ,dQn]=affine2(b,P,varargin)
%AFFINE2 2D affine transform for the DBAT projection model
%
%   Q=AFFINE2(B,P) applies the linear transform AFFINE2MAT(B) to the
%   2D points in the 2-by-N array P, i.e. computes AFFINE2MAT(B)*P.
%
%   [Q,dQ]=... also returns a struct dQ with the analytical Jacobians
%   with respect to B and P in the fields dB and dP. For more details,
%   see DBAT_BUNDLE_FUNCTIONS.
%
%SEE ALSO: DBAT_BUNDLE_FUNCTIONS

% Treat selftest call separately.
if nargin>=1 && ischar(b), Q=selftest(nargin>1 && P); return; end

% Otherwise, verify number of parameters.
narginchk(1,4);

Q=[];
dQ=[];
dQn=[];

if nargout>1
    % Construct empty Jacobian struct.
    dQ=struct('dP',[],...
              'dB',[]);
    dQn=dQ;
end

% What Jacobians to compute?
cB=nargout>1 && (length(varargin)<1 || varargin{1});
cP=nargout>1 && (length(varargin)<2 || varargin{2});

%% Test parameters
[m,n]=size(P);
if m~=2 || any(size(b)~=[2,1])
    error([mfilename,': bad size']);
end

%% Actual function code
A=affine2mat(b);
Q=A*P;

if nargout>2
    %% Numerical Jacobian

    % FMT is function handle to repackage vector argument to what
    % the function expects.
    if cP
        fmt=@(P)reshape(P,2,[]);
        fun=@(P)feval(mfilename,b,fmt(P));
        dQn.dP=jacapprox(fun,P);
    end
    if cB
        fun=@(b)feval(mfilename,b,P);
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
        dQ.dP=sparse(ii,jj,vv,2*n,2*n);
    end
    if cB
        dQ.dB=zeros(2*n,2);
        dQ.dB(1:2:end,:)=P';
    end    
end


function fail=selftest(verbose)

% Set up test data.
m=2;
n=5;
b=rand(m,1);
P=rand(m,n);

fail=full_self_test(mfilename,{b,P},1e-8,1e-8,verbose);
