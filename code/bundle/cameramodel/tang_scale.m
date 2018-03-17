function [v,dv,dvn]=tang_scale(u,p,varargin)
%TANG_SCALE Radial scaling function used by DBAT lens distortion functions.
%
%   V=TANG_SCALE(U,P) returns the tangential scaling between the 2-vector P
%   and 2-by-N point array U. Each column V(:,I) is
%   V(I,:)=U(:,I)'*U(:,I)*P+2*U(:,I)*U(:,I)'*P.
%
%   [V,dV]=... also returns a struct dV with the analytical Jacobians
%   with respect to U and P in the field dU and dP, respectively. For
%   more details, see DBAT_BUNDLE_FUNCTIONS.
%
%SEE ALSO: DBAT_BUNDLE_FUNCTIONS

% Treat selftest call separately.
if nargin>=1 && ischar(u), v=selftest(nargin>1 && p); return; end

% Otherwise, verify number of parameters.
narginchk(1,4);

v=[]; %#ok<NASGU>
dv=[];
dvn=[];

if nargout>1
    % Construct empty Jacobian struct.
    dv=struct('dP',[],...
              'dU',[]);
    dvn=dv;
end

% What Jacobians to compute?
cU=nargout>1 && (length(varargin)<1 || varargin{1});
cP=nargout>1 && (length(varargin)<2 || varargin{2});

%% Test parameters
[um,un]=size(u);
if um~=2 || any(size(p)~=[2,1])
    error([mfilename,': bad size']);
end

%% Actual function code
% v=u'*u*p + 2*u*u'*p
uTu=sum(u.^2,1);
pTu=p'*u;
v=p*uTu+2*repmat(pTu,2,1).*u;

if nargout>2
    %% Numerical Jacobian

    % FMT is function handle to repackage vector argument to what
    % the function expects.
    if cP
        fun=@(p)feval(mfilename,u,p);
        dvn.dP=jacapprox(fun,p);
    end
    if cU
        fmt=@(u)reshape(u,2,[]);
        fun=@(u)feval(mfilename,fmt(u),p);
        dvn.dU=jacapprox(fun,u);
    end
end

if nargout>1
    %% Analytical Jacobian
    if cP
        % Each block is u'*u*eye(2) + 2*u*u'
        u12=u(1,:).^2;
        u22=u(2,:).^2;
        u1u2=u(1,:).*u(2,:);
        % This creates [dP(1,1); dP(1,2); dP(2,1); dP(2,2)]
        dPT=[uTu+2*u12; 2*u1u2; 2*u1u2; uTu+2*u22];
        % Reshape to 2-by-2 blocks and transpose.
        dv.dP=reshape(dPT,2,[])';
    end
    if cU
        % Each block is 2*(p*u'+u*p'+p'u*eye(2))
        p1u1=p(1)*u(1,:);
        p1u2=p(1)*u(2,:);
        p2u1=p(2)*u(1,:);
        p2u2=p(2)*u(2,:);
        % Values in memory order.
        vv=reshape(2*[2*p1u1+pTu;p1u2+p2u1;p1u2+p2u1;2*p2u2+pTu],2,[]);
        % Row indices.
        ii=repmat([1,2,1,2]',1,un)+repmat(2*(0:un-1),4,1);
        % Column indices.
        jj=repmat([1,1,2,2]',1,un)+repmat(2*(0:un-1),4,1);
        dv.dU=sparse(ii,jj,vv,2*un,2*un);
    end
end


function fail=selftest(verbose)

% Set up test data.
p=rand(2,1);
m=7;
u=rand(2,m);

fail=full_self_test(mfilename,{u,p},1e-8,1e-8,verbose);

