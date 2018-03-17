function [v,dv,dvn]=rad_scale(u,c,varargin)
%RAD_SCALE Radial scaling function used by DBAT lens distortion functions.
%
%   V=RAD_SCALE(U,C) returns the radial scaling between the M-vector C
%   with coefficients 2-by-N point array U. Each element V(I) is
%   V(I)=C(1)*R(U)^2 + C(2)*R(U)^4 + ....
%
%   [V,dV]=... also returns a struct dV with the analyticals Jacobian
%   with respect to C and U in the field dC and dU, respectively. For
%   more details, see DBAT_BUNDLE_FUNCTIONS.
%
%SEE ALSO: LENS_RAD2, POWER_VEC, DBAT_BUNDLE_FUNCTIONS

% Treat selftest call separately.
if nargin>=1 && ischar(u), v=selftest(nargin>1 && c); return; end

% Otherwise, verify number of parameters.
narginchk(1,4);

v=[]; %#ok<NASGU>
dv=[];
dvn=[];

if nargout>1
    % Construct empty Jacobian struct.
    dv=struct('dC',[],...
              'dU',[]);
    dvn=dv;
end

% What Jacobians to compute?
cU=nargout>1 && (length(varargin)<1 || varargin{1});
cC=nargout>1 && (length(varargin)<2 || varargin{2});

%% Test parameters
[~,cn]=size(c);
[um,un]=size(u);
if cn~=1 || um~=2
    error([mfilename,': bad size']);
end

%% Actual function code
if nargout<2
    r2=lens_rad2(u);
    pv=power_vec(r2',length(c));
    v=pv'*c;
else
    [r2,dr2]=lens_rad2(u);
    [pv,dpv]=power_vec(r2',length(c));
    v=pv'*c;
end    

if nargout>2
    %% Numerical Jacobian

    % FMT is function handle to repackage vector argument to what
    % the function expects.
    if cC
        fun=@(c)feval(mfilename,u,c);
        dvn.dC=jacapprox(fun,c);
    end
    if cU
        fmt=@(u)reshape(u,2,[]);
        fun=@(u)feval(mfilename,fmt(u),c);
        dvn.dU=jacapprox(fun,u);
    end
end

if nargout>1
    %% Analytical Jacobian
    if cC
        dv.dC=pv';
    end
    if cU
        dv.dU=kron(speye(un),c')*dpv.dX*dr2.dU;
    end
end


function fail=selftest(verbose)

% Set up test data.
c=rand(5,1);
m=7;
u=rand(2,m);

fail=full_self_test(mfilename,{u,c},1e-8,1e-8,verbose);

