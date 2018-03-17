function [v,dv,dvn]=brown_tang(u,P,varargin)
%BROWN_TANG Radial scaling function used by DBAT lens distortion functions.
%
%   V=BROWN_TANG(U,P) returns the tangential distortion of Brown
%   (1971) for each 2D point in the 2-by-N array U. The vector P
%   contain the tangential coefficients. P cannot have any length
%   except 1.
%
%   [V,dV]=... also returns a struct dV with the analytical Jacobians
%   with respect to U and P in the field dU and dP, respectively. For
%   more details, see DBAT_BUNDLE_FUNCTIONS.
%
%SEE ALSO: BROWN_RAD, BROWN_DIST_ABS, DBAT_BUNDLE_FUNCTIONS.

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
[pm,pn]=size(P);
if um~=2 || (pn~=1 && ~isempty(P)) || pm==1
    error([mfilename,': bad size']);
end

if isempty(P)
    % Treat this special to get simpler code below.
    v=zeros(size(u));
    if nargout>1
        dv.dP=sparse(2*un,0);
        dv.dU=sparse(2*un,2*un);
    end
    if nargout>2
        dvn.dP=sparse(2*un,0);
        dvn.dU=sparse(2*un,2*un);
    end
    return;
end

%% Actual function code
if nargout<2
    % Only need the function value.
    v=tang_scale(u,P(1:2));
    if length(P)>2
        v=v.*repmat(1+rad_scale(u,P(3:end))',2,1);
    end
else
    [ts,dts]=tang_scale(u,P(1:2));
    v=ts;
    if pm>2
        [rs,drs]=rad_scale(u,P(3:end));
        v=v.*repmat(1+rs',2,1);
    end
end

if nargout>2
    %% Numerical Jacobian

    % FMT is function handle to repackage vector argument to what
    % the function expects.
    if cP
        fun=@(P)feval(mfilename,u,P);
        dvn.dP=jacapprox(fun,P);
    end
    if cU
        fmt=@(u)reshape(u,2,[]);
        fun=@(u)feval(mfilename,fmt(u),P);
        dvn.dU=jacapprox(fun,u);
    end
end

if nargout>1
    %% Analytical Jacobian
    if cP
        dPt=dts.dP;
        if pm<=2
            dP=dPt;
        else
            % Each 2-by-2 block of Pt should be multiplied by 1+rs. 
            % Expand each scalar rs to 2-by-2.
            RS=repmat(rs(reshape(repmat(1:un,2,1),[],1)),1,2);
            dPt=(1+RS).*dPt;
            % Each 2-by-(pn-2) block is ts*drs.dC
            % Expand ts and dC
            TS=repmat(ts(:),1,pm-2);
            dC=drs.dC(reshape(repmat(1:un,2,1),[],1),:);
            dPr=TS.*dC;
            dP=[dPt,dPr];
        end
        dv.dP=dP;
    end
    if cU
        dU=dts.dU;
        if pm>2
            % dU is block-diagonal with 2-by-2 blocks.

            % First, scale each block by 1+rs
            vv=repmat(1+rs',4,1);
            % Row indices.
            ii=repmat([1,2,1,2]',1,un)+repmat(2*(0:un-1),4,1);
            % Column indices.
            jj=repmat([1,1,2,2]',1,un)+repmat(2*(0:un-1),4,1);
            % Build block-diagonal matrix.
            RS=sparse(ii,jj,vv,2*un,2*un);
            % Scale.
            dU=dU.*RS;
            
            % Now, compute ts*drs.dU for each block.
            
            % Expand drs.dU row-wise
            dU21=drs.dU(reshape(repmat(1:un,2,1),[],1),:);
            % Expand ts column-wise.
            vv=repmat(ts,2,1);
            % Build block-diagonal matrix.
            dU22=sparse(ii,jj,vv,2*un,2*un);
            
            dU=dU+dU21.*dU22;
        end
        dv.dU=dU;
    end
end


function fail=selftest(verbose)

% Set up test data.
m=7;
u=rand(2,m);
for l=[0,2,5]
    p=rand(l,1);
    fail=full_self_test(mfilename,{u,p},1e-8,1e-8,verbose);
end
