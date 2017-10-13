function [p,dp0,df,dP]=proj(p0,f,P,cp0,cf,cP)
%PROJ Projection of points in the camera coordinate system.
%
%[p,dp0,df,dP]=proj(p0,f,P[,cp0,cf,cP])
%p0   - principal point [xp;yp].
%f    - focal length.
%P    - 3xn matrix of object points in the camera coordinate system.
%c[x] - should we calculate partial derivative w.r.t. x?
%p    - projected points.
%d[x] - jacobian w.r.t x.

if nargin<4, cp0=(nargout>1); end
if nargin<5, cf=(nargout>2); end
if nargin<6, cP=(nargout>3); end

dp0=[];
df=[];
dP=[];

% Number of points.
n=size(P,2);

U=P(1,:);
V=P(2,:);
W=P(3,:);

UVdivW=[U./W;V./W];
fUVdivW=f*UVdivW;

p=repmat(p0,1,n)-fUVdivW;

if cp0
    dp0=repmat(eye(2),n,1);
end

if cf
    df=-UVdivW(:);
end

if cP
    v12=-f./W;
    v3=fUVdivW./[W;W];

    vv=reshape([v12;v12;v3],2,[]);
    i0=reshape(repmat(0:2:2*n-1,2,1),1,[]);
    ii=[i0+1;i0+2];
    jj=repmat(0:3:3*n-1,4,1)+repmat([1;2;3;3],1,n);
    dP=sparse(ii,jj,vv,2*n,3*n);
end
