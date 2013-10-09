function [p,dp0,df,dP]=pm_proj1(p0,f,P,cp0,cf,cP)
%PM_PROJ1 Projection of points in the camera coordinate system.
%
%[p,dp0,df,dP]=pm_proj1(p0,f,P[,cp0,cf,cP])
%p0   - principal point [xp;yp].
%f    - focal length.
%P    - 3xn matrix of object points in the camera coordinate system.
%c[x] - should we calculate partial derivative w.r.t. x?
%p    - projected points.
%d[x] - jacobian w.r.t x.

% v1.0  2003-10-21. Niclas Borlin, niclas@cs.umu.se.

if (nargin<4), cp0=(nargout>1); end
if (nargin<5), cf=(nargout>2); end
if (nargin<6), cP=(nargout>3); end

dp0=[];
df=[];
dP=[];

% Number of points.
n=size(P,2);

U=P(1,:);
V=P(2,:);
W=P(3,:);

UVdivW=[U./W;V./W];

p=repmat(p0,1,n)-f*UVdivW;

if (cp0)
	dp0=kron(ones(n,1),speye(2));
end

if (cf)
	df=-UVdivW(:);
end

if (cP)
	i1=1:2:2*n;
	j1=1:3:3*n;
	v1=-f./W;
	
	i2=i1+1;
	j2=j1+1;
	v2=v1;
	
	i3=1:2*n;
	j3=kron(3:3:3*n,ones(1,2));
	v3=f*UVdivW./[W;W];
	
	dP=sparse([i1,i2,i3],[j1,j2,j3],[v1(:);v2(:);v3(:)],2*n,3*n);
end