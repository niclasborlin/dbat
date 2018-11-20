function [delta,dp,dp0,dK,dP]=pm_lens1(p,p0,K,P,cp,cp0,cK,cP)
%PM_LENS1 Lens distortion for camera model.
%
%[delta,dp,dp0,dK,dP]=pm_lens1(p,p0,K,P[,cp,cp0,cK,cP])
%p     - 2xn matrix of projected points.
%p0    - [xp;yp] principal point.
%K     - radial distortion coefficients. May be of any length.
%P     - tangential distortion coefficients. May be of length 0, 2, or 3.
%c[x]  - should we calculate partial derivative w.r.t. x?
%delta - [deltax;deltay] shift caused by lens distortion.
%d[x]  - jacobian w.r.t. x.

% v1.0  2003-10-21. Niclas Borlin, niclas@cs.umu.se.

if (nargin<5), cp=(nargout>1); end
if (nargin<6), cp0=(nargout>1); end
if (nargin<7), cK=(nargout>2); end
if (nargin<8), cP=(nargout>3); end

if (isempty(K) && isempty(P))
	delta=sparse(size(p,1),size(p,2));
	dp=sparse(numel(p),numel(p));
	dp0=sparse(numel(p),2);
	dK=sparse(numel(p),0);
	dP=sparse(numel(p),0);
	return;
end

dp=[];
dp0=[];
dK=[];
dP=[];

% Number of points.
n=size(p,2);

% Difference between x,y and principal point.
xBar=p(1,:)-p0(1);
yBar=p(2,:)-p0(2);

% Radial distance squared.
r2=xBar.^2+yBar.^2;

% Create r2.^[1..nK], where nK is the number of K coefficients.

% Create r2 exponent matrix.
nK=length(K);
r2e=repmat((1:nK)',1,n);
			
r2k=repmat(r2,nK,1).^r2e;

% Inner product K(1)*r^2+K(2)*r^4+K(3)*r^6+...
Kr=(r2k'*K)';

deltaRx=xBar.*Kr;
deltaRy=yBar.*Kr;

if (isempty(P))
	deltaTx=0;
	deltaTy=0;
else
	if (length(P)>2)
		P3=P(3);
	else
		P3=0;
	end
	deltaTx=(P(1)*(r2+2*xBar.^2)+2*P(2)*xBar.*yBar)*(1+P3);
	deltaTy=(P(2)*(r2+2*yBar.^2)+2*P(1)*xBar.*yBar)*(1+P3);
end

delta=[deltaRx+deltaTx;
	   deltaRy+deltaTy];

if (any([cp,cp0]) && ~isempty(K))
	% Inner product K(1)+2*K(2)*r^2+3*K(3)*r^4+...
	kK=(1:nK)'.*K;
	Krm1=([ones(1,n);r2k(1:end-1,:)]'*kK)';
else
	Krm1=[];
end

if (cp)
	if (isempty(K))
		drxdx=0;
		drxdy=0;
		drydx=0;
		drydy=0;
	else
		drxdx=Kr+2*xBar.^2.*Krm1;
		drxdy=2*xBar.*yBar.*Krm1;
		drydx=drxdy;
		drydy=Kr+2*yBar.^2.*Krm1;
	end
	
	if (isempty(P))
		dtxdx=0;
		dtxdy=0;
		dtydx=0;
		dtydy=0;
	else
		if (length(P)>2)
			P3=P(3);
		else
			P3=0;
		end
		dtxdx=(6*P(1)*xBar+2*P(2)*yBar)*(1+P3);
		dtxdy=(2*P(1)*yBar+2*P(2)*xBar)*(1+P3);
		dtydx=dtxdy;
		dtydy=(6*P(2)*yBar+2*P(1)*xBar)*(1+P3);
	end
	
	i1=1:2:2*n;
	j1=1:2:2*n;
	v1=drxdx+dtxdx;
	
	i2=i1+1;
	j2=j1;
	v2=drydx+dtydx;
	
	i3=i1;
	j3=j1+1;
	v3=drxdy+dtxdy;
	
	i4=i1+1;
	j4=j1+1;
	v4=drydy+dtydy;
	
	dp=sparse([i1;i2;i3;i4],[j1;j2;j3;j4],[v1;v2;v3;v4],2*n,2*n);
end

if (cp0)
	if (isempty(K))
		drxdxp=0;
		drxdyp=0;
		drydxp=0;
		drydyp=0;
	else
		drxdxp=-(Kr+2*xBar.^2.*Krm1);
		drxdyp=-(2*xBar.*yBar.*Krm1);
		drydxp=drxdyp;
		drydyp=-(Kr+2*yBar.^2.*Krm1);
	end
	
	if (isempty(P))
		dtxdxp=0;
		dtxdyp=0;
		dtydxp=0;
		dtydyp=0;
	else
		if (length(P)>2)
			P3=P(3);
		else
			P3=0;
		end
		dtxdxp=-(6*P(1)*xBar+2*P(2)*yBar)*(1+P3);
		dtxdyp=-(2*P(1)*yBar+2*P(2)*xBar)*(1+P3);
		dtydxp=dtxdyp;
		dtydyp=-(6*P(2)*yBar+2*P(1)*xBar)*(1+P3);
	end
	
	ddxp=[drxdxp+dtxdxp;drydxp+dtydxp];
	ddyp=[drxdyp+dtxdyp;drydyp+dtydyp];
	dp0=[ddxp(:),ddyp(:)];
end

if (cK)
	dxdK=(repmat(xBar,nK,1).*r2k)';
	dydK=(repmat(yBar,nK,1).*r2k)';
	dxy=cat(3,dxdK,dydK);
	dK=reshape(permute(dxy,[3,1,2]),2*n,nK);
end

if (cP)
	if (isempty(P))
		dP=zeros(2*n,0);
	else
		if (length(P)>2)
			P3=P(3);
			dtxdP3=P(1)*(r2+2*xBar.^2)+2*P(2)*xBar.*yBar;
			dtydP3=P(2)*(r2+2*yBar.^2)+2*P(1)*xBar.*yBar;
		else
			P3=0;
			dtxdP3=zeros(n,0);
			dtydP3=zeros(n,0);
		end
		dtxdP1=(r2+2*xBar.^2)*(1+P3);
		dtxdP2=2*xBar.*yBar*(1+P3);
		dtydP1=dtxdP2;
		dtydP2=(r2+2*yBar.^2)*(1+P3);
		dtdP1=[dtxdP1;dtydP1];
		dtdP2=[dtxdP2;dtydP2];
		dtdP3=[dtxdP3;dtydP3];
		dP=[dtdP1(:),dtdP2(:),dtdP3(:)];
	end
end
