function [p,dp0,df,dP,dC,dAng]=eulerpinhole(p0,f,P,C,ang,ax,cp0,cf,cP,cC,cAng)
%EULERPINHOLE General pinhole camera projection with Euler angles.
%
%[p,dp0,df,dP,dC,dAng]=eulerpinhole(p0,f,P,C,ang,ax[,cp0,cf,cP,cC,cAng])
%p0   - principal point [xp;yp].
%f    - focal length.
%P    - 3xn matrix of object points in global coordinate system.
%C    - 3-vector with camera center.
%ang  - Euler angles.
%ax   - Euler axis sequence.
%c[x] - calculate partial derivative w.r.t. x? cP may be matrix for
%       elementwise indication.
%p    - projected points.
%d[x] - jacobian w.r.t [x].

if nargin<7, cp0=(nargout>1); end
if nargin<8, cf=(nargout>2); end
if nargin<9, cP=(nargout>3); end
if nargin<10, cC=(nargout>4); end
if nargin<11, cAng=(nargout>5); end

dp0=[];
df=[];
dP=[];
dC=[];
dAng=[];

if nargout==1
    % Rotate point to camera coordinate system.
    T=roteuler(ang,ax,P,C);

    % Project rotated point into camera.
    p=proj(p0,f,T);
else
    % Rotate point to camera coordinate system.
    [T,dTdAng,dTdP,dTdC]=roteuler(ang,ax,P,C,cAng,cP,cC);

    % Project rotated point into camera.
    [p,dpdp0,dpdf,dpdT]=proj(p0,f,T,cp0,cf,any(cP(:)) | cC | cAng);
    
    if cp0
        dp0=dpdp0;
    end
    if cf
        df=dpdf;
    end
    if any(cP(:))
        dP=dpdT*dTdP;
    end
    if cC
        dC=dpdT*dTdC;
    end
    if cAng
        dAng=dpdT*dTdAng;
    end
end
