function [T,dAng,dP,dC]=pm_roteuler1(ang,ax,P,C,cAng,cP,cC)
%PM_ROTEULER1 Camera model point rotation using Euler angles.
%
%[T,dAng,dP,dC]=pm_roteuler1(ang,ax,P,C[,cAng,cP,cC])
%ang  - Euler angles.
%ax   - Euler axis sequence.
%P    - 3xn matrix of object points [X;Y;Z].
%C    - 3-vector with camera center [XC;YC;ZC].
%c[x] - calculate partial derivative w.r.t. x? cP may be matrix for
%       elementwise indication.
%T    - 3xn matrix of transformed object points [U;V;W].
%d[x] - jacobian of T w.r.t. x.

if (nargin<5), cAng=(nargout>1); end
if (nargin<6), cP=(nargout>2); end
if (nargin<7), cC=(nargout>3); end

dAng=[];
dP=[];
dC=[];

% Number of points.
n=size(P,2);
% Subtract camera center.
deltaP=P-repmat(C,1,n);

[M,da1,da2,da3]=pm_eulerrotmat(ang,ax);

T=M*deltaP;

if (cAng)
    dA1=da1*deltaP;
    dA2=da2*deltaP;
    dA3=da3*deltaP;
    dAng=[dA1(:),dA2(:),dA3(:)];
end

if (any(cP(:)))
    dP=kron(speye(n),M);
    if (length(cP)>1 & any(~cP(:)))
        % Remove unwanted columns.
        dP=dP(:,cP(:));
    end
end

if (cC)
    dC=kron(ones(n,1),-M);
end
