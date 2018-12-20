function [T,R,d,alpha]=rigidalign(X,Y,scale)
%RIGIDALIGN Compute rigid-body transformation between two point sets.
%
%   [T,R,d,ALPHA]=RIGIDALIGN(X,Y,TRUE), where X and Y are M-by-N
%   arrays of M-D points, computes the rigid-body transformation that
%   minimizes
%
%      sum  norm( alpha * R * X(:,i) + d - Y(:,i) )^2,
%       i
%
%   where R is an M-by-M rotation matrix, d is an M-by-1 translation
%   vector and alpha is a scale factor. The M+1-by-M+1 matrix T is the
%   homogeneous transformation matrix
%
%       T=[ alpha * R  d
%                   0  1 ]
%
%   that transforms X to Y.
%
%   RIGIDALIGN(X,Y) or RIGIDALIGN(X,Y,FALSE) assumes ALPHA=1.
%
%   References: Soderkvist and Wedin, (1993), Determining the
%   movement of the skeleton using well-configured markers. Journal
%   of Biomechanics 26(12):1473-7.

if nargin<3, scale=false; end

% Verify sizes.
if any(size(X)~=size(Y))
    error('Input size mismatch: X and Y must have the same size.');
end

[m,n]=size(X);

% Compute center of mass for point clouds.
xm=mean(X,2);
ym=mean(Y,2);

% Shift point clouds to have center of mass at origin.
A=X-repmat(xm,1,n);
B=Y-repmat(ym,1,n);

C=B*A';

% Perform singular value decomposition of C.
[P,~,Q]=svd(C);

% Compute optimal rotation matrix.
R=P*diag([ones(1,m-1),det(P*Q')])*Q';

% Compute scaling factor.
if scale
    gamma=trace((R*A)'*B);
    alpha=gamma/trace(A'*A);
else
    alpha=1;
end
% Compute optimal shift.
d=ym-alpha*R*xm;

T=[alpha*R,d;zeros(1,m),1];
