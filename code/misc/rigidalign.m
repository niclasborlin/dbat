function [T,R,d]=rigidalign(X,Y)
%RIGIDALIGN Compute rigid-body transformation between two point sets.
%
%   [T,R,d]=RIGIDALIGN(X,Y), where X and Y are M-by-N arrays of M-D
%   points, computes the rigid-body transformation that minimizes
%
%      sum  norm( R * X(:,i) + d - Y(:,i) )^2,
%       i
%
%   where R is an M-by-M rotation matrix and d is an M-by-1
%   translation vector. The M+1-by-M+1 matrix T is the homogeneous
%   transformation matrix
%
%       T=[ R  d
%           0  1 ]
%
%   that transforms X to Y.
%
%   References: Soderkvist and Wedin, (1993), Determining the
%   movement of the skeleton using well-configured markers. Journal
%   of Biomechanics 26(12):1473-7.

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

% Compute optimal shift.
d=ym-R*xm;

T=[R,d;zeros(1,m),1];
