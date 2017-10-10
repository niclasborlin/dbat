function T=rigidalign(X,Y)
%RIGIDALIGN Compute rigid-body transformation between two point sets.
%
%   T=RIGIDALIGN(X,Y), where X and Y are M-by-N arrays of M-D
%   points, computes the rigid-body transformation the minimizes
%
%      sum  norm( R * X(:,i) + d - Y(:,i) )^2,
%       i
%
%   where R is an M-by-M rotation matrix and d is an M-by-1
%   translation vector.
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



