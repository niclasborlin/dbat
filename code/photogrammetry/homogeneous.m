function Y=homogeneous(X)
%HOMOGENEOUS Convert Euclidean coordinates to homogeneous.
%
%   HOMOGENEOUS(X) converts the K-by-N array X with K-dimensional Euclidean
%   coordinates to a (K+1)-by-N array with the corresponding homogeneous
%   coordinates. X can have multiple layers.
%
%See also: EUCLIDEAN.


% Append unity.
Y=X;
Y(end+1,:)=1;
