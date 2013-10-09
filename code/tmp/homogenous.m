function Y=homogenous(X)
%HOMOGENOUS Convert Euclidean coordinates to homogenous.
%
%   HOMOGENOUS(X) converts the K-by-N array X with K-dimensional
%   Euclidean coordinates to a (K+1)-by-N array with the corresponding
%   homogenous coordinates. X can have multiple layers.
%
%See also: EUCLIDEAN.

% $Id: e5c2450808b5148e865a913f5453dbdf106f7c6b $

% Append unity.
Y=X;
Y(end+1,:)=1;
