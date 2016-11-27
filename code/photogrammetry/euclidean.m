function Y=euclidean(X)
%EUCLIDEAN Convert homogeneous coordinates to Euclidean.
%
%   EUCLIDEAN(X) converts the (K+1)-by-N array X with K-dimensional
%   homogeneous coordinates to a K-by-N array with the corresponding
%   Euclidean coordinates. X can have multiple layers.
%
%See also: HOMOGENEOUS.


% Divide by element in last row.
Y=X(1:end-1,:,:)./X(repmat(end,1,size(X,1)-1),:,:);
