function Y=euclidean(X)
%EUCLIDEAN Convert homogenous coordinates to Euclidean.
%
%   EUCLIDIAN(X) converts the (K+1)-by-N array X with K-dimensional
%   homogenous coordinates to a K-by-N array with the corresponding
%   Euclidean coordinates. X can have multiple layers.
%
%See also: HOMOGENOUS.

% $Id: 499bb09231193bac991cec397d38c9370a674c23 $

% Divide by element in last row.
Y=X(1:end-1,:,:)./X(repmat(end,1,size(X,1)-1),:,:);
