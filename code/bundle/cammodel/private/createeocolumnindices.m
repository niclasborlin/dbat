function [ixC,ixAng]=createeocolumnindices(cEO);
% Create arrays of columns indices for EO derivatives.
% A zero element means that the partial derivative should not be
% calculated/stored.

% How many photos do we have?
nPhotos=size(cEO,2);

ix=reshape(cumsum(cEO(:)),6,nPhotos).*cEO;

ixC=ix(1:3,:);
ixAng=ix(4:6,:);
