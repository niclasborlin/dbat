function [ixC,ixAng]=createeocolumnindices(cEO);
% Create arrays of columns indices for EO derivatives.
% A zero element means that the partial derivative should not be
% calculated/stored.

% $Id$

ix=reshape(cumsum(cEO(:)),size(cEO));

ixC=ix(1:3,:);
ixAng=ix(4:end-1,:);
