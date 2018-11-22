function [ixpp,ixf,ixK,ixP,ixa,ixu]=createiocolumnindices(cIO,nK,nP)
% Create arrays of columns indices for IO derivatives.
% A zero element means that the partial derivative should not be
% calculated/stored.


% How many cameras do we have?
nCams=size(cIO,2);

ix=reshape(cumsum(cIO(:)),5+nK+nP,nCams).*cIO;

ixf=ix(1,:);
ixpp=ix(2:3,:);
ixa=ix(4:5,:);
ixK=ix(5+(1:nK),:);
ixP=ix(5+nK+(1:nP),:);

