function [ixpp,ixf,ixK,ixP,ixa,ixu]=createiocolumnindices(cIO,nK,nP);
% Create arrays of columns indices for IO derivatives.
% A zero element means that the partial derivative should not be
% calculated/stored.


% How many cameras do we have?
nCams=size(cIO,2);

ix=reshape(cumsum(cIO(:)),11+nK+nP,nCams).*cIO;

ixpp=ix(1:2,:);
ixf=ix(3,:);
ixK=ix(3+(1:nK),:);
ixP=ix(3+nK+(1:nP),:);
ixa=ix(3+nK+nP+(1:2),:);
ixu=ix(3+nK+nP+6+(1:2),:);
