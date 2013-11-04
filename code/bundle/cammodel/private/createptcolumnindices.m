function ixp=createptcolumnindices(cp);
% Create arrays of columns indices for OP derivatives.
% A zero element means that the partial derivative should not be
% calculated/stored.

cp=[cp;cp];

% How many points do we have?
nPts=size(cp,2);

ixp=reshape(cumsum(cp(:)),2,nPts).*cp;
