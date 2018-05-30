function ixOP=createopcolumnindices(cOP);
% Create arrays of columns indices for OP derivatives.
% A zero element means that the partial derivative should not be
% calculated/stored.


% How many points do we have?
nObjs=size(cOP,2);

ixOP=reshape(cumsum(cOP(:)),3,nObjs).*cOP;
