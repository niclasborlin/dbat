function [i,j,k,v,COP]=high_op_correlations(s,e,thres,cross)
%HIGH_OP_CORRELATIONS Return list of OP parameters with high correlations.
%
%   [I,J,K,V]=HIGH_OP_CORRELATIONS(S,E,T) returns a list of parameter pairs
%   (I,J) whose estimates for point K have an error correlation V that is
%   higher than the threshold T. The structures S and E are given and
%   returned by BUNDLE, respectively. The I,J indices will be 1..3.
%
%   [I,J,V]=HIGH_OP_CORRELEATIONS(S,E,T,TRUE) also considers cross-camera
%   correlations, i.e. correlations between all point coordinates. In this
%   case, the I,J indices will be 1..3M, where M is the number of points.
%   Warning! May require a lot of memory for large point clouds!
%
%   [...,COP]=... also returns the estimated covariance matrix COP or COPF
%   returned from BUNDLE_COV.
%
%   Example
%      [I,J,K,V]=HIGH_OP_CORRELATIONS(S,E,0.95) returns all intra-point
%      correlations higher than 95%.
%
%See also: BUNDLE, BUNDLE_COV, HIGH_IO_CORRELATIONS, HIGH_OP_CORRELATIONS.

% $Id$

if nargin<4, cross=false; end

if cross
    % Extract OP covariances and correlations.
    COP=bundle_cov(s,e,'COPF');
    COPC=tril(corrmat(COP,true));
    [i,j]=find(abs(COPC)>thres);
    % Output parameters are shifted, i.e. k is v, v is CIO.
    k=full(COPC(sub2ind(size(COPC),i,j)));
    v=COP;
else
    % Ditto but sparse.
    COP=bundle_cov(s,e,'COP');
    COPC=tril(corrmat(COP,true));
    [i,j]=find(abs(COPC)>thres);
    v=full(COPC(sub2ind(size(COPC),i,j)));
    % All index pairs will be on the block diagonal. Renumber them to fall
    % within each camera.
    k=floor((i-1)/size(s.OP,1))+1;
    i=rem(i-1,size(s.OP,1))+1;
    j=rem(j-1,size(s.OP,1))+1;
end
