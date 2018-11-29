function [i,j,k,v,COP]=high_op_correlations(s,e,thres,varargin)
%HIGH_OP_CORRELATIONS Return list of OP parameters with high correlations.
%
%   [I,J,K,V]=HIGH_OP_CORRELATIONS(S,E,T) returns a list of parameter pairs
%   (I,J) whose estimates for point K have an error correlation V that is
%   higher than the threshold T. The structures S and E are given and
%   returned by BUNDLE, respectively. The I,J indices will be 1..3.
%
%   [I,J,K,V]=HIGH_OP_CORRELATIONS(S,E,T,COP) supplies the OP covariance
%   matrix COP.
%
%   [I,J,V]=HIGH_OP_CORRELEATIONS(S,E,T,TRUE) also considers cross-camera
%   correlations, i.e. correlations between all point coordinates. In this
%   case, the I,J indices will be 1..3M, where M is the number of points.
%   Warning! May require a lot of memory for large point clouds!
%
%   [I,J,K,V]=HIGH_OP_CORRELATIONS(S,E,T,TRUE,COPF) supplies the OP
%   covariance matrix COPF.
%
%   [...,COP]=... also returns the estimated covariance matrix COP or COPF
%   returned from BUNDLE_COV.
%
%   Example
%      [I,J,K,V]=HIGH_OP_CORRELATIONS(S,E,0.95) returns all intra-point
%      correlations higher than 95%.
%
%See also: BUNDLE, BUNDLE_COV, HIGH_IO_CORRELATIONS, HIGH_OP_CORRELATIONS.


cross=false;
COP=[];

for i=1:length(varargin)
    if islogical(varargin{i})
        cross=varargin{i};
    else
        COP=varargin{i};
    end
end

if isempty(COP)
    if cross
        COP=bundle_cov(s,e,'COPF');
    else
        COP=bundle_cov(s,e,'COP');
    end
end

% Extract OP covariances and correlations.
COPC=tril(corrmat(COP,true));
[i,j]=find(abs(COPC)>thres);
v=full(COPC(sub2ind(size(COPC),i,j)));

if cross
    % Output parameters are shifted, i.e. k is v, v is CIO.
    k=v;
    v=COP;
else
    % All index pairs will be on the block diagonal. Renumber them to be
    % within each camera.
    k=floor((i-1)/size(s.OP.val,1))+1;
    i=rem(i-1,size(s.OP.val,1))+1;
    j=rem(j-1,size(s.OP.val,1))+1;
end
