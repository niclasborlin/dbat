function [i,j,k,v,CEO]=high_eo_correlations(s,e,thres,cross)
%HIGH_EO_CORRELATIONS Return list of EO parameters with high correlations.
%
%   [I,J,K,V]=HIGH_EO_CORRELATIONS(S,E,T) returns a list of parameter pairs
%   (I,J) whose estimates for camera station K have an error correlation V
%   that is higher than the threshold T. The structures S and E are given to
%   and returned by BUNDLE, respectively. The I,J indices will be 1..N,
%   where N is the number of EO parameters for a single camera station.
%
%   [I,J,V]=HIGH_EO_CORRELEATIONS(S,E,T,TRUE) also considers cross-camera
%   station correlations, i.e. correlations between EO parameters of all
%   camera stations. In this case, the I,J indices will be 1..M*N, where M
%   is the number of camera stations.
%
%   [...,CEO]=... also returns the estimated covariance matrix CEO or CEOF
%   returned from BUNDLE_COV.
%
%   Example
%      [I,J,K,V]=HIGH_EO_CORRELATIONS(S,E,0.95) returns all intra-camera
%      station correlations higher than 95%.
%
%See also: BUNDLE, BUNDLE_COV, HIGH_OP_CORRELATIONS, HIGH_OP_CORRELATIONS.


if nargin<4, cross=false; end

if cross
    % Extract EO covariances and correlations.
    CEO=bundle_cov(s,e,'CEOF');
    CEOC=tril(corrmat(CEO,true));
    [i,j]=find(abs(CEOC)>thres);
    % Output parameters are shifted, i.e. k is v, v is CEO.
    k=full(CEOC(sub2ind(size(CEOC),i,j)));
    v=CEO;
    % Only keep the correlations that correspond to unique cameras.
    [~,ia]=unique(s.EO.struct.block','rows');
    keep=ismember(k,ia);
    i=i(keep);
    j=j(keep);
    k=k(keep);
    v=v(keep);
else
    % Ditto but sparse.
    CEO=bundle_cov(s,e,'CEO');
    CEOC=tril(corrmat(CEO,true));
    [i,j]=find(abs(CEOC)>thres);
    v=full(CEOC(sub2ind(size(CEOC),i,j)));
    % All index pairs will be on the block diagonal. Renumber them to fall
    % within each camera.
    k=floor((i-1)/(size(s.EO.val,1)))+1;
    i=rem(i-1,size(s.EO.val,1))+1;
    j=rem(j-1,size(s.EO.val,1))+1;
    % Only keep the correlations that correspond to unique cameras.
    [~,ia]=unique(s.EO.struct.block','rows');
    keep=ismember(k,ia);
    i=i(keep);
    j=j(keep);
    k=k(keep);
    v=v(keep);
end
