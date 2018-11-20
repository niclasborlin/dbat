function [i,j,k,v,CIO]=high_io_correlations(s,e,thres,cross)
%HIGH_IO_CORRELATIONS Return list of IO parameters with high correlations.
%
%   [I,J,K,V]=HIGH_IO_CORRELATIONS(S,E,T) returns a list of parameter pairs
%   (I,J) whose estimates for camera K have an error correlation V that is
%   higher than the threshold T. The structures S and E are given to and
%   returned by BUNDLE, respectively. The I,J indices will be 1..N, where N
%   is the number of IO parameters for a single camera.
%
%   [I,J,V]=HIGH_IO_CORRELEATIONS(S,E,T,TRUE) also considers cross-camera
%   correlations, i.e. correlations between IO parameters of all cameras.
%   In this case, I,J will have two columns [R,C], indicating the
%   row and column of each parameter.
%
%   [...,CIO]=... also returns the estimated covariance matrix CIO or CIOF
%   returned from BUNDLE_COV.
%
%   Example
%      [I,J,K,V]=HIGH_IO_CORRELATIONS(S,E,0.95) returns all intra-camera
%      correlations higher than 95%.
%
%See also: BUNDLE, BUNDLE_COV, HIGH_EO_CORRELATIONS, HIGH_OP_CORRELATIONS.


if nargin<4, cross=false; end

if cross
    % Extract IO covariances and correlations.
    CIO=bundle_cov(s,e,'CIOF');
    CIOC=tril(corrmat(CIO,true));
    % Remove any high correlations due to parameter blocks.
    CIOC(~s.IO.struct.leading,:)=0;
    CIOC(:,~s.IO.struct.leading)=0;
    [i,j]=find(abs(CIOC)>thres);
    % Output parameters are shifted, i.e. k is v, v is CIO.
    k=full(CIOC(sub2ind(size(CIOC),i,j)));
    v=CIO;
    % Convert to row, column within IO.
    [ri,ci]=ind2sub(size(s.IO.val),i);
    [rj,cj]=ind2sub(size(s.IO.val),j);
    i=[ri,ci];
    j=[rj,cj];
else
    % Ditto but sparse.
    CIO=bundle_cov(s,e,'CIO');
    CIOC=tril(corrmat(CIO,true));
    [i,j]=find(abs(CIOC)>thres);
    v=full(CIOC(sub2ind(size(CIOC),i,j)));
    % All index pairs will be on the block diagonal. Renumber them to fall
    % within each camera.
    k=floor((i-1)/size(s.IO.val,1))+1;
    i=rem(i-1,size(s.IO.val,1))+1;
    j=rem(j-1,size(s.IO.val,1))+1;
    % Only keep the correlations that correspond to unique cameras.
    [~,ia]=unique(s.IO.struct.block','rows');
    keep=ismember(k,ia);
    i=i(keep);
    j=j(keep);
    k=k(keep);
    v=v(keep);
end
