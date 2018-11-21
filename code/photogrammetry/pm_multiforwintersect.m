function [OP,R]=pm_multiforwintersect(IO,EO,ei,colPos,pts,i)
%PM_MULTIFORWINTERSECT Camera network forward intersection.
%
%   OP=PM_MULTIFORWINTERSECT(IO,EO,EI,COLPOS,PTS,I) calculates object points
%   by forward intersection. The 6-by-N array EO contains the external
%   orientation of each camera. The N-array (or scalar) EI contains indices
%   into the IO array. The M-by-N array COLPOS contains column indices into
%   the 2-by-P array PTS of observations. The Q-vector I contains the
%   indices of which points to calculate. The estimated points are returned
%   in the 3-by-Q array OP.
%
%   [OP,R]=... also returns the squared residual in R, normalized by the
%   number of rays.

if any(~isfinite(EO(:))), error('Bad or uninitialized EO data'); end
if any(~isfinite(IO(:))), error('Bad or uninitialized IO data'); end
    
if isscalar(ei)
    ei=repmat(ei,1,size(EO,2));
end

% Construct camera matrices for active cameras only.
P=nan(3,4,size(EO,2));
for j=find(any(colPos(i,:),1))
    RR=pm_eulerrotmat(EO(4:6,j));
    CC=EO(1:3,j);
    K=[-IO(1,ei(j))*eye(2),IO(2:3,ei(j));0,0,1];
    P(:,:,j)=K*RR*[eye(3),-CC];
end

% Pre-allocate points and residuals.
OP=nan(3,length(i));
R=nan(1,length(i));
% Find out which camera combinations we have.
[camComb,~,ui]=unique(colPos(i,:)~=0,'rows');

% Do one forward intersection for each camera combination.
for ii=1:size(camComb,1)
    % Active cameras.
    camIx=find(camComb(ii,:));
    if nnz(camIx)>1
        % Which object point have this camera combination?
        objIx=i(ui==ii);
        % Corresponding measurement columns.
        colIx=full(colPos(objIx,camIx));
        % Extract measurements.
        xy=reshape(pts(:,colIx'),2,nnz(camIx),nnz(objIx));
        % Do forward intersection.
        [OP(:,ui==ii),R(ui==ii)]=pm_forwintersect3(P(:,:,camIx),xy);
    end
end
