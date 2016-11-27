function d=pm_multidepth(IO,EO,OP,vis,camNo)
%PM_MULTIDEPTH Object point depth for complete network.
%
%   D=pm_multidepth(IO,EO,OP,vis,camNo)
%   IO    - 7+nK+nP-by-K interior orientation [pp;f;K;P;a;u].
%   EO    - 7-by-N exterior orientation [C;ang;0]
%   OP    - 3-by-M array with object points.
%   vis   - M-by-N sparse visibility matrix.
%   camNo - scalar or N-vector with column indices into IO for each EO.
%   D     - M-by-N array with depths for each OP w.r.t. each camera is it
%           visible in.


[m,n]=size(vis);

if isscalar(camNo), camNo=repmat(camNo,1,n); end

% First create camera calibration matrices.
KK=zeros(3,3,size(IO,2));
for i=1:size(IO,2)
    KK(:,:,i)=[-IO(3,i)*eye(2),IO(1:2,i);0,0,1];
end

% Next, create camera matrices.
PP=zeros(3,4,n);
for i=1:n
    if all(~isnan(EO(:,i)))
        K=KK(:,:,camNo(i));
        PP(:,:,i)=K*rotmat(EO(4:6,i))*[eye(3),-EO(1:3,i)];
    end
end

d=nan(m,n);
% For each camera with a visible point.
for i=find(any(vis,1))
    j=vis(:,i);
    d(j,i)=-ptdepth(PP(:,:,i),OP(:,j));
end
