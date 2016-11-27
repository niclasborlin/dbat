function [OP,r]=pm_forwintersect3(P,xy)
%PM_FORWINTERSECT3 Estimate world coordinates by forward intersection.
%
%[OP,r]=PM_FORWINTERSECT3(P,xy)
%P      - 3-by-4-by-N array with camera matrices.
%xy     - 2-by-N-by-K array with measured coordinates.
%OP     - 3-by-K array with estimated OP coordinates.
%r      - 1-by-K vector with averaged squared residual.


n=size(P,3);
k=size(xy,3);

% Each camera center is a point on the ray in space.
C=zeros(3,n);
for i=1:n
    C(:,i)=euclidean(null(P(:,:,i)));
end

% Expand pts to homogeneous.
xy1=homogeneous(xy);

% Get a second point on each ray and ray direction vector.
Ppx=zeros(4,k,n);
for i=1:n
    % Pre-multiply with psuedo-inverse to get a second point on the ray
    % in space (since P*(P^+*xy)=xy).
    Ppx(:,:,i)=pinv(P(:,:,i))*reshape(xy1(:,i,:),3,k);
    % Move any points "close" to infinity closer.
    far=abs(Ppx(end,:,i))<1e-8;
    % Replace with closer points.
    Ppx(:,far,i)=Ppx(:,far,i)+repmat(homogeneous(C(:,i)),1,nnz(far));
end
    
% Direction vector for each ray.
t=euclidean(Ppx)-repmat(reshape(C,3,1,n),1,k);
% Normalize.
t=t./repmat(sqrt(sum(t.^2,1)),[3,1,1]);

OP=nan(3,k);
OPc=nan(1,k);

% Ray equation:
% Ci - ti*alphai = p
%
% Stack to obtain (for two rays).
%
% [ I t1   0] * [p         [C1
% [ I  0  t2]    alpha1  =  C2]
%                alpha2] 
%
%      A            x    =  b

% For each OP...

% Set up A. First 3 columns will remain the same.
A=[repmat(eye(3),n,1),zeros(3*n,n)];
% Subscripts of elements of last n columns of A that will change for each
% point.
ind=sub2ind(size(A),1:3*n,reshape(repmat(3+(1:n),3,1),1,[]));
% Right-hand-side will be the same.
b=reshape(C,3*n,1);
for j=1:k
    A(ind)=t(:,j,:);
    % Solve for optimal point.
    x=A\b;
    % Optimal point.
    OP(:,j)=x(1:3);
    if (nargout>1)
        % Residual for this point.
        r(j)=norm(b-A*x)/n;
    end
end
