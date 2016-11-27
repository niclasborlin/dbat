function d=ptdepth(P,X)
%PTDEPTH Determine depth of a point w.r.t. a camera.
%
%d=ptdepth(P,X)
%P - 3x4 camera matrix.
%X - 3xn (Euclidean) or 4xn (homogeneous) matrix with object coordinates.
%d - n-vector with depth of each point.


if size(X,1)<4, X=homogeneous(X); end

x=P*X;
M=P(:,1:3);
d=sign(det(M))*x(3,:)./X(4,:)/norm(M(3,:));
