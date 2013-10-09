function [P1,P2,X,inFront,sp]=camsfrome2(E,x,xp,frontDir)
%CAMSFROME Calculate cameras from an essential matrix.
%
%[P1,P2,X,inFront,sp]=camsfrome(E,x,xp,frontDir)
%E        - 3x3 essential matrix.
%x        - 3xn matrix with normalized 2d points in image 1.
%xp       - 3xn matrix with normalized 2d points in image 2.
%frontDir - Which Z direction +1/-1 is in front of the camera?
%P1       - Canonic camera [I,0].
%P2       - Other camera [R,t] with baseline length |t|=1.
%X        - 4xn matrix of reconstructed 3d points with respect to P1.
%inFront  - logical n-vector specifying if point X(:,i) is in front of
%           both cameras.
%sp       - sorted vector with number of points in front of each camera config.

[U,S,V]=svd(E);
if (det(U)<0)
	U(:,3)=-U(:,3);
end
if (det(V)<0)
	V(:,3)=-V(:,3);
end

W=[0,-1,0;1,0,0;0,0,1];
Z=[0,1,0;-1,0,0;0,0,0];

% Camera matrices. Second camera will change during iterations.
% Canonical camera.
P1=eye(3,4);
P=cat(3,P1,zeros(3,4));

% Four choices for camera 2
C2=cat(3,[U*W*V',U(:,3)],[U*W*V',-U(:,3)],[U*W'*V',U(:,3)],[U*W'*V',-U(:,3)]);
% Reconstructed points for each choice.
XX=cell(1,size(C2,3));
% Is the z coordinate positive?
zp=zeros(size(C2,3),size(x,2));
% Convert points to Euclidean and stack for pm_forwintersect3.
xE=euclidean(x);
xpE=euclidean(xp);
xy=cat(2,reshape(xE,2,1,[]),reshape(xpE,2,1,[]));
for j=1:length(C2)
	P(:,:,2)=C2(:,:,j);
	% Reconstruct points.
    X=pm_forwintersect3(P,xy);
	% Check what points are in front of both cameras.
	d1=ptdepth(P(:,:,1),X);
	d2=ptdepth(P(:,:,2),X);
	zp(j,:)=sign(d1)==sign(frontDir) & sign(d2)==sign(frontDir);
    XX{j}=X;
end

% Select the camera pair that has the largest number of points in front
% of both cameras.
[sp,ii]=sort(-sum(zp,2));
sp=-sp;
P2=C2(:,:,ii(1));
inFront=zp(ii(1),:);
X=XX{ii(1)};
