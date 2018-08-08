function [tri,area,T,A]=largesttriangle(pts,cHull)
%LARGESTTRIANGLE Select largest triangle(s) from a point set.
%
%[tri,area,T,A]=largesttriangle(pts[,cHull])
%pts   - 2xn or 3xn matrix of points.
%cHull - should we only include points in the convex hull? Default: true.
%        Only affects T and A return parameters.
%tri   - indices of points of largest triangle pts(tri,:).
%area  - area of pts(tri,:).
%T     - mx3 index matrix with triangles sorted descending by area.
%A     - m-vector with areas corresponding to T.
%
%Uses NCHOOSEK, so large point sets might produce memory problems.


if (nargin<2), cHull=true; end

if (cHull)
    % Only choose among points on the convex hull.
    i=convhulln(pts');
    i=unique(i(:));
else
    % Choose among all points.
    i=(1:size(pts,2))';
end

% Get all possible triplets.
T=nchoosek(i,3);

% Calculate area for each triplet.
A=zeros(size(T,1),1);

for j=1:size(T,1)
    A(j)=polyarea(pts(1,T(j,:)),pts(2,T(j,:)));
end

[~,i]=sort(-A);
T=T(i,:);
A=A(i);
tri=T(1,:);
area=A(1);
