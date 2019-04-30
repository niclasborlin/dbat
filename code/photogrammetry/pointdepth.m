function d=pointdepth(s,cams,pts)
%POINTDEPTH Compute point depth with respect to cameras from DBAT struct.
%
%   D=POINTDEPTH(S,CAMS,PTS) computes the point depth for each point
%   with index PTS with respect to each camera with index CAMS. The
%   string 'all' can be used to mean all points and/or cameras. The
%   rows of the returned array D correspond to points, the column
%   correspond to cameras. A value of D(I,J)=NaN is returned for a
%   point I that has not been observed in camera J.
%
%   D=POINTDEPTH(S,CAMS) and D=POINTDEPTH(S) uses PTS='all' and
%   CAMS='all'.

% Extract visibility matrix.
vis=s.IP.vis;

if nargin<2 || (ischar(cams) && strcmp(cams,'all'))
    cams=1:size(vis,2);
end

if nargin<3 || (ischar(pts) && strcmp(pts,'all'))
    pts=1:size(vis,1);
end

% Create camera matrices.
PP=zeros(3,4,length(cams));
for i=1:length(cams)
    % Camera calibration matrix.
    K=getcamvals(s,s.EO.cam(cams(i)),'KCAM');
    % Pose matrix.
    P=geteovals(s,cams(i),'P');
    % Full camera matrix.
    PP(:,:,i)=K*P;
end

% Extract points/cameras of interest.
vis=vis(pts,cams);

d=nan(size(vis));

% For each camera with a visible point.
for i=find(any(vis,1))
    j=vis(:,i);
    d(j,i)=-ptdepth(PP(:,:,i),s.OP.val(:,pts(j)));
end
