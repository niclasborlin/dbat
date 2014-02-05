function a=angles(s)
%ANGLES Object point angles for a project.
%
%   A=ANGLES(S), where S is a struct returned by PROB2DBATSTRUCT with N
%   object points, returns an N-vector A with the maximum angle in radians
%   between rays for each object point. The maximum angle is the angle
%   closest to being orthogonal between pairs of rays for each object point.

% $Id$

a=nan(size(s.OPid));

for i=1:length(s.OPid)
    % Point
    p=s.OP(:,i);
    % Camera centers.
    cc=s.EO(1:3,s.vis(i,:));
    % Direction vectors.
    d=repmat(p,1,size(cc,2))-cc;
    % Normalize
    dn=d./repmat(sqrt(sum(d.^2,1)),3,1);
    % Compute all inner products. Guard for round-off errors.
    ip=max(min(dn'*dn,1),-1);
    % Angle is acos of inner product of normalized vectors.
    a(i)=max(acos(abs(ip(:))));
end
