function remove2objpts(pszFile)
%REMOVE2OBJPTS Remove 2-ray obj pts from PSZ PLY files.
%
%    REMOVE2OBJPTS(PSZFILE) unpacks the .psz file PSZFILE and
%    removes the 2-ray object points from the .ply files in the
%    unpacked dir.

psz=loadpsz(pszFile,true);

[prob,pmReport,pts3d,pts2d]=ps2pmstruct(psz);

s0=prob2dbatstruct(prob);

allIds=[psz.raw.ctrlPts(:,1);psz.raw.objPts(:,1)];

obj3rays=sum(s0.vis,2)>=3 & ~s0.isCtrl;

ids2keep=allIds(obj3rays);

fprintf('Keeping %d pts...\n',length(ids2keep));

[~,~,points,~]=ply_read(psz.raw.paths.points,'tri');

keep=ismember(points.vertex.id,ids2keep);

points.vertex.x=points.vertex.x(keep);
points.vertex.y=points.vertex.y(keep);
points.vertex.z=points.vertex.z(keep);
points.vertex.id=points.vertex.id(keep);

ply_write(points,psz.raw.paths.points,'binary_little_endian')

for i=1:length(psz.raw.paths.projections)
    [~,~,pts,~]=ply_read(psz.raw.paths.projections{i},'tri');

    keep=ismember(pts.vertex.id,ids2keep);

    pts.vertex.x=pts.vertex.x(keep);
    pts.vertex.y=pts.vertex.y(keep);
    pts.vertex.size=pts.vertex.size(keep);
    pts.vertex.id=pts.vertex.id(keep);

    ply_write(pts,psz.raw.paths.projections{i},'binary_little_endian')
end
