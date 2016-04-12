function s=loadpsz2(psFile)
%LOADPSZ2 Load psz file

psDir=fileparts(psFile);
unpackpsz(psFile);
psUnpackedDir=fullfile(psDir,'unpacked');
fName=fullfile(psUnpackedDir,'doc.xml');
s=xml2struct2(fName);

ptCloud=s.document.chunks.chunk.frames.frame.point_cloud;

[~,~,points,~]=ply_read(fullfile(psUnpackedDir,ptCloud.points.Attributes.path),'tri');
s.points=points;
s.objPts=[points.vertex.id,points.vertex.x,points.vertex.y,points.vertex.z];

[~,~,tracks,~]=ply_read(fullfile(psUnpackedDir,ptCloud.tracks.Attributes.path),'tri');
s.tracks=tracks;

projections=cell(size(ptCloud.projections));

cameraIds=cellfun(@(x)sscanf(x.Attributes.camera_id,'%d')+1,ptCloud.projections);
for i=1:length(projections)
    j=find(cameraIds==i);
    [~,~,proj,~]=ply_read(fullfile(psUnpackedDir, ...
                                   ptCloud.projections{j}.Attributes.path),...
                          'tri');
    projections{j}=proj;
end
s.projections=projections;

markPts=[];

for i=1:length(projections)
    ni=length(projections{i}.vertex.id);
    markPts=[markPts;repmat(i-1,ni,1),projections{i}.vertex.id,...
             projections{i}.vertex.x,projections{i}.vertex.y];
end
s.markPts=markPts;

camera=s.document.chunks.chunk.frames.frame.cameras.camera;
cameraIds=cellfun(@(x)sscanf(x.Attributes.camera_id,'%d')+1,camera);

imNames=cell(1,length(camera));
for i=1:length(camera)
    j=find(cameraIds==i);
    imNames{j}=fullfile(psDir,camera{j}.photo.Attributes.path);
end
s.imNames=imNames;