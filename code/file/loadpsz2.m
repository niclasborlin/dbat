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

ctrlPts=nan(length(s.document.chunks.chunk.markers.marker),4);
for i=1:size(ctrlPts,1);
    m=s.document.chunks.chunk.markers.marker{i};
    id=sscanf(m.Attributes.id,'%d');
    x=sscanf(m.reference.Attributes.x,'%g');
    y=sscanf(m.reference.Attributes.y,'%g');
    z=sscanf(m.reference.Attributes.z,'%g');
    ctrlPts(i,:)=[id,x,y,z];
end

s.ctrlPts=ctrlPts;

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

camera=s.document.chunks.chunk.cameras.camera;
cameraIds=cellfun(@(x)sscanf(x.Attributes.id,'%d')+1,camera);
xforms=nan(4,4,length(cameraIds));
for i=1:length(cameraIds)
    j=find(cameraIds==i);
    T=reshape(sscanf(camera{j}.transform.Text,'%g '),4,4)';
    xforms(:,:,i)=T;
end
s.transforms=xforms;

xform=s.document.chunks.chunk.transform;
R=blkdiag(reshape(sscanf(xform.rotation.Text,'%g '),3,3)',1);
T=[eye(3),sscanf(xform.translation.Text,'%g ');0,0,0,1];
S=diag([repmat(sscanf(xform.scale.Text,'%g '),1,3),1]);

s.transform.R0=R;
s.transform.T0=T;
s.transform.S0=S;
s.transform.TSR=T*S*R;

camera=s.document.chunks.chunk.frames.frame.cameras.camera;
cameraIds=cellfun(@(x)sscanf(x.Attributes.camera_id,'%d')+1,camera);

imNames=cell(1,length(camera));
for i=1:length(camera)
    j=find(cameraIds==i);
    imNames{j}=fullfile(psDir,camera{j}.photo.Attributes.path);
end
s.imNames=imNames;

cal=s.document.chunks.chunk.sensors.sensor.calibration;
fx=sscanf(cal.fx.Text,'%g');
fy=sscanf(cal.fy.Text,'%g');
cx=sscanf(cal.cx.Text,'%g');
cy=sscanf(cal.cy.Text,'%g');
K=[fx,0,cx;0,fy,cy;0,0,1];

s.K=K;

s.ctrlPtsLocal=[s.ctrlPts(:,1),euclidean(inv(s.transform.TSR)*homogeneous(s.ctrlPts(:,2:4)'))'];

s.objPtsGlobal=[s.objPts(:,1),euclidean(s.transform.TSR*homogeneous(s.objPts(:,2:4)'))'];

globalXforms=nan(size(xforms));
for i=1:size(globalXforms,3)
    globalXforms(:,:,i)=s.transform.TSR*xforms(:,:,i);

end

s.transformsGlobal=globalXforms;