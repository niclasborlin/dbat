function s=loadpsz2(psFile)
%LOADPSZ2 Load psz file

psDir=fileparts(psFile);
unpackpsz(psFile);
psUnpackedDir=fullfile(psDir,'unpacked');
fName=fullfile(psUnpackedDir,'doc.xml');
s=xml2struct2(fName);

% Extract local-to-global transformation.
xform=s.document.chunks.chunk.transform;
R=blkdiag(reshape(sscanf(xform.rotation.Text,'%g '),3,3)',1);
T=[eye(3),sscanf(xform.translation.Text,'%g ');0,0,0,1];
S=diag([repmat(sscanf(xform.scale.Text,'%g '),1,3),1]);

% Transformations between Local and global coordinate systems.
L2G=T*S*R;
G2L=R'*inv(S)*inv(T);

s.transform.R0=R;
s.transform.T0=T;
s.transform.S0=S;
s.transform.L2G=L2G;
s.transform.G2L=G2L;

ptCloud=s.document.chunks.chunk.frames.frame.point_cloud;

% Object points are in local coordinates.
[~,~,points,~]=ply_read(fullfile(psUnpackedDir,ptCloud.points.Attributes.path),'tri');
s.raw.points=points;
s.raw.objPts=[points.vertex.id,points.vertex.x,points.vertex.y,points.vertex.z];

% Ctrl points are in global coordinates.
ctrlPts=nan(length(s.document.chunks.chunk.markers.marker),7);
for i=1:size(ctrlPts,1);
    m=s.document.chunks.chunk.markers.marker{i};
    id=sscanf(m.Attributes.id,'%d');
    x=nan;
    y=nan;
    z=nan;
    sx=nan;
    sy=nan;
    sz=nan;
    if isfield(m,'reference')
        if isfield(m.reference.Attributes,'x')
            x=sscanf(m.reference.Attributes.x,'%g');
        end
        if isfield(m.reference.Attributes,'y')
            y=sscanf(m.reference.Attributes.y,'%g');
        end
        if isfield(m.reference.Attributes,'z')
            z=sscanf(m.reference.Attributes.z,'%g');
        end
        if isfield(m.reference.Attributes,'sxyz')
            sx=sscanf(m.reference.Attributes.sxyz,'%g');
            sy=sx;
            sz=sx;
        elseif isfield(m.reference.Attributes,'sxy')
            sx=sscanf(m.reference.Attributes.sxy,'%g');
            sy=sx;
        end
        if isfield(m.reference.Attributes,'sx')
            sx=sscanf(m.reference.Attributes.sx,'%g');
        end    
        if isfield(m.reference.Attributes,'sy')
            sy=sscanf(m.reference.Attributes.sy,'%g');
        end    
        if isfield(m.reference.Attributes,'sz')
            sz=sscanf(m.reference.Attributes.sz,'%g');
        end    
        ctrlPts(i,:)=[id,x,y,z,sx,sy,sz];
    end
end
ctrlPts(any(isnan(ctrlPts(:,1:4)),2),:)=[];
s.raw.ctrlPts=ctrlPts;

% Make local/global ctrl pt ids 1-based.
if isempty(ctrlPts)
    ctrlIdShift=0;
else
    ctrlIdShift=1-min(ctrlPts(:,1));
end

s.global.ctrlPts=s.raw.ctrlPts;
s.global.ctrlPts(:,1)=s.global.ctrlPts(:,1)+ctrlIdShift;
s.local.ctrlPts=XformPtsi(s.global.ctrlPts,G2L);

% Shift object point ids to above control point ids.
if isempty(s.raw.objPts) || isempty(s.raw.ctrlPts);
    objIdShift=0;
else
    objIdShift=max(s.local.ctrlPts(:,1))+1-min(s.raw.objPts(:,1));
end

s.local.objPts=s.raw.objPts;
s.local.objPts(:,1)=s.local.objPts(:,1)+objIdShift;
s.global.objPts=XformPtsi(s.local.objPts,L2G);

[~,~,tracks,~]=ply_read(fullfile(psUnpackedDir,ptCloud.tracks.Attributes.path),'tri');
s.raw.tracks=tracks;

% Image coordinates.
projections=cell(size(ptCloud.projections));

cameraIds=cellfun(@(x)sscanf(x.Attributes.camera_id,'%d')+1, ...
                  ptCloud.projections);

s.cameraIds=cameraIds;

for i=1:length(projections)
    [~,~,proj,~]=ply_read(fullfile(psUnpackedDir, ...
                                   ptCloud.projections{i}.Attributes.path),...
                          'tri');
    projections{i}=proj;
end
s.raw.projections=projections;

objMarkPts=[];

for i=1:length(projections)
    ni=length(projections{i}.vertex.id);
    objMarkPts=[objMarkPts;repmat(i-1,ni,1),projections{i}.vertex.id,...
             projections{i}.vertex.x,projections{i}.vertex.y];
end
objMarkPts=msort(objMarkPts);
s.raw.objMarkPts=objMarkPts;
s.markPts.obj=objMarkPts;
s.markPts.obj(:,2)=s.markPts.obj(:,2)+objIdShift;

marker=s.document.chunks.chunk.frames.frame.markers.marker;

ids=cellfun(@(x)sscanf(x.Attributes.marker_id,'%d'),marker);

ctrlMarkPts=[];

for i=1:length(marker)
    if isfield(marker{i},'location')
        location=marker{i}.location;
        if ~iscell(location)
            location={location};
        end
        camIds=cellfun(@(x)sscanf(x.Attributes.camera_id,'%d'),location);
        x=cellfun(@(m)sscanf(m.Attributes.x,'%g'),location);
        y=cellfun(@(m)sscanf(m.Attributes.y,'%g'),location);
        ctrlMarkPts=[ctrlMarkPts;[camIds;repmat(ids(i),size(camIds));x;y]'];
    end
end
s.raw.ctrlMarkPts=msort(ctrlMarkPts);
s.markPts.ctrl=s.raw.ctrlMarkPts;
s.markPts.ctrl(:,2)=s.markPts.ctrl(:,2)+ctrlIdShift;

s.markPts.all=msort([s.markPts.ctrl;s.markPts.obj]);

% Transformations are from "image" coordinate system to local.
camera=s.document.chunks.chunk.cameras.camera;
cameraIds=cellfun(@(x)sscanf(x.Attributes.id,'%d')+1,camera);
xforms=nan(4,4,length(cameraIds));
% Camera matrices from local coordinates.
P=nan(3,4,length(cameraIds));
% Camera centers in local coordinates.
CC=nan(3,length(cameraIds));
for i=1:length(cameraIds)
    T=reshape(sscanf(camera{i}.transform.Text,'%g '),4,4)';
    xforms(:,:,i)=T;
    P(:,:,i)=eye(3,4)*inv(T*diag([1,-1,-1,1]));
    CC(:,i)=euclidean(null(P(:,:,i)));
end
s.raw.transforms=xforms;
s.raw.P=P;
s.raw.CC=CC;

s.local.P=s.raw.P;
s.local.CC=s.raw.CC;
s.local.R=nan(3,3,size(s.local.P,3));
for i=1:size(s.local.R,3)
    R=s.local.P(:,1:3,i);
    s.local.R(:,:,i)=R/det(R)^(1/3);
end

s.global.P=XformCams(s.local.P,L2G);
s.global.CC=XformPts(s.local.CC,L2G);
s.global.R=nan(3,3,size(s.global.P,3));
for i=1:size(s.global.R,3)
    R=s.global.P(:,1:3,i);
    s.global.R(:,:,i)=R/det(R)^(1/3);
end

camera=s.document.chunks.chunk.frames.frame.cameras.camera;
cameraIds=cellfun(@(x)sscanf(x.Attributes.camera_id,'%d')+1,camera);

imNames=cell(1,length(camera));
for i=1:length(camera)
    imNames{i}=fullfile(psDir,camera{i}.photo.Attributes.path);
end
s.imNames=imNames;

cal=s.document.chunks.chunk.sensors.sensor.calibration;
fx=sscanf(cal.fx.Text,'%g');
fy=sscanf(cal.fy.Text,'%g');
cx=sscanf(cal.cx.Text,'%g');
cy=sscanf(cal.cy.Text,'%g');
K=[fx,0,cx;0,fy,cy;0,0,1];

s.K=K;

sensor=s.document.chunks.chunk.sensors.sensor;
imSz=[sscanf(sensor.resolution.Attributes.width,'%d'),...
      sscanf(sensor.resolution.Attributes.height,'%d')];

sensorProps=cellfun(@(x)x.Attributes.name,sensor.property,'uniformoutput',false);
pixelWidth=sscanf(sensor.property{strcmp(sensorProps,'pixel_width')}.Attributes.value,'%g');
pixelHeight=sscanf(sensor.property{strcmp(sensorProps,'pixel_height')}.Attributes.value,'%g');

s.camera.imSz=imSz;
s.camera.pixelSz=[pixelWidth,pixelHeight];
s.camera.sensorFormat=s.camera.imSz.*s.camera.pixelSz;
s.camera.focal=fx*s.camera.pixelSz(1);
s.camera.pp=[cx,cy].*s.camera.pixelSz;

function q=XformPts(p,M,forceHomogeneous)
%Apply 4-by-4 point transformation matrix M to 3D euclidean or
%homogeneous points. Output has same format as input. Use
%XformPts(p,M,true) to force homogeneous output.

if nargin<3, forceHomogeneous=false; end

if size(p,1)==3
    q=M*homogeneous(p);
    if ~forceHomogeneous
        q=euclidean(q);
    end
else
    q=M*p;
end

function q=XformPtsi(p,M)
%Apply 4-by-4 point transformation matrix M to points with ids p=[id,x,y,z].

q=[p(:,1),XformPts(p(:,2:4)',M)'];


function Q=XformCams(P,M)
%Apply 4-by-4 point transformation matrix M to 3-by-4-by-K array P
%of camera matrices.

Q=nan(size(P));
for i=1:size(P,3)
    Q(:,:,i)=P(:,:,i)*inv(M);
end

