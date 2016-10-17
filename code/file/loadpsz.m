function s=loadpsz(psFile,varargin)
%LOADPSZ Load Photoscan .PSZ file.
%
%   S=LOADPSZ(FILE) loads the PhotoScan .PSZ file in FILE into a struct S.
%   If the .PSZ file has multiple chunks, chunk #1 is processed.
%
%   S=LOADPSZ(FILE,N) loads chunk N instead of chunk #1.
%
%   The struct S has fields
%   document  - recursive struct that correspond to the XML
%               structure in FILE.
%   transform - struct with fields
%               R, T, S - 4x4 matrices with the individual rotation,
%                         translation, and scaling transformations,
%                         respectively, 
%               L2G     - 4x4 matrix with the composite local-to-global
%                         transformation,
%               G2L     - 4x4 matrix with the composite global-to-local
%                         transformation.
%   raw       - struct with all 3D information in raw coordinates,
%   global,
%   local     - struct with 3D info in global/local coordinates with fields
%               ctrlPts - MC-by-7 array with [id,x,y,z,sx,sy,sz] for ctrl pts,
%               objPts  - MO-by-4 array with [id,x,y,z] for object pts,
%               P       - 3-by-4-by-N array with camera matrices,
%               CC      - 3-by-N array with camera centers,
%               R       - 3-by-3-by-N array with camera rotation matrices.
%   cameraIds - N vector with camera ids,
%   imNames   - N-vector with image file names,
%   markPts   - struct with fields
%               obj  - MMO-by-4 array with [imNo,id,x,y] for object points,
%               ctrl - MMC-by-4 array with [imNo,id,x,y] for ctrl points,
%               all  - concatenation of obj and ctrl.
%   K         - 3-by-3 camera calibration matrix,
%   camera    - struct with camera information with fields
%               imSz         - 2-vector with image size [w,h] in pixels,
%               sensorFormat - 2-vector with sensor size [w,h] in mm,
%               pixelSz      - 2-vector with pixel size [pw,ph] in mm,
%               focal        - scalar with focal length in mm,
%               pp           - 2-vector with principal point in mm.
%
%   By default, LOADPSZ unpacks the .PSZ file (a .ZIP archive) into a
%   directory in TEMPDIR and deletes the unpacked files after loading.
%   LOADPSZ(FILE,...,TRUE) instead unpacks the files into a FILE-local
%   subdir called 'unpacked' and does not delete the unpacked files.
%   Additionally, LOADPSZ(FILE,TRUE,TRUE) creates ascii versions of
%   each .PLY file in a further 'ascii' subdir.
%
%See also: TEMPDIR, TEMPNAME, UNPACKPSZ.

% Default values.
chunkNo=1;
unpackLocal=false;
asciiToo=false;

% Deal with numeric chunk number.
if ~isempty(varargin) && isnumeric(varargin{1})
    chunkNo=varargin{1};
    varargin(1)=[];
end

if length(varargin)>=1
    unpackLocal=varargin{1};
end

if length(varargin)>=2
    asciiToo=varargin{2};
end

psDir=fileparts(psFile);
if unpackLocal
    % Unpack to 'unpacked' subdir.
    unpackDir=fullfile(psDir,'unpacked');
else
    % Unpack to temporary dir.
    unpackDir=tempname;
end

% Unpack the .psz file.
dirs=unpackpsz(psFile,unpackDir,asciiToo);
% Load project data from the main xml file.
fName=fullfile(unpackDir,'doc.xml');
s=dbatxml2struct(fName);

% Chunks can be rearranged in PS so just get the first one.
if iscell(s.document.chunks.chunk)
    if length(s.document.chunks.chunk)>=chunkNo
        chnk = s.document.chunks.chunk{chunkNo};
    else
        error('LOADPSZ: Chunk number out of bounds.');
    end
else
    if chunkNo==1
        chnk = s.document.chunks.chunk
    else
        error('LOADPSZ: Chunk number out of bounds.');
    end        
end
% Extract local-to-global transformation.
if ~isfield(chnk,'transform')
    warning('No local-to-global transform. Using defaults.');
    R=eye(4);
    T=eye(4);
    S=eye(4);
else
    % Use actual transform.
    xform=chnk.transform;
    R=blkdiag(reshape(sscanf(xform.rotation.Text,'%g '),3,3)',1);
    T=[eye(3),sscanf(xform.translation.Text,'%g ');0,0,0,1];
    S=diag([repmat(sscanf(xform.scale.Text,'%g '),1,3),1]);
end

% Transformations between Local and global coordinate systems.
L2G=T*S*R;
G2L=R'*inv(S)*inv(T);

s.transform.R=R;
s.transform.T=T;
s.transform.S=S;
s.transform.L2G=L2G;
s.transform.G2L=G2L;

ptCloud=chnk.frames.frame.point_cloud;

if unpackLocal
    s.raw.paths.points=fullfile(unpackDir,ptCloud.points.Attributes.path);
else
    s.raw.paths.points='';
end

% Object points are in local coordinates.
[~,~,points,~]=ply_read(fullfile(unpackDir,ptCloud.points.Attributes.path),'tri');
s.raw.points=points;
s.raw.objPts=[points.vertex.id,points.vertex.x,points.vertex.y,points.vertex.z];

% Ctrl points are in global coordinates.
if isfield(chnk,'markers')
    markers=chnk.markers.marker;
    if ~iscell(markers)
        markers={markers};
    end
else
    markers=cell(0);
end

ctrlPts=nan(length(markers),7);

% Get default marker std.
settingsProps=chnk.settings.property;
if ~iscell(settingsProps), settingsProps={settingsProps}; end
settingsPropNames=cellfun(@(x)x.Attributes.name,settingsProps,...
                          'uniformoutput',false);
accMarkersIx=find(strcmp(settingsPropNames,'accuracy_markers'));
if length(accMarkersIx)==1
    accMarkers=sscanf(settingsProps{accMarkersIx}.Attributes.value,'%g');
else
    accMarkers=nan;
end

for i=1:size(ctrlPts,1);
    m=markers{i};
    id=sscanf(m.Attributes.id,'%d');
    x=nan;
    y=nan;
    z=nan;
    % Use default marker std setting.
    sx=accMarkers;
    sy=accMarkers;
    sz=accMarkers;
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
    objIdShift=max(s.local.ctrlPts(:,1))+1;
end

s.local.objPts=s.raw.objPts;
s.local.objPts(:,1)=s.local.objPts(:,1)+objIdShift;
s.global.objPts=XformPtsi(s.local.objPts,L2G);

if unpackLocal
    s.raw.paths.tracks=fullfile(unpackDir,ptCloud.tracks.Attributes.path);
else
    s.raw.paths.tracks='';
end

[~,~,tracks,~]=ply_read(fullfile(unpackDir,ptCloud.tracks.Attributes.path),'tri');
s.raw.tracks=tracks;

% Image coordinates.
projections=cell(size(ptCloud.projections));
s.raw.paths.projections=cell(size(ptCloud.projections));

cameraIds=cellfun(@(x)sscanf(x.Attributes.camera_id,'%d')+1, ...
                  ptCloud.projections);

s.cameraIds=cameraIds;

for i=1:length(projections)
    if unpackLocal
        s.raw.paths.projections{i}=...
            fullfile(unpackDir,ptCloud.projections{i}.Attributes.path);
    else
        s.raw.paths.points='';
    end
    [~,~,proj,~]=ply_read(fullfile(unpackDir, ...
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

if isfield(chnk.frames.frame,'markers')
    marker=chnk.frames.frame.markers.marker;
    if ~iscell(marker)
        marker={marker};
    end
else
    marker=cell(1,0);
end

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
if ~isempty(s.markPts.ctrl)
    s.markPts.ctrl(:,2)=s.markPts.ctrl(:,2)+ctrlIdShift;
end

s.markPts.all=msort([s.markPts.ctrl;s.markPts.obj]);

% Transformations are from "image" coordinate system to local.
camera=chnk.cameras.camera;
if ~iscell(camera)
    camera={camera};
end
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

camera=chnk.frames.frame.cameras.camera;
if ~iscell(camera)
    camera={camera};
end
cameraIds=cellfun(@(x)sscanf(x.Attributes.camera_id,'%d')+1,camera);

imNames=cell(1,length(camera));
for i=1:length(camera)
    imNames{i}=fullfile(psDir,camera{i}.photo.Attributes.path);
end
s.imNames=imNames;

cal=chnk.sensors.sensor.calibration;
fx=sscanf(cal.fx.Text,'%g');
fy=sscanf(cal.fy.Text,'%g');
cx=sscanf(cal.cx.Text,'%g');
cy=sscanf(cal.cy.Text,'%g');
K=[fx,0,cx;0,fy,cy;0,0,1];

s.K=K;

sensor=chnk.sensors.sensor;
imSz=[sscanf(sensor.resolution.Attributes.width,'%d'),...
      sscanf(sensor.resolution.Attributes.height,'%d')];

sProps=sensor.property;
if ~iscell(sProps), sProps={sProps}; end

sensorProps=cellfun(@(x)x.Attributes.name,sProps,'uniformoutput',false);
pwProp=sProps(strcmp(sensorProps,'pixel_width'));
phProp=sProps(strcmp(sensorProps,'pixel_height'));
if isempty(pwProp) || isempty(phProp)
    warning('No pixel size specified, using unity.');
    pixelWidth=1;
    pixelHeight=1;
else
    pixelWidth=sscanf(pwProp{1}.Attributes.value,'%g');
    pixelHeight=sscanf(phProp{1}.Attributes.value,'%g');
end

s.camera.imSz=imSz;
s.camera.pixelSz=[pixelWidth,pixelHeight];
s.camera.sensorFormat=s.camera.imSz.*s.camera.pixelSz;
s.camera.focal=fx*s.camera.pixelSz(1);
s.camera.pp=[cx,cy].*s.camera.pixelSz;

% Delete unpacked files unless they should be kept.
if ~unpackLocal
    for i=length(dirs):-1:1
        delete(fullfile(dirs{i},'*'));
        rmdir(dirs{i});
    end
end


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
