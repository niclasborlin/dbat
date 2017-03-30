function s=loadpsz(psFile,varargin)
%LOADPSZ Load Photoscan .PSZ file.
%
%   S=LOADPSZ(FILE) loads the PhotoScan .PSZ file in FILE into a struct S.
%   If the .PSZ file has multiple chunks, chunk #1 is processed.
%
%   S=LOADPSZ(FILE,N) loads chunk N instead of chunk #1.
%
%   The struct S has fields
%   fileName  - string with file name.
%   document  - recursive struct that correspond to the XML
%               structure in FILE.
%   transform - struct with fields
%               R, T, S - 4x4 matrices with the individual rotation,
%                         translation, and scaling transformations,
%                         respectively, 
%               G2L     - 4x4 matrix with the composite global-to-local
%                         transformation.
%               L2G     - 4x4 matrix with the composite local-to-global
%                         transformation,
%               G2SL    - 4x4 matrix with the composite global-to-semilocal
%                         (no rotation) transformation.
%               SL2G    - 4x4 matrix with the composite semilocal-to-global
%                         transformation,
%   raw       - struct with all 3D information in raw coordinates,
%   global,
%   local,
%   semilocal - structs with 3D info in global/local coordinates with fields
%               ctrlPts - MC-by-7 array with [id,x,y,z,sx,sy,sz] for ctrl pts,
%               objPts  - MO-by-4 array with [id,x,y,z] for object pts,
%               P       - 3-by-4-by-N array with camera matrices,
%               CC      - 3-by-N array with camera centers,
%               R       - 3-by-3-by-N array with camera rotation matrices.
%   cameraIds - N-vector with camera ids,
%   cameraLabels - N-cell vector with camera labels,
%   cameraEnabled - logical N-vector indicating which cameras are enabled,
%   imNames   - N-vector with image file names,
%   markPts   - struct with fields
%               obj  - MMO-by-4 array with [imNo,id,x,y] for object points,
%               sz   - MO-by-1 array with key point size
%               ctrl - MMC-by-4 array with [imNo,id,x,y] for ctrl points,
%               ctrlPinned - MMC-by-1 logical vector with pinned status,
%               IsPinned - function handle to test if a control
%                          point is pinned.
%               all  - concatenation of obj and ctrl.
%   DBATCamId - function handle to convert from PS to DBAT camera id.
%   PSCamId   - function handle to convert from DBAT to PS camera id.
%   DBATCPid  - function handle to convert from PS to DBAT ctrl pt id.
%   PSCPid    - function handle to convert from DBAT to PS ctrl pt id.
%   DBATOPid  - function handle to convert from PS to DBAT object pt id.
%   PSOPid    - function handle to convert from DBAT to PS object pt id.
%   K         - 3-by-3 camera calibration matrix,
%   vis       - nOP-by-nImages sparse logical visibility matrix,
%   camera    - struct with camera information with fields
%               imSz           - 2-vector with image size [w,h] in pixels,
%               sensorFormat   - 2-vector with sensor size [w,h] in mm,
%               pixelSz        - 2-vector with pixel size [pw,ph] in mm,
%               focal          - scalar with focal length in mm,
%               pp             - 2-vector with principal point in mm.
%               nominalFocal   - nominal focal length in mm.
%               k              - radial lens distortion coefficients in px.
%               p              - tangential lens distortion coefficients in px.
%               isFixed        - true if the camera is fixed.
%               isAdjusted     - true if the camera parameters have been  
%                                adjusted by Photoscan.
%               adjustedParams - struct with fields indicating which
%                                parameters have been adjusted by Photoscan:
%                                f      - scalar
%                                cxcy   - 2-vector
%                                aspect - scalar
%                                skew   - scalar
%                                b      - 2-vector
%                                k      - 4-vector
%                                p      - 4-vector
%   defStd    - struct with default standard deviations
%               tiePoints   - std for automatically detected tie points [pix]
%               projections - std for manually measured markers [pix]
%               markers     - std for marker positions [m]
%               camPos      - std for camera position [m]
%               camAng      - std for camera angles [deg]
%               scaleBars   - std for scale bar lengths [m]
%
%   By default, LOADPSZ unpacks the .PSZ file (a .ZIP archive) into a
%   directory in TEMPDIR and deletes the unpacked files after loading.
%   LOADPSZ(FILE,...,TRUE) instead unpacks the files into a FILE-local
%   subdir called 'unpacked' and does not delete the unpacked files.
%   Additionally, LOADPSZ(FILE,TRUE,TRUE) creates ascii versions of
%   each .PLY file in a further 'ascii' subdir.
%
%   If both an initial and an adjusted camera is available in the .psz
%   file, the adjusted camera is loaded.
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

% Initialize waitbar to delay for 1s and update every 1s.
DelayedWaitBar('init',1,1,'Loading Photoscan project file...');

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
DelayedWaitBar(0.25);

% Load project data from the main xml file.
fName=fullfile(unpackDir,'doc.xml');
s=dbatxml2struct(fName);
DelayedWaitBar(0.3);

s.fileName=psFile;

% Get the requested chunk.
if iscell(s.document.chunks.chunk)
    if length(s.document.chunks.chunk)>=chunkNo
        chnk = s.document.chunks.chunk{chunkNo};
    else
        error('LOADPSZ: Chunk number out of bounds.');
    end
else
    if chunkNo==1
        chnk = s.document.chunks.chunk;
    else
        error('LOADPSZ: Chunk number out of bounds.');
    end        
end

% Extract default standard deviations.
s.defStd=getdefstd(chnk);

% Default transformation matrices.
R=eye(4);
T=eye(4);
S=eye(4);

% Extract local-to-global transformation.
if ~isfield(chnk,'transform')
    warning('No local-to-global transform. Using defaults.');
else
    % Use actual transform.
    xform=chnk.transform;
    if isfield(xform,'rotation')
        R=blkdiag(reshape(sscanf(xform.rotation.Text,'%g '),3,3)',1);
    end
    if isfield(xform,'translation')
        T=[eye(3),sscanf(xform.translation.Text,'%g ');0,0,0,1];
    end
    if isfield(xform,'scale')
        S=diag([repmat(sscanf(xform.scale.Text,'%g '),1,3),1]);
    end
end

% Transformations between local/semilocal and global coordinate systems.
L2G=T*S*R;
SL2G=T*S;
% Avoid explicit inverse for numerical reasons.
G2L=(R'/S)/T; % = R'*inv(S)*inv(T)
G2SL=inv(S)/T;
L2SL=R;
SL2L=R';

s.transform.R=R;
s.transform.T=T;
s.transform.S=S;
s.transform.G2L=G2L;
s.transform.L2G=L2G;
s.transform.G2SL=G2SL;
s.transform.SL2G=SL2G;
s.transform.L2SL=L2SL;
s.transform.SL2L=SL2L;

% Determine what cameras we have.
camera=chnk.cameras.camera;
if ~iscell(camera)
    camera={camera};
end
cameraIds=cellfun(@(x)sscanf(x.Attributes.id,'%d'),camera);
cameraLabels=cellfun(@(x)x.Attributes.label,camera,'uniformoutput',false);
sensorIds=cellfun(@(x)sscanf(x.Attributes.sensor_id,'%d'),camera);
cameraEnabled=cellfun(@(x)strcmp(x.Attributes.enabled,'true'),camera);
if length(unique(sensorIds(cameraEnabled)))>1
    error('Handling of cameras for multiple sensor ids not implemented yet');
end

s.cameraIds=cameraIds;
s.cameraLabels=cameraLabels;
s.cameraEnabled=cameraEnabled;

invCameraIds=nan(max(cameraIds)+1,1);
invCameraIds(cameraIds+1)=1:length(cameraIds);

% Functions to convert between Photoscan camera id and DBAT camera number.
PSCamId=@(id)cameraIds(id);
DBATCamId=@(id)invCameraIds(id+1);
s.PSCamId=PSCamId;
s.DBATCamId=DBATCamId;

% Extract transformation.
% Transformations are from "image" coordinate system to local.
xforms=nan(4,4,length(cameraIds));
% Camera matrices from local coordinates.
P=nan(3,4,length(cameraIds));
% Camera centers in local coordinates.
CC=nan(3,length(cameraIds));
for i=1:length(cameraIds)
    T=reshape(sscanf(camera{i}.transform.Text,'%g '),4,4)';
    xforms(:,:,i)=T;
    if 1
        % TODO: Check this "mirroring"...
        P(:,:,i)=eye(3,4)/(T*diag([1,-1,-1,1]));
    else
        warning('Untested non-mirroring');
        P(:,:,i)=eye(3,4)/T; %#ok<UNRCH> % *inv(T)
    end
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
    if det(R)<0
        warning('Loaded rotation matrix has det(R)!=1')
        disp(det(R))
    end
    s.global.R(:,:,i)=R/det(R)^(1/3);
end

s.semilocal.P=XformCams(s.local.P,L2SL);
s.semilocal.CC=XformPts(s.local.CC,L2SL);
s.semilocal.R=nan(3,3,size(s.semilocal.P,3));
for i=1:size(s.semilocal.R,3)
    R=s.semilocal.P(:,1:3,i);
    if det(R)<0
        warning('Loaded rotation matrix has det(R)!=1')
        disp(det(R))
    end
    s.semilocal.R(:,:,i)=R/det(R)^(1/3);
end


% Process the point cloud.
if isfield(chnk.frames.frame,'point_cloud')
    ptCloud=chnk.frames.frame.point_cloud;
else
    ptCloud=[];
end

if unpackLocal && ~isempty(ptCloud)
    s.raw.paths.points=fullfile(unpackDir,ptCloud.points.Attributes.path);
else
    s.raw.paths.points='';
end

% Object points are in local coordinates.
if ~isempty(ptCloud) && ~isempty(ptCloud.points.Attributes.path)
    [~,~,points,~]=ply_read(fullfile(unpackDir,ptCloud.points.Attributes.path),'tri');
else
    points=[];
end
s.raw.points=points;
   
DelayedWaitBar(0.35);

if ~isempty(points)
    s.raw.objPts=[points.vertex.id,points.vertex.x,points.vertex.y,points.vertex.z];
else
    s.raw.objPts=zeros(0,4);
end

% Process the markers - ctrl points.
if isfield(chnk,'markers')
    markers=chnk.markers.marker;
    if ~iscell(markers)
        markers={markers};
    end
else
    markers=cell(0);
end

% Pre-allocate.
ctrlPts=nan(length(markers),7);
ctrlPtsEnabled=false(length(markers),1);
ctrlPtsLabels=cell(1,length(markers));

for i=1:size(ctrlPts,1);
    m=markers{i};
    id=sscanf(m.Attributes.id,'%d');
    if isfield(m.Attributes,'label')
        ctrlPtsLabels{i}=m.Attributes.label;
    end
    x=nan;
    y=nan;
    z=nan;
    % Use default marker std setting.
    sx=s.defStd.markers;
    sy=s.defStd.markers;
    sz=s.defStd.markers;
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
        if isfield(m.reference.Attributes,'enabled')
            ctrlPtsEnabled(i)=strcmp(m.reference.Attributes.enabled,'true');
        end
    end
    ctrlPts(i,:)=[id,x,y,z,sx,sy,sz];
end
s.raw.ctrlPts=ctrlPts;
s.raw.ctrlPtsLabels=ctrlPtsLabels;
s.raw.ctrlPtsEnabled=ctrlPtsEnabled;

DelayedWaitBar(0.4);

% Make local/global ctrl pt ids 1-based.
rawCPids=ctrlPts(:,1);

invCPids=nan(max(rawCPids)+1,1);
invCPids(rawCPids+1)=1:length(rawCPids);

% Will convert a zero-based id to a one-based id. Generate NaN's for
% out-of-bounds CP ids.
PSCPid=@(id)rawCPids(id);
DBATCPid=@(id)invCPids(id+1);

s.DBATCPid=DBATCPid;
s.PSCPid=PSCPid;

% Copy raw ctrl pts and adjust id.
s.global.ctrlPtsLabels=s.raw.ctrlPtsLabels;
s.global.ctrlPtsEnabled=s.raw.ctrlPtsEnabled;
s.global.ctrlPts=s.raw.ctrlPts;
s.global.ctrlPts(:,1)=DBATCPid(s.global.ctrlPts(:,1));
% Transform ctrl pts from global to local and semilocal coordinate systems.
s.local.ctrlPts=XformPtsi(s.global.ctrlPts,G2L);
s.semilocal.ctrlPts=XformPtsi(s.global.ctrlPts,G2SL,true);

s.semilocal.ctrlPtsLabels=s.raw.ctrlPtsLabels;
s.semilocal.ctrlPtsEnabled=s.raw.ctrlPtsEnabled;

% Highest DBAT CP id.
maxDBATCPid=length(rawCPids);

% Map object point ids to above control point ids. Generate NaN's for
% out-of-bounds OP ids.
DBATOPid=@(id)id+1+maxDBATCPid+0./(id>=0);
PSOPid=@(id)id-1-maxDBATCPid+0./(id>=1+maxDBATCPid);
s.DBATOPid=DBATOPid;
s.PSOPid=PSOPid;

% Copy raw object points and adjust id.
s.local.ctrlPtsLabels=s.raw.ctrlPtsLabels;
s.local.ctrlPtsEnabled=s.raw.ctrlPtsEnabled;
s.local.objPts=s.raw.objPts;
s.local.objPts(:,1)=DBATOPid(s.local.objPts(:,1));
% Transform obj pts from local to global and semilocal coordinate systems.
s.global.objPts=XformPtsi(s.local.objPts,L2G);
s.semilocal.objPts=XformPtsi(s.global.objPts,G2SL);

if unpackLocal && ~isempty(ptCloud) && ~isempty(ptCloud.tracks.Attributes.path)
    s.raw.paths.tracks=fullfile(unpackDir,ptCloud.tracks.Attributes.path);
else
    s.raw.paths.tracks='';
end

if ~isempty(ptCloud) && ~isempty(ptCloud.tracks.Attributes.path)
    [~,~,tracks,~]=ply_read(fullfile(unpackDir,ptCloud.tracks.Attributes.path),'tri');
else
    tracks=[];
end
s.raw.tracks=tracks;

DelayedWaitBar(0.45);

% Load all measured image coordinates.
if ~isempty(ptCloud)
    projs=ptCloud.projections;
else
    projs={};
end
if ~iscell(projs), projs={projs}; end

% Projections will be sorted to match s.cameraId.
projections=cell(size(projs));
s.raw.paths.projections=cell(size(projs));
projCameraIds=cellfun(@(x)sscanf(x.Attributes.camera_id,'%d'),projs);

for i=1:length(projections)
    % Where to store these measured points.
    j=DBATCamId(projCameraIds(i));
    if unpackLocal
        s.raw.paths.projections{j}=fullfile(unpackDir,projs{i}.Attributes.path);
    end
    [~,~,proj,~]=ply_read(fullfile(unpackDir,projs{i}.Attributes.path),'tri');
    projections{j}=proj;
    DelayedWaitBar(0.45+i/length(projections)*0.5);
end
s.raw.projections=projections;

% Process all measured image coordinates.
nPts=cellfun(@(x)length(x.vertex.id),projections);
ptIx=cumsum([0,nPts]);
objMarkPts=nan(sum(nPts),4);
objKeyPtSize=nan(sum(nPts),1);

for i=1:length(projections)
    % Index for where to put the points.
    ni=nPts(i);
    ix=ptIx(i)+1:ptIx(i+1);
    % Store object points with PS ids.
    objMarkPts(ix,:)=[repmat(s.cameraIds(i),ni,1),projections{i}.vertex.id,...
                      projections{i}.vertex.x,projections{i}.vertex.y];
    if isfield(projections{i}.vertex,'size')
        objKeyPtSize(ix)=projections{i}.vertex.size;
    end
end
% Ensure that the mark points are sorted by image, then id.
[objMarkPts,i]=sortrows(objMarkPts,[1,2]);
objKeyPtSize=objKeyPtSize(i);

s.raw.objMarkPts=objMarkPts;
s.raw.objKeyPtSize=objKeyPtSize;

% Create copy with DBAT ids.
s.markPts.obj=objMarkPts;
if ~isempty(s.markPts.obj)
    s.markPts.obj(:,1)=DBATCamId(s.markPts.obj(:,1));
    s.markPts.obj(:,2)=DBATOPid(s.markPts.obj(:,2));
end
% Ensure that the mark points are sorted by image, then id.
[s.markPts.obj,i]=sortrows(s.markPts.obj,[1,2]);
s.markPts.sz=s.raw.objKeyPtSize(i);

% Process measurements of 'markers' - control points.
if isfield(chnk.frames.frame,'markers')
    marker=chnk.frames.frame.markers.marker;
    if ~iscell(marker)
        marker={marker};
    end
else
    marker=cell(1,0);
end

% [camId, markerId, x, y]
ctrlMarkPts=zeros(0,4);
ctrlMarkPtsPinned=false(0,1);

for i=1:length(marker)
    % Extract marker id and convert to marker number.
    markerId=sscanf(marker{i}.Attributes.marker_id,'%d');
    if isfield(marker{i},'location')
        location=marker{i}.location;
        if ~iscell(location)
            location={location};
        end
        % Extract camera ids for each measured point.
        camIds=cellfun(@(x)sscanf(x.Attributes.camera_id,'%d'),location);
        x=cellfun(@(m)sscanf(m.Attributes.x,'%g'),location);
        y=cellfun(@(m)sscanf(m.Attributes.y,'%g'),location);
        pinned=cellfun(@(m)isfield(m.Attributes,'pinned') && ...
                       strcmp(m.Attributes.pinned,'true'),location);
        % What does 'pinned' mean? For now, just warn if a marker measurement
        % is not pinned.
        if ~all(pinned)
            warning('Warning: Unpinned marker measurements!');
        end
        ctrlMarkPts=[ctrlMarkPts;[camIds;repmat(markerId,size(camIds));x;y]']; %#ok<AGROW>
        ctrlMarkPtsPinned=[ctrlMarkPtsPinned;pinned']; %#ok<AGROW>
    end
end
% Sort ctrlMarkPts by marker id, then camera id, to match order in xml file.
[s.raw.ctrlMarkPts,i]=sortrows(ctrlMarkPts,[2,1]);
s.raw.ctrlMarkPtsPinned=ctrlMarkPtsPinned(i);

% Create copy with DBAT ids.
s.markPts.ctrl=s.raw.ctrlMarkPts;
if ~isempty(s.markPts.ctrl)
    s.markPts.ctrl(:,1)=DBATCamId(s.markPts.ctrl(:,1));
    s.markPts.ctrl(:,2)=DBATCPid(s.markPts.ctrl(:,2));
end
% Ensure that the mark points are sorted by image, then id.
[s.markPts.ctrl,i]=sortrows(s.markPts.ctrl,[1,2]);
s.markPts.ctrlPinned=s.raw.ctrlMarkPtsPinned(i);

% Returns true for every [imNo, markerId] that are pinned ctrl pt
% measurements.
s.markPts.IsPinned=@(x)ismember(x(:,1:2),s.markPts.ctrl(s.markPts.ctrlPinned,1:2),'rows');

% Pack all measured points together.
s.markPts.all=sortrows([s.markPts.ctrl;s.markPts.obj],[1,2]);

% Create visibility matrix.
vis=sparse(s.markPts.all(:,2),s.markPts.all(:,1),true);

s.vis=vis;

DelayedWaitBar(1);

% Process image files names. File names will be sorted as s.cameraId.
camera=chnk.frames.frame.cameras.camera;
if ~iscell(camera)
    camera={camera};
end

imNames=cell(1,length(camera));
for i=1:length(camera)
    % Extract camera id.
    camId=sscanf(camera{i}.Attributes.camera_id,'%d');
    % Convert to DBAT camera number.
    j=DBATCamId(camId);
    imNames{j}=fullfile(psDir,camera{i}.photo.Attributes.path);
end
s.imNames=imNames;

% Collect calibrated camera parameters.
cal=chnk.sensors.sensor.calibration;
if iscell(cal)
    % If we have multiple cameras, prefer the adjusted.
    camTypes=cellfun(@(x)x.Attributes.class,cal,'uniformoutput',false);
    adjustedCam=find(strcmp(camTypes,'adjusted'));
    if any(adjustedCam)
        cal=cal{adjustedCam};
    else
        error(['Unknown camera types: ',sprintf(' %s',camTypes{:})]);
    end
end
isAdjusted=strcmp(cal.Attributes.class,'adjusted');
fx=sscanf(cal.fx.Text,'%g');
fy=sscanf(cal.fy.Text,'%g');
cx=sscanf(cal.cx.Text,'%g');
cy=sscanf(cal.cy.Text,'%g');
% Dynamic for lens distortion parameters.
fn=fieldnames(cal);
% Get all 'kN' fields.
kFields=regexp(fn,'^k(\d+)$','tokens');
i=find(~cellfun(@isempty,kFields));
k=zeros(1,0);
for j=1:length(i)
    % Get index.
    n=sscanf(kFields{i(j)}{1}{1},'%d');
    k(n)=sscanf(cal.(fn{i(j)}).Text,'%g');
end
% Get all 'pN' fields.
pFields=regexp(fn,'^p(\d+)$','tokens');
i=find(~cellfun(@isempty,pFields));
p=zeros(1,0);
for j=1:length(i)
    % Get index.
    n=sscanf(pFields{i(j)}{1}{1},'%d');
    p(n)=sscanf(cal.(fn{i(j)}).Text,'%g');
end

K=[fx,0,cx;0,fy,cy;0,0,1];

s.K=K;

sensor=chnk.sensors.sensor;
imSz=[sscanf(sensor.resolution.Attributes.width,'%d'),...
      sscanf(sensor.resolution.Attributes.height,'%d')];

camName=sensor.Attributes.label;
camType=sensor.Attributes.type;

sProps=sensor.property;
if ~iscell(sProps), sProps={sProps}; end

% Specified sensor pixel height/width.
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
% Specified sensor focal length.
fProp=sProps(strcmp(sensorProps,'focal_length'));
nominalFocal=NaN;
if ~isempty(fProp)
    nominalFocal=sscanf(fProp{1}.Attributes.value,'%g');
end
% Is the sensor fixed?
fProp=sProps(strcmp(sensorProps,'fixed'));
sensorFixed=true;
if ~isempty(fProp)
    sensorFixed=strcmp(fProp{1}.Attributes.value,'true');
end

s.camera.name=camName;
s.camera.type=camType;
s.camera.imSz=imSz;
s.camera.pixelSz=[pixelWidth,pixelHeight];
s.camera.sensorFormat=s.camera.imSz.*s.camera.pixelSz;
s.camera.focal=fx*s.camera.pixelSz(1);
s.camera.pp=[cx,cy].*s.camera.pixelSz;
s.camera.k=k; % TODO: Fix conversion to mm.
s.camera.p=p; % TODO: Fix conversion to mm.
s.camera.isFixed=sensorFixed;
s.camera.isAdjusted=isAdjusted;
s.camera.nominalFocal=nominalFocal;

% See what camera parameters have been estimated.
adjustedParams=struct('aspect',false,'cxcy',false(1,2),'f',false,...
                      'k',false(1,4),'p',false(1,4),'skew',false,...
                      'b',false(1,2));
% Parameters to warn about.
warnNotSupported={};
warnUsePhotoModeler={};
if isfield(chnk,'meta') && isfield(chnk.meta,'property')
    p=chnk.meta.property;
    if ~iscell(p)
        p={p};
    end
    % Extract property names.
    flds=cellfun(@(x)x.Attributes.name,p,'uniformoutput',false);
    
    % We should deal with all optimize/fit_XXX flds.
    stub='optimize/fit_';
    i=strncmp(stub,flds,length(stub));
    flds=flds(i);
    p=p(i);
    for i=1:length(flds)
        value=p{i}.Attributes.value;
        fullName=flds{i};
        param=fullName(14:end);
        switch param
          case 'flags'
            % Photoscan 1.2.4+ lists all estimated parameters in
            % optimize/fit_flags.
            
            % Split on whitespace
            j=[0,find(isspace(value)),length(value)+1];
            params=arrayfun(@(i)value(j(i)+1:j(i+1)-1),1:length(j)-1,...
                            'uniformoutput',false);
            % Value is true for each listed parameter.
            value='1';
          case 'cxcy'
            params={'cx','cy'};
          case 'k1k2k3'
            params={'k1','k2','k3'};
          case 'p1p2'
            params={'p1','p2'};
          otherwise
            % Photoscan pre-1.2.4 has multiple optimize/fit_PARAM fields.
            params={param};
        end
        for j=1:length(params)
            param=params{j};
            switch param
              case {'f','aspect','skew'}
                adjustedParams.(param)(:)=strcmp(value,'1');
              case 'cx'
                adjustedParams.cxcy(1)=strcmp(value,'1');
              case 'cy'
                adjustedParams.cxcy(2)=strcmp(value,'1');
              otherwise
                % Generic [bkp][1234]
                if length(param)==2 && ...
                        (param(1)=='b' && ismember(param(2),'12') || ...
                         ismember(param(1),'kp') && ismember(param(2),'1234'))
                    ix=sscanf(param(2),'%d');
                    adjustedParams.(param(1))(ix)=strcmp(value,'1');
                else
                    warning('Unknown camera calibration parameter: %s.',param);
                end
            end
            % Should we warn that a parameter not supported?
            switch param
              case {'aspect','skew','b1','b2','k4','p3','p4'}
                warnNotSupported{end+1}=param; %#ok<AGROW>
              case {'k1','k2','k3','p1','p2'}
                warnUsePhotoModeler{end+1}=param; %#ok<AGROW>
            end
        end
    end
end

s.camera.adjustedParams=adjustedParams;

if ~isempty(warnNotSupported)
    warning(['The following camera calibration parameters are not ' ...
             'supported yet:',sprintf(' %s',warnNotSupported{:})]);
end

if ~isempty(warnUsePhotoModeler)
    warning(['The following camera calibration parameters are supported, ' ...
             'but currently the Photomodeler lens distortion model ' ...
             'will be used:',sprintf(' %s',warnUsePhotoModeler{:})]);
end

% Delete unpacked files unless they should be kept.
if ~unpackLocal
    for i=length(dirs):-1:1
        delete(fullfile(dirs{i},'*'));
        rmdir(dirs{i});
    end
end

DelayedWaitBar('close');

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

function q=XformPtsi(p,M,stdToo)
%Apply 4-by-4 point transformation matrix M to points with ids p=[id,x,y,z].
%If stdToo is true, also transform stdx,stdy,stdz.

if nargin<3,stdToo=false; end

q=[p(:,1),XformPts(p(:,2:4)',M)'];
if stdToo
    std=M(1:3,1:3)*p(:,5:7)';
    q=[q,std'];
end

function Q=XformCams(P,M)
%Apply 4-by-4 point transformation matrix M to 3-by-4-by-K array P
%of camera matrices.

Q=nan(size(P));
for i=1:size(P,3)
    Q(:,:,i)=P(:,:,i)/M; % =*inv(M)
end


function defStd=getdefstd(chnk)
% Get default standard deviations from chunk settings.

tiePoints=nan;   % std for automatically detected tie points [pix]
projections=nan; % std for manually measured markers [pix]
markers=nan;     % std for marker positions [m]
camPos=nan;      % std for camera position [m]
camAng=nan;      % std for camera angles [deg]
scaleBars=nan;    % std for scale bar lengths [m]

defStd=struct('tiePoints',tiePoints,'projections',projections,...
              'markers',markers,'camPos',camPos,'camAng',camAng,...
              'scaleBars',scaleBars);

% Collect attribute names.
settingsProps=chnk.settings.property;
if ~iscell(settingsProps), settingsProps={settingsProps}; end
settingsPropNames=cellfun(@(x)x.Attributes.name,settingsProps,...
                          'uniformoutput',false);

% Field conversion table.
tbl={'tiepoints','tiePoints'
     'cameras','camPos'
     'cameras_ypr','camAng'
     'markers','markers'
     'scalebars','scaleBars'
     'projections','projections'};

for i=1:size(tbl,1)
    fld=['accuracy_',tbl{i,1}];
    ix=find(strcmp(settingsPropNames,fld));
    if length(ix)==1
        val=sscanf(settingsProps{ix}.Attributes.value,'%g');
    end
    defStd.(tbl{i,2})=val;
end


function DelayedWaitBar(varargin)
% Use:
% DELAYEDWAITBAR('init',delay,update,msg)
% Loop:
%    DELAYEDWAITBAR(val)
% DELAYEDWAITBAR('close')

persistent START LAPTIME H DELAY UPDATESTEP MSG

if ischar(varargin{1})
    switch varargin{1}
      case 'init'
        varargin(1)=[];
        % Defaults.
        DELAY=1;
        UPDATESTEP=1;
        MSG='';
        if length(varargin)>0 %#ok<ISMT>
            DELAY=varargin{1};
        end
        if length(varargin)>1
            UPDATESTEP=varargin{2};
        end
        if length(varargin)>2
            MSG=varargin{3};
        end
        START=clock;
        LAPTIME=START;
        H=[];
      case 'close'
        % Guard against user close.
        if ishandle(H)
            close(H)
        end
      otherwise
        error('Bad action');
    end
else
    % Numeric action, update or initialize bar.
    val=varargin{1};
    
    % Initialize if execution has taken more than DELAY seconds and
    % this is not the last update.
    if isempty(H) && etime(clock,START)>DELAY && val<1
        H=waitbar(val,MSG);
        LAPTIME=clock;
    elseif etime(clock,LAPTIME)>UPDATESTEP
        % Update dialog every UPDATESTEP s.
        if ishandle(H) % Guard against window close.
            waitbar(val,H);
        end
        LAPTIME=clock;
    end
end
