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
%               controlPts - struct with arrays
%                 id      - 1-by-MC with ids
%                 pos     - 3-by-MC with positions
%                 std     - 3-by-MC with standard deviations
%                 cov     - 3-by-3-by-MC with full covariance matrices
%                 enabled - 1-by-MC with logical
%                 labels  - 1-by-MC cell array with strings
%               objPts  - MO-by-4 array with [id,x,y,z] for object pts,
%               P       - 3-by-4-by-N array with camera matrices,
%               CC      - 3-by-N array with camera centers,
%               R       - 3-by-3-by-N array with camera rotation matrices.
%   cameraIds - N-vector with camera ids,
%   cameraLabels - N-cell vector with camera labels,
%   cameraEnabled - logical N-vector indicating which cameras are enabled,
%   cameraOriented - logical N-vector indicating which cameras are oriented,
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
%               isAdjusted     - true if the camera parameters are
%                                marked as adjusted by Photoscan.
%               givenParams    - struct with fields indicating what
%                                parameters were given in the psz file:
%                                f      - scalar
%                                cxcy   - 2-vector
%                                aspect - scalar
%                                skew   - scalar
%                                b      - 2-vector
%                                k      - 4-vector
%                                p      - 4-vector
%               optimizedParams - struct with same fields as givenParams,  
%                                indicating which parameters have
%                                been optimized by Photoscan. 
%   defStd    - struct with default standard deviations
%               tiePoints   - std for automatically detected tie points [pix]
%               projections - std for manually measured markers [pix]
%               markers     - std for marker positions [m]
%               camPos      - std for camera position [m]
%               camAng      - std for camera angles [deg]
%               scaleBars   - std for scale bar lengths [m]
%
%   By default, LOADPSZ loads enabled and oriented images. Use
%   LOADPSZ(FILE,TRUE) to also load unoriented images.
%
%   By default, LOADPSZ unpacks the .PSZ file (a .ZIP archive) into a
%   directory in TEMPDIR and deletes the unpacked files after loading.
%   LOADPSZ(FILE,...,'keep') instead unpacks the files into a local
%   subdir and does not delete the unpacked files. If FILE is called
%   PROJECT.PSZ, the subdir is called PROJECT_unpacked. Additionally,
%   LOADPSZ(FILE,...,'ascii') creates ascii versions of each .PLY file
%   in a further 'ascii' subdir.
%
%   If both an initial and an adjusted camera is available in the .psz
%   file, the adjusted camera is loaded.
%
%See also: TEMPDIR, TEMPNAME, UNPACKPSZ.

% Highest tested version.
highestTestedVersion='1.4.0'; %#ok<NASGU>

% Default values.
chunkNo=1;
unpackLocal=false;
asciiToo=false;
keepUnoriented=false;

% Handle arguments.
while ~isempty(varargin)
    if islogical(varargin{1})
        keepUnoriented=varargin{1};
    elseif isnumeric(varargin{1})
        chunkNo=varargin{1};
    else
        switch varargin{1}
          case 'keep'
            unpackLocal=true;
          case 'ascii'
            unpackLocal=true;
            asciiToo=true;
          otherwise
            error('%s: Bad argument %s',mfilename,varargin{1});
        end
    end
    varargin(1)=[];
end

% Initialize waitbar to delay for 1s and update every 1s.
DelayedWaitBar('init',1,1,'Loading Photoscan project file...');

if ~exist(psFile,'file')
    error('File %s does not exists',psFile);
end

[psDir,psName,~]=fileparts(psFile);
if unpackLocal
    % Unpack to 'unpacked' subdir.
    unpackDir=fullfile(psDir,[psName,'_unpacked']);
else
    % Unpack to temporary dir.
    unpackDir=tempname;
end

% Unpack the .psz file.
dirs=unpackpsz(psFile,unpackDir,asciiToo);
if unpackLocal
    fprintf('psz files unpacked to directory %s.\n',unpackDir);
end
DelayedWaitBar(0.25);

% Load project data from the main xml file.
fName=fullfile(unpackDir,'doc.xml');
s=dbatxml2struct(fName);
DelayedWaitBar(0.3);

% Extract document version.
s.version='0.0.0';
if isfield(s,'document') && ...
        isfield(s.document,'Attributes') && ...
        isfield(s.document.Attributes,'version')
    s.version=s.document.Attributes.version;
end

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

if isfield(chnk,'reference') && isfield(chnk.reference,'Text')
    % Warn unless the coordinate system is local.
    prefix='local_cs[';
    isLocal=strncmpi(chnk.reference.Text,prefix,length(prefix));
    if ~isLocal
        warning('%s: Non-local coordinate system. Bundle may fail: %s',...
                mfilename,chnk.reference.Text);
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
cameraEnabled=false(size(camera));
for i=1:length(cameraEnabled)
    cameraEnabled(i)=ParseTF(camera{i}.Attributes.enabled,...
                             sprintf('Unknown status for camera %d',...
                                     cameraIds(i)));
end

if length(unique(sensorIds(cameraEnabled)))>1
    error('Handling of cameras for multiple sensor ids not implemented yet');
end

% Extract transformation.
% Transformations are from "image" coordinate system to local.
xforms=nan(4,4,length(cameraIds));
% Camera matrices from local coordinates.
P=nan(3,4,length(cameraIds));
% Camera centers in local coordinates.
CC=nan(3,length(cameraIds));

% Prior observations of camera centers
priorCamPos=struct('pos',nan(3,length(cameraIds)),...
                   'std',nan(3,length(cameraIds)),...
                   'cov',nan(3,3,length(cameraIds)),...
                   'enabled',false(1,length(cameraIds)));
for i=1:length(cameraIds)
    if isfield(camera{i},'transform')
        T=reshape(sscanf(camera{i}.transform.Text,'%g '),4,4)';
        xforms(:,:,i)=T;
        if 1
            % TODO: Check this "mirroring"...
            P(:,:,i)=eye(3,4)/(T*diag([1,-1,-1,1]));
        else
            warning('Untested non-mirroring'); %#ok<UNRCH>
            P(:,:,i)=eye(3,4)/T; % *inv(T)
        end
        CC(:,i)=euclidean(null(P(:,:,i)));
    end
    
    % Check if we have reference EO coordinates.
    if isfield(camera{i},'reference')
        attr=camera{i}.reference.Attributes;
        % TODO: Always load, recognize enable/disable status.
        priorCamPos.enabled(i)=...
            ParseTF(attr.enabled,...
                    sprintf('Unknown status for reference for camera %d (%s)',...
                            cameraIds(i),cameraLabels{i}));
        [priorCamPos.pos(:,i),...
         priorCamPos.std(:,i),...
         priorCamPos.cov(:,:,i)]=ParseReferencePos(attr,s.defStd.camPos);
    end
end

cameraOriented=all(isfinite(CC),1);

% Keep enabled cameras.
keep=cameraEnabled;
if ~keepUnoriented
    keep=keep & cameraOriented;
end

s.cameraIds=cameraIds(keep);
s.cameraLabels=cameraLabels(keep);
s.cameraEnabled=cameraEnabled(keep);
s.cameraOriented=cameraOriented(keep);

camIds=cameraIds(keep);

% Functions to convert between Photoscan camera id and DBAT camera number.
PSCamId=@(dbatId)IDLookup(camIds,dbatId);
DBATCamId=@(psCamId)IDInvLookup(camIds,psCamId);
s.PSCamId=PSCamId;
s.DBATCamId=DBATCamId;

s.raw.transforms=xforms(:,:,keep);
s.raw.P=P(:,:,keep);
s.raw.CC=CC(:,keep);
priorCamPos.pos=priorCamPos.pos(:,keep);
priorCamPos.std=priorCamPos.std(:,keep);
priorCamPos.cov=priorCamPos.cov(:,:,keep);
priorCamPos.enabled=priorCamPos.enabled(keep);
s.raw.priorCamPos=priorCamPos;

s.local.P=s.raw.P;
s.local.CC=s.raw.CC;
s.local.priorCamPos=XformPos(s.raw.priorCamPos,G2L);
s.local.R=nan(3,3,size(s.local.P,3));
for i=1:size(s.local.R,3)
    R=s.local.P(:,1:3,i);
    s.local.R(:,:,i)=R/det(R)^(1/3);
end

s.global.P=XformCams(s.local.P,L2G);
s.global.CC=XformPts(s.local.CC,L2G);
s.global.priorCamPos=s.raw.priorCamPos;
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
s.semilocal.priorCamPos=XformPos(s.global.priorCamPos,G2SL);
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

% Pre-allocate control point information.
controlPts=struct('id',nan(1,length(markers)),...
                  'pos',nan(3,length(markers)),...
                  'std',nan(3,length(markers)),...
                  'cov',nan(3,3,length(markers)),...
                  'enabled',false(1,length(markers)),...
                  'labels',{cell(1,length(markers))});

for i=1:length(markers)
    m=markers{i};
    id=sscanf(m.Attributes.id,'%d');
    controlPts.id(i)=id;
    if isfield(m.Attributes,'label')
        controlPts.labels{i}=m.Attributes.label;
    end
    pos=nan(3,1);
    st=nan(3,1);
    cc=nan(3);
    if isfield(m,'reference') && isfield(m.reference,'Attributes')
        [pos,st,cc]=ParseReferencePos(m.reference.Attributes,s.defStd.markers);
        if isfield(m.reference.Attributes,'enabled')
            controlPts.enabled(i)=ParseTF(m.reference.Attributes.enabled,...
                                          sprintf('Unknown status for ctrl pt %d',id));
        end
    end
    controlPts.pos(:,i)=pos;
    controlPts.std(:,i)=st;
    controlPts.cov(:,:,i)=cc;
end
s.raw.controlPts=controlPts;

DelayedWaitBar(0.4);

% Make local/global ctrl pt ids 1-based.
rawCPids=controlPts.id(:);

invCPids=nan(max(rawCPids)+1,1);
invCPids(rawCPids+1)=1:length(rawCPids);

% Will convert a zero-based id to a one-based id. Generate NaN's for
% out-of-bounds CP ids.
PSCPid=@(id)rawCPids(id);
DBATCPid=@(id)invCPids(id+1);

s.DBATCPid=DBATCPid;
s.PSCPid=PSCPid;

% Copy raw ctrl pts and adjust id.
s.global.controlPts=s.raw.controlPts;
s.global.controlPts.id=DBATCPid(s.global.controlPts.id);
% Transform ctrl pts from global to local and semilocal coordinate systems.
s.local.controlPts=XformPos(s.global.controlPts,G2L);
s.semilocal.controlPts=XformPos(s.global.controlPts,G2SL);

% Highest DBAT CP id.
maxDBATCPid=length(rawCPids);

% Map object point ids to above control point ids. Generate NaN's for
% out-of-bounds OP ids.
DBATOPid=@(id)id+1+maxDBATCPid+0./(id>=0);
PSOPid=@(id)id-1-maxDBATCPid+0./(id>=1+maxDBATCPid);
s.DBATOPid=DBATOPid;
s.PSOPid=PSOPid;

% Copy raw object points and adjust id.
s.local.controlPts=s.raw.controlPts;
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

% Remove any projections not to load.
projCameraIds=cellfun(@(x)sscanf(x.Attributes.camera_id,'%d'),projs);
keepProjs=ismember(projCameraIds,s.cameraIds);
projs=projs(keepProjs);
projCameraIds=projCameraIds(keepProjs);

% Projections will be sorted to match s.cameraId.
projections=cell(size(projs));
s.raw.paths.projections=cell(size(projs));

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
hasProj=~cellfun(@isempty,projections);
if all(hasProj)
    nPts=cellfun(@(x)length(x.vertex.id),projections);
else
    nPts=zeros(size(hasProj));
    nPts(hasProj)=cellfun(@(x)length(x.vertex.id),projections(hasProj));
end
ptIx=cumsum([0,nPts]);
objMarkPts=nan(sum(nPts),4);
objKeyPtSize=nan(sum(nPts),1);

% Collect all measured points.
for i=find(nPts)
    ni=nPts(i);
    % Index for where to put the points.
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
        % Filter out markers in ignored images.
        keep=~isnan(DBATCamId(camIds));
        
        if ~any(keep)
            % No measurements in kept images.
            continue;
        end
        
        location=location(keep);
        camIds=camIds(keep);

        % Extract coordinates.
        x=cellfun(@(m)sscanf(m.Attributes.x,'%g'),location);
        y=cellfun(@(m)sscanf(m.Attributes.y,'%g'),location);
        pinned=false(size(location));
        for j=1:length(pinned)
            if isfield(location{j}.Attributes,'pinned')
                pinned(j)=ParseTF(location{j}.Attributes.pinned,...
                                  sprintf(['Unknown pinned status ' ...
                                    'for marker %d in camera %d'],...
                                          markerId,camIds(j)));
            end
        end
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

% Filter out unwanted images.
camId=cellfun(@(x)sscanf(x.Attributes.camera_id,'%d'),camera);
keep=~isnan(DBATCamId(camId));
camId=camId(keep);
camera=camera(keep);
imNames=cell(1,length(camera));
for i=1:length(camera)
    % Convert to DBAT camera number.
    j=DBATCamId(camId(i));
    % Guess  if image path is absolute or relative.
    p=camera{i}.photo.Attributes.path;
    if isempty(p)
        % Empty paths (should never happen) are defined to be absolute.
        pathIsAbsolute=true;
    elseif ismember(p(1),'/\')
        % Path starting with / or \ are absolute.
        pathIsAbsolute=true;
    else
        % Path did not start with / or \. Determine if its Z:/
        if length(p)>1 && p(2)==':'
            pathIsAbsolute=true;
        else
            pathIsAbsolute=false;
        end
    end
    if pathIsAbsolute
        imNames{j}=p;
    else
        imNames{j}=fullfile(psDir,p);
    end
end
s.imNames=imNames;

% What camera parameters are given?
givenParams=struct('aspect',false,'cxcy',false(1,2),'f',false,...
                   'k',false(1,4),'p',false(1,4),'skew',false,...
                   'b',false(1,2));
% What camera parameters are listed as optimized?
optimizedParams=givenParams;

% Determine which sensor we want
wantedSensorId=unique(sensorIds);
sensors=chnk.sensors.sensor;
if ~iscell(sensors)
    sensors={sensors};
end
sensorId=cellfun(@(x)sscanf(x.Attributes.id,'%d'),sensors);
keep=sensorId==wantedSensorId;
sensor=sensors{keep};

% Collect calibrated camera parameters.
cal=sensor.calibration;
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

% cx, cy are relative to the image center if 'f' is specified and
% aboslute if 'fx', 'fy' are specified.
ppIsAbsolute=false;

% Parse focal length(s)
if isfield(cal,'fx')
    fx=sscanf(cal.fx.Text,'%g');
    ppIsAbsolute=true;
else
    fx=nan;
end
if isfield(cal,'fy')
    fy=sscanf(cal.fy.Text,'%g');
    ppIsAbsolute=true;
else
    fy=nan;
end
if isfield(cal,'f')
    ppIsAbsolute=false;
    fy=sscanf(cal.f.Text,'%g');
    if isfield(cal,'b1')
        b1=sscanf(cal.b1.Text,'%g');
    else
        b1=0;
    end
    fx=fy+b1;
end
if isfield(cal,'cx')
    cx=sscanf(cal.cx.Text,'%g');
else
    cx=0;
end
if isfield(cal,'cy')
    cy=sscanf(cal.cy.Text,'%g');
else
    cy=0;
end
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

% f and cxcy are required.
givenParams.f=true;
givenParams.cxcy(:)=true;
% Mark supplied K and P as given.
givenParams.k(1:length(k))=true;
givenParams.p(1:length(p))=true;

skew=0;
% Is skew given?
if isfield(cal,'skew')
    skew=sscanf(cal.skew.Text,'%g');
    givenParams.skew=true;
elseif isfield(cal,'b2')
    skew=sscanf(cal.b2.Text,'%g');
    givenParams.skew=true;
end

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

% File version before 1.4.0 had raw cx/cy values. From 1.4.0, cx
% and cy are w.r.t the image center.
if ~ppIsAbsolute
    cx=cx+imSz(1)/2;
    cy=cy+imSz(2)/2;
end

% Construct camera calibration matrix.
K=[fx,skew,cx;0,fy,cy;0,0,1];

s.K=K;

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
    sensorFixed=ParseTF(fProp{1}.Attributes.value,...
                        'Unknown sensorfixed property');
end

if (fx*pixelWidth~=fy*pixelHeight)
    % We always have fx, fy. Assume that if fx==fy, aspect is
    % locked at unity. Otherwise assume it has been estimated.
    givenParams.aspect=true;
end

s.camera.name=camName;
s.camera.type=camType;
s.camera.imSz=imSz;
s.camera.pixelSz=[pixelWidth,pixelHeight];
s.camera.sensorFormat=s.camera.imSz.*s.camera.pixelSz;
s.camera.focal=fx*s.camera.pixelSz(1);
s.camera.pp=[cx,cy].*s.camera.pixelSz;
% Deal with scaling of distortion coefficients
s.camera.k=-k.*s.camera.focal.^(-2*(1:length(k)));
pScaled=p;
if ~isempty(pScaled)
    pScaled(1:2)=pScaled(1:2)/s.camera.focal;
    pScaled(2)=-pScaled(2);
    pScaled(1:2)=pScaled([2,1]);
end
s.camera.p=pScaled;
s.camera.isFixed=sensorFixed;
s.camera.isAdjusted=isAdjusted;
s.camera.nominalFocal=nominalFocal;

if fx*s.camera.pixelSz(1)~=fy*s.camera.pixelSz(2)
    warning('Non-square pixel size not currently supported')
end

s.camera.givenParams=givenParams;

% Parameters to warn about.
warnNotSupported={};

% Have we found any optimize/fit_XXX fields?
optFitTagsFound=false;

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
    optFitTagsFound=~isempty(flds);
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
                optimizedParams.(param)(:)=strcmp(value,'1');
              case 'cx'
                optimizedParams.cxcy(1)=strcmp(value,'1');
              case 'cy'
                optimizedParams.cxcy(2)=strcmp(value,'1');
              otherwise
                % Generic [bkp][1234]
                if length(param)==2 && ...
                        (param(1)=='b' && ismember(param(2),'12') || ...
                         ismember(param(1),'kp') && ismember(param(2),'1234'))
                    ix=sscanf(param(2),'%d');
                    optimizedParams.(param(1))(ix)=strcmp(value,'1');
                else
                    warning('Unknown camera calibration parameter: %s.',param);
                end
            end
            % Should we warn that a parameter not supported?
            switch param
              case {'aspect','skew','b1','b2','k4','p3','p4'}
                warnNotSupported{end+1}=param; %#ok<AGROW>
            end
        end
    end
end

s.camera.optimizedParamsFound=optFitTagsFound;
s.camera.optimizedParams=optimizedParams;

if ~isempty(warnNotSupported)
    warning(['The following camera calibration parameters are not ' ...
             'supported yet:',sprintf(' %s',warnNotSupported{:})]);
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


function q=XformPos(p,M)
%Apply 4-by-4 point transformation matrix M to point structure p.
%The structure has fields pos, std, and cov.

R=M(1:3,1:3);

q=p;
q.pos=euclidean(M*homogeneous(p.pos));
for i=1:size(p.pos,2)
    q.std(:,i)=diag(sqrt(R*diag(p.std(:,i).^2)*R'));
    q.cov(:,:,i)=R*p.cov(:,:,i)*R';
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


function j=IDLookup(tbl,i)
% Do a forward id lookup in the vector TBL, i.e., return TBL(I). If I
% is outside the vector, return NaN.

j=nan(size(i));
iOk=i>=1 & i<=length(tbl);
j(iOk)=tbl(i(iOk));


function i=IDInvLookup(tbl,j)
% Do an inverse id lookup in the vector tbl, i.e., the position i
% such that tbl(i)==j. If j is not found, return NaN.

i=nan(size(j));
for jj=1:length(j)
    ii=find(tbl==j(jj),1);
    if isscalar(ii)
        i(jj)=ii;
    end
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


function z=CompareVersion(s,t) %#ok<DEFNU>
%Compare version strings. Return -1, 0 or +1 if s is lower, equal
%to, or higher than t.
%
%S and T must be in digits separated by dots, e.g. 1.4.0.
%The comparison is done section by section.

sDigits=sscanf(s,'%d.');
tDigits=sscanf(t,'%d.');

% Make versions have the same number of digits.
l=max(length(sDigits),length(tDigits));

z=0;

for i=1:l
    if i<=length(sDigits)
        sd=sDigits(i);
    else
        sd=0;
    end
    if i<=length(tDigits)
        td=tDigits(i);
    else
        td=0;
    end
    if sd>td
        z=1;
        return;
    end
    if sd<td
        z=-1;
        return;
    end
end


function [p,st,cc]=ParseReferencePos(s,defStd)
%Parse x/y/z + std values for a reference struct.
%
%s contains the struct to be parsed. Fields that are parsed:
%  x, y, z, sx, sy, sz, sxy, sxyz
%defStd is the default standard deviation if none is specified in
%  the struct. Can be scalar or 3-by-1.
%3-by-1 p returns the position.
%3-by-1 st returns the standard deviations.
%3-by-3 cc returns any covariance matrix, if specified.

p=nan(3,1);
st=nan(3,1);
if isscalar(defStd)
    st(:)=defStd;
elseif length(defStd)==3
    st(:)=defStd(1:3);
end

if isfield(s,'x')
    p(1)=sscanf(s.x,'%g');
end
if isfield(s,'y')
    p(2)=sscanf(s.y,'%g');
end
if isfield(s,'z')
    p(3)=sscanf(s.z,'%g');
end
if isfield(s,'sxyz')
    st(:)=sscanf(s.sxyz,'%g');
elseif isfield(s,'sxy')
    st(1:2)=sscanf(s.sxy,'%g');
end
if isfield(s,'sx')
    st(1)=sscanf(s.sx,'%g');
end    
if isfield(s,'sy')
    st(2)=sscanf(s.sy,'%g');
end    
if isfield(s,'sz')
    st(3)=sscanf(s.sz,'%g');
end
cc=diag(st.^2);


function t=ParseTF(s,msg)
%Parse the string s containing 0/1 or false/true. Return false or true.
%Uses msg to generate an error message if s is bad.

switch s
  case {'1','true'}
    t=true;
  case {'0','false'}
    t=false;
  otherwise
    warning('%s: Trying to parse true/false string, got %s.',msg,s);
end
