function s=loadlnz(lnzFile,varargin)
%LOADLNZ Load Photoscan .LNZ file.
%
%   S=LOADLNZ(FILE) loads the PhotoScan .LNZ file in FILE into a struct S.
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
%               isAdjusted     - true if the camera parameters are
%                                marked as adjusted by Photoscan.
%               givenParams    - struct with fields indicating what
%                                parameters were given in the lnz file:
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
%   By default, LOADLNZ unpacks the .LNZ file (a .ZIP archive) into a
%   directory in TEMPDIR and deletes the unpacked files after loading.
%   LOADLNZ(FILE,...,TRUE) instead unpacks the files into a local
%   subdir and does not delete the unpacked files. If FILE is
%   called PROJECT.LNZ, the subdir is called PROJECT_unpacked.
%   Additionally, LOADLNZ(FILE,TRUE,TRUE) creates ascii versions of
%   each .PLY file in a further 'ascii' subdir.
%
%   If both an initial and an adjusted camera is available in the .lnz
%   file, the adjusted camera is loaded.
%
%See also: TEMPDIR, TEMPNAME.

% Default values.
unpackLocal=false;

if length(varargin)>=1
    unpackLocal=varargin{1};
end

% Initialize waitbar to delay for 1s and update every 1s.
DelayedWaitBar('init',1,1,'Loading Photoscan lens project file...');

[lnzDir,lnzName,~]=fileparts(lnzFile);
if unpackLocal
    % Unpack to 'unpacked' subdir.
    unpackDir=fullfile(lnzDir,[lnzName,'_unpacked']);
else
    % Unpack to temporary dir.
    unpackDir=tempname;
end

% Unpack the .lnz file.
dirs=unpackpsz(lnzFile,unpackDir,false);
if unpackLocal
    fprintf('lnz files unpacked to directory %s.\n',unpackDir);
end
DelayedWaitBar(0.25);

% Load project data from the main xml file.
fName=fullfile(unpackDir,'doc.xml');
s=dbatxml2struct(fName);
DelayedWaitBar(0.3);

s.fileName=lnzFile;

s.defStd=struct('tiePoints',1,...
                'projections',0.1,...
                'markers',0.005,...
                'camPos',10,...
                'camAng',2,...
                'scaleBars',1e-3);

% Get the group holding the data.
group=s.document.group;

% Determine what photos we have.
photo=group.photo;
if ~iscell(photo)
    photo={photo};
end

% Extract transformation.
% Transformations are from "image" coordinate system to local.
xforms=nan(4,4,length(photo));
% Photo matrices from local coordinates.
P=nan(3,4,length(photo));
% Photo centers in local coordinates.
CC=nan(3,length(photo));
priorCC=CC;
for i=1:length(photo)
    T=reshape(sscanf(photo{i}.transform.Text,'%g '),4,4)';
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
s.raw.priorCC=priorCC;

s.local.P=s.raw.P;
s.local.CC=s.raw.CC;
s.local.priorCC=s.raw.priorCC;
s.local.R=nan(3,3,size(s.local.P,3));
for i=1:size(s.local.R,3)
    R=s.local.P(:,1:3,i);
    s.local.R(:,:,i)=R/det(R)^(1/3);
end

DelayedWaitBar(0.35);

% Process image files names. File names will be sorted as s.photoId.
imNames=cell(1,length(photo));
for i=1:length(photo)
    % Guess  if image path is absolute or relative.
    p=photo{i}.location.Attributes.path;
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
        imNames{i}=p;
    else
        imNames{i}=fullfile(lnzDir,p);
    end
end
s.imNames=imNames;

cameraIds=0:length(imNames)-1;;
s.cameraIds=cameraIds;
s.cameraLabels=cell(size(imNames));
for i=1:length(imNames)
    [p,n,e]=fileparts(imNames{i});
    s.cameraLabels{i}=[n,e];
end
s.cameraEnabled=true(size(imNames));

invCameraIds=nan(max(cameraIds)+1,1);
invCameraIds(cameraIds+1)=1:length(cameraIds);

% Functions to convert between Photoscan camera id and DBAT camera number.
PSCamId=@(id)cameraIds(id);
DBATCamId=@(id)invCameraIds(id+1);
s.PSCamId=PSCamId;
s.DBATCamId=DBATCamId;


% Get image size, resolution, focal length.
w=ParseMeta(photo,'width','%d');
h=ParseMeta(photo,'height','%d');
f=ParseMeta(photo,'flength','%g');
xRes=ParseMeta(photo,'fplane_xres','%g');
yRes=ParseMeta(photo,'fplane_yres','%g');

if ~isscalar(w)
    warning('No or no unique image width');
end
if ~isscalar(h)
    warning('No or no unique image height');
end
if ~isscalar(f)
    warning('No or no unique focal length');
end
if ~isscalar(xRes)
    warning('No or no unique x resolution');
end
if ~isscalar(yRes)
    warning('No or no unique y resolution');
end
% Compute sensor size.
sw=w/xRes;
sh=h/yRes;

% Package camera information.
imSz=[w,h];
sensorFormat=[sw,sh];
pixelSz=1./[xRes,yRes];
nominalFocal=f;

s.camera=struct('imSz',imSz,'sensorFormat',sensorFormat,'pixelSz',pixelSz,...
                'nominalFocal',nominalFocal,'focal',nan,'pp',nan(1,2),'k',[],'p',[]);

% Process measured corners.
cc=cell(size(photo));

for i=1:length(photo)
    p=photo{i};
    m=nan(length(p.corner),6);
    m(:,1)=i;
    m(:,2)=cellfun(@(x)sscanf(x.Attributes.img_x,'%g'),p.corner)';
    m(:,3)=cellfun(@(x)sscanf(x.Attributes.img_y,'%g'),p.corner)';
    m(:,4)=cellfun(@(x)sscanf(x.Attributes.obj_x,'%g'),p.corner)';
    m(:,5)=cellfun(@(x)sscanf(x.Attributes.obj_y,'%g'),p.corner)';
    m(:,6)=cellfun(@(x)strcmp(x.Attributes.valid,'true'),p.corner)';
    cc{i}=m;
end

% Concatenate corners from all images.
corners=cat(1,cc{:});
% Only keep valid corners.
corners=corners(corners(:,6)~=0,:);

% Create unique ID based on object coordinates.
[uc,~,ic1]=unique(corners(:,4:5),'rows');

% Create visibility matrix.
s.vis=sparse(ic1,corners(:,1),corners(:,6));

ctrlPts=[(1:size(uc,1))',uc,zeros(size(uc,1),1),zeros(size(uc,1),3)];
s.local.ctrlPts=ctrlPts;
s.local.ctrlPtsLabels=cell(size(ctrlPts,1),1);
for i=1:size(ctrlPts,1)
    s.local.ctrlPtsLabels{i}=sprintf('(%d,%d)',ctrlPts(i,2:3));
end
s.local.ctrlPtsEnabled=sum(s.vis,2)>0;
s.local.objPts=zeros(0,4);
s.raw.ctrlPtsLabels=s.local.ctrlPtsLabels;
s.raw.ctrlPts=s.local.ctrlPts;

s.global=s.local;

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

% Highest DBAT CP id.
maxDBATCPid=length(rawCPids);

% Map object point ids to above control point ids. Generate NaN's for
% out-of-bounds OP ids.
DBATOPid=@(id)id+1+maxDBATCPid+0./(id>=0);
PSOPid=@(id)id-1-maxDBATCPid+0./(id>=1+maxDBATCPid);
s.DBATOPid=DBATOPid;
s.PSOPid=PSOPid;

s.markPts.obj=zeros(0,4);
s.markPts.sz=nan(0,1);
s.markPts.ctrl=[corners(:,1),ic1,corners(:,2:3)];
[~,i]=msort(s.markPts.ctrl(:,1:2));
s.markPts.ctrl=s.markPts.ctrl(i,:);
s.markPts.ctrlPinned=true(size(s.markPts.ctrl,1),1);
s.markPts.IsPinned=@(x)true;
s.markPts.all=s.markPts.ctrl;

% Delete unpacked files unless they should be kept.
if ~unpackLocal
    for i=length(dirs):-1:1
        delete(fullfile(dirs{i},'*'));
        rmdir(dirs{i});
    end
end

DelayedWaitBar('close');

function v=ParseMeta(photo,name,fmt)
% Filter photo{:}.meta{:}.Attributes.name for name and parse
% ...Attributes.value using fmt. Return vector of unique values.

v=nan(size(photo));

for i=1:length(photo)
    % Extract meta fields. Guard against single field.
    m=photo{i}.meta;
    if ~iscell(m)
        m={m};
    end
    % Extract field names.
    metaNames=cellfun(@(x)x.Attributes.name,m,'uniformoutput',false);

    % Find wanted field.
    j=find(strcmp(metaNames,name));
    if isscalar(j)
        v(i)=sscanf(m{j}.Attributes.value,fmt);
    end
end

v=unique(v(~isnan(v)));

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
