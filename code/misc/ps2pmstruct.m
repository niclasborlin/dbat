function [prob,pmReport,pts3d,pts2d]=ps2pmstruct(s,useSemiLocal)
%PS2PMSTRUCT Convert PhotoScan structure to PhotoModeler structure.
%
%   PROB=PS2PMSTRUCT(S) converts the PhotoScan struct S, as returned
%   by LOADPSZ, to a PhotoModeler struct PROB, as loaded by
%   LOADPM. All EO, OP parameters are in global coordinates. Use
%   PS2PMSTRUCT(S,TRUE) to use semilocal parameters instead
%   (translation and scaling from global, but no rotation).
%
%   [PROB,PTS3D,PTS2D]=PS2PMSTRUCT(S) furthermore returns 3D and 2D
%   points.
%
%   Zero-based object point ids in S will be converted to one-based
%   ids in PROB.
%
%   Camera positions, object and control points are set to global
%   coordinates. Object point uncertainties are set to unknown. Mark
%   points are given a 1 pixel standard deviation.
%
%See also: LOADPSZ, LOADPM.

if nargin<2, useSemiLocal=false; end

% Create a fake job header with default camera.
imSz=s.camera.imSz(:);
% Lens distortion parameters.
k1k3=zeros(3,1);
p1p2=zeros(2,1);

if length(s.camera.k)<=length(k1k3)
    k1k3(1:length(s.camera.k))=s.camera.k(:);
else
    k1k3(1:end)=s.camera.k(1:length(k1k3));
    warning('Ignoring K4');
end

if length(s.camera.p)<=length(p1p2)
    p1p2(1:length(s.camera.p))=s.camera.p(:);
else
    p1p2(1:end)=s.camera.p(1:length(p1p2));
    warning('Ignoring P3, ...');
end

defCam=[s.camera.focal;s.camera.pp(:);s.camera.sensorFormat(:);k1k3;p1p2];

job=struct('fileName',s.fileName,'title','Photoscan import','defCam',defCam,'defCamStd',zeros(size(defCam)),'imSz',imSz);

if useSemiLocal
    pos=s.semilocal;
else
    pos=s.global;
end

% Create EO parameters.

% Copy camera positions.
CC=pos.CC;

% Convert rotation matrices to omega-phi-kappa.
ang=nan(size(CC));
for i=1:size(ang,2)
    RR=pos.R(:,:,i);
    ang(:,i)=derotmat3d(RR)';
end
angPM=ang([3,2,1],:)*180/pi;

images=repmat(struct('imName','','outer',nan(1,6),'outerStd',zeros(1,6),...
                     'imSz',imSz,'id',nan,'label',''),size(CC,2),1);
for i=1:length(images)
    images(i).imName=s.imNames{i};
    images(i).outer=[CC(:,i);angPM(:,i)]';
    images(i).id=s.cameraIds(i);
    images(i).label=s.cameraLabels{i};
end

% Enabled control points are control points. Disabled control
% points are check points.
enabled=pos.controlPts.enabled;
ctrlPts=[pos.controlPts.id(enabled),pos.controlPts.pos(:,enabled)',...
         pos.controlPts.std(:,enabled)'];
checkPts=[pos.controlPts.id(~enabled),pos.controlPts.pos(:,~enabled)',...
          pos.controlPts.std(:,~enabled)'];
% Remove check points with <2 rays.
ok=true(size(checkPts,1),1);
for i=1:size(checkPts,1)
    if nnz(s.markPts.ctrl(:,2)==checkPts(i,1))<2
        ok(i)=false;
    end
end
checkPts=checkPts(ok,:);

% Create prior camera positions.
hasPriorCamPos=~any(isnan(pos.priorCamPos.pos),1);
priorCamPos=[s.cameraIds(hasPriorCamPos);
             pos.priorCamPos.pos(:,hasPriorCamPos);
             pos.priorCamPos.std(:,hasPriorCamPos)]';

% Track original ids and labels.
rawCPids=s.PSCPid(ctrlPts(:,1));
CPlabels=pos.controlPts.labels(enabled);
rawCCPids=s.PSCPid(checkPts(:,1));
CCPlabels=pos.controlPts.labels(~enabled);
CCPlabels=CCPlabels(ok);

% Merge and sort ctrl and check pts.
cPts=[ctrlPts;checkPts];
rawCids=[rawCPids;rawCCPids];
cLabels=cat(2,CPlabels,CCPlabels);
[~,i]=sort(cPts(:,1));
cPts=cPts(i,:);
rawCids=rawCids(i);
cLabels=cLabels(i);

% Copy global object points. Set posterior uncertainty to unknown.
objPts=[cPts;[pos.objPts,nan(size(pos.objPts,1),3)]];

% Track original object point IDs.
rawOPids=[rawCids;s.PSOPid(pos.objPts(:,1))];
OPlabels=[cLabels,repmat({''},1,size(pos.objPts,1))];

% Copy mark points and set std.
ctrlStd=s.defStd.projections;
objStd=s.defStd.tiePoints;
markPts=[s.markPts.ctrl,repmat(ctrlStd,size(s.markPts.ctrl,1),2);
         s.markPts.obj,repmat(objStd,size(s.markPts.obj,1),2)];

% Sort by image, then id.
markPts=msort(markPts);

% Remove any mark points not among 3D points.
if ~isempty(markPts)
    keep=ismember(markPts(:,2),objPts(:,1));
    markPts=markPts(keep,:);
end

if ~isempty(markPts)
    % Convert image indices to zero-based for PhotoModeler compatibility.
    markPts(:,1)=markPts(:,1)-1;
end

% Construct PM structure.
prob=struct('job',job,'images',images,'ctrlPts',ctrlPts,'checkPts', checkPts, ...
            'objPts',objPts,'priorCamPos',priorCamPos, 'rawOPids',rawOPids, ...
            'OPlabels',{OPlabels}, 'markPts',markPts);

EO=[CC;ang;zeros(1,size(CC,2))];
EOstd=nan(size(EO));
    
pmReport=struct('EO',EO,'EOstd',EOstd);

% Create tables of 3D and 2D info.
id3d=zeros(1,0);
if ~isempty(pos.controlPts.id)
    id3d=[id3d,pos.controlPts.id'];
end
if ~isempty(pos.objPts)
    id3d=[id3d,pos.objPts(:,1)'];
end
name=repmat({''},size(id3d));
name(s.DBATCPid(s.raw.controlPts.id))=s.raw.controlPts.labels;
pos3d=zeros(3,0);
if ~isempty(pos.controlPts.pos)
    pos3d=[pos3d,pos.controlPts.pos];
end
if size(pos.objPts,2)>=4
    pos3d=[pos3d,pos.objPts(:,2:4)'];
end
std3d=zeros(3,0);
if ~isempty(pos.controlPts.std)
    std3d=[std3d,pos.controlPts.std];
else
    std3d=[std3d,nan(size(pos.controlPts.pos))];
end
if size(pos.objPts,2)>=7
    std3d=[std3d,pos.objPts(:,5:7)'];
else
    std3d=[std3d,nan(size(pos.objPts,1),3)'];
end
% Create a 'raw' visibility map, i.e. for all object points.
rawVis=sparse(s.markPts.all(:,2),s.markPts.all(:,1),1);
% Keep only the object points that were actually computed.
vis=rawVis(id3d,:);

% Sort 3d pts by id.
[~,i]=sort(id3d);
id3d=id3d(i);
name=name(i);
vis=vis(i,:);
pos3d=pos3d(:,i);
std3d=std3d(:,i);
pts3d=struct('id',id3d,'name',{name},'vis',logical(vis),'pos',pos3d,...
             'std',std3d);

% Sort 2d pts by image first, then id.
s.markPts.all=sortrows(s.markPts.all,[1,2]);

% Only keep 2d measurements that corresponds to 3D points.
keep=ismember(s.markPts.all(:,2),pts3d.id);

id2d=s.markPts.all(keep,2)';
imNo=s.markPts.all(keep,1)';
pos2d=s.markPts.all(keep,3:4)';
res2d=nan(size(pos2d));

pts2d=struct('id',id2d,'imNo',imNo,'pos',pos2d,'res',res2d,'kept',keep);
