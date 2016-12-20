function [prob,pmReport,pts3d,pts2d]=ps2pmstruct(s,useSemiLocal)
%PS2PMSTRUCT Convert PhotoScan structure to PhotoModeler structure.
%
%   PROB=PS2PMSTRUCT(S) converts the PhotoScan struct S, as returned
%   by LOADPSZ, to a PhotoModeler struct PROB, as loaded by
%   LOADPM. All EO, OP parameters are global. Use PS2PMSTRUCT(S,TRUE)  
%   to use semilocal parameters instead (translation and scaling
%   from global, but no rotation).
%
%   [PROB,PTS3D,PTS2D]=PS2PMSTRUCT(S) furthermore returns 3D and 2D
%   points.
%
%   Zero-based object point ids in S will be converted to one-based
%   ids in PROB.
%
%   Camera positions, object and control points are set to global
%   coordinates. Object point uncertainties are set to
%   unknown. Mark points are given 1 pixel standard deviation.
%
%See also: LOADPSZ, LOADPM.

if nargin<2, useSemiLocal=false; end

% Create a fake job header with default camera.
imSz=s.camera.imSz(:);
defCam=[s.camera.focal;s.camera.pp(:);s.camera.sensorFormat(:);zeros(5,1)];

job=struct('title','Photoscan import','defCam',defCam,'defCamStd',zeros(size(defCam)),'imSz',imSz);

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
                     'imSz',imSz),size(CC,2),1);
for i=1:length(images)
    images(i).imName=s.imNames{i};
    images(i).outer=[CC(:,i);angPM(:,i)]';
end

% Copy enabled control points.
ctrlPts=pos.ctrlPts(pos.ctrlPtsEnabled,:);

% Copy global object points. Set posterior uncertainty to unknown.
objPts=[ctrlPts;[pos.objPts,nan(size(pos.objPts,1),3)]];

% Copy mark points and set std.
ctrlStd=s.defStd.projections;
objStd=s.defStd.tiePoints;
markPts=[s.markPts.ctrl,repmat(ctrlStd,size(s.markPts.ctrl,1),2);
         s.markPts.obj,repmat(objStd,size(s.markPts.obj,1),2)];

% Sort by image, then id.
markPts=msort(markPts);

if ~isempty(markPts)
    % Convert image indices to zero-based for PhotoModeler compatibility.
    markPts(:,1)=markPts(:,1)-1;
end

% Construct PM structure.
prob=struct('job',job,'images',images,'ctrlPts',ctrlPts,'objPts',objPts,...
             'markPts',markPts);

EO=[CC;ang;zeros(1,size(CC,2))];
EOstd=nan(size(EO));
    
pmReport=struct('EO',EO,'EOstd',EOstd);

% Create tables of 3D and 2D info.
id3d=zeros(1,0);
if ~isempty(pos.ctrlPts)
    id3d=[id3d,pos.ctrlPts(:,1)'];
end
if ~isempty(pos.objPts)
    id3d=[id3d,pos.objPts(:,1)'];
end
name=repmat({''},size(id3d));
name(s.DBATCPid(s.raw.ctrlPts(:,1)))=s.raw.ctrlPtsLabels;
pos3d=zeros(3,0);
if size(pos.ctrlPts,2)>=4
    pos3d=[pos3d,pos.ctrlPts(:,2:4)'];
end
if size(pos.objPts,2)>=4
    pos3d=[pos3d,pos.objPts(:,2:4)'];
end
std3d=zeros(3,0);
if size(pos.ctrlPts,2)>=7
    std3d=[std3d,pos.ctrlPts(:,5:7)'];
else
    std3d=[std3d,nan(size(pos.ctrlPts,1),3)'];
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
