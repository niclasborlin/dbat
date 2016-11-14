function [prob,pmReport,pts3d,pts2d]=ps2pmstruct(s)
%PS2PMSTRUCT Convert PhotoScan structure to PhotoModeler structure.
%
%   PROB=PS2PMSTRUCT(S) converts the PhotoScan struct S, as returned
%   by LOADPSZ, to a PhotoModeler struct PROB, as loaded by LOADPM.
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

% Create a fake job header with default camera.
imSz=s.camera.imSz(:);
defCam=[s.camera.focal;s.camera.pp(:);s.camera.sensorFormat(:);zeros(5,1)];

job=struct('title','Photoscan import','defCam',defCam,'defCamStd',zeros(size(defCam)),'imSz',imSz);

% Create EO parameters.

% Copy global camera positions.
CC=s.global.CC;

% Convert rotation matrices to omega-phi-kappa.
ang=nan(size(CC));
for i=1:size(ang,2)
    RR=s.global.R(:,:,i);
    ang(:,i)=derotmat3d(RR)';
end
angPM=ang([3,2,1],:)*180/pi;

images=repmat(struct('imName','','outer',nan(1,6),'outerStd',zeros(1,6),...
                     'imSz',imSz),size(CC,2),1);
for i=1:length(images)
    images(i).imName=s.imNames{i};
    images(i).outer=[CC(:,i);angPM(:,i)]';
end

% Copy global control points.
ctrlPts=s.global.ctrlPts;

% Copy global object points. Set posterior uncertainty to unknown.
objPts=[ctrlPts;[s.global.objPts,nan(size(s.global.objPts,1),3)]];

% Copy mark points and set std.
ctrlStd=s.defStd.projections;
objStd=s.defStd.tiePoints;
markPts=[s.markPts.ctrl,repmat(ctrlStd,size(s.markPts.ctrl,1),2);
         s.markPts.obj,repmat(objStd,size(s.markPts.obj,1),2)];

% Sort by image, then id.
markPts=msort(markPts);

% Convert ids to one-based.
if ~isempty(ctrlPts)
    ctrlPts(:,1)=ctrlPts(:,1)+1;
end
if ~isempty(objPts)
    objPts(:,1)=objPts(:,1)+1;
end
if ~isempty(markPts)
    % Leave image indices zero-based for PhotoModeler compatibility.
    markPts(:,2)=markPts(:,2)+1;
end

% Construct PM structure.
prob=struct('job',job,'images',images,'ctrlPts',ctrlPts,'objPts',objPts,...
             'markPts',markPts);

EO=[CC;ang;zeros(1,size(CC,2))];
EOstd=nan(size(EO));
    
pmReport=struct('EO',EO,'EOstd',EOstd);

% Create tables of 3D and 2D info.
id3d=zeros(1,0);
if ~isempty(s.global.ctrlPts)
    id3d=[id3d,s.global.ctrlPts(:,1)'];
end
if ~isempty(s.global.objPts)
    id3d=[id3d,s.global.objPts(:,1)'];
end
name=repmat({''},size(id3d));
pos3d=zeros(3,0);
if size(s.global.ctrlPts,2)>=4
    pos3d=[pos3d,s.global.ctrlPts(:,2:4)'];
end
if size(s.global.objPts,2)>=4
    pos3d=[pos3d,s.global.objPts(:,2:4)'];
end
std3d=zeros(3,0);
if size(s.global.ctrlPts,2)>=7
    std3d=[std3d,s.global.ctrlPts(:,5:7)'];
else
    std3d=[std3d,nan(size(s.global.ctrlPts,1),3)'];
end
if size(s.global.objPts,2)>=7
    std3d=[std3d,s.global.objPts(:,5:7)'];
else
    std3d=[std3d,nan(size(s.global.objPts,1),3)'];
end
% Convert to one-based.
id3d=id3d+1;
% Create a 'raw' visibility map, i.e. for all object points.
rawVis=sparse(s.markPts.all(:,2)+1,s.markPts.all(:,1)+1,1);
% Keep only the object points that were actually computed.
vis=rawVis(id3d,:);

% Sort 3d pts by id.
[~,i]=sort(id3d);
id3d=id3d(i);
name=name(i);
vis=vis(i,:);
pos3d=pos3d(:,i);
std3d=std3d(:,i);
pts3d=struct('id',id3d,'name',{name},'vis',logical(vis),'pos',pos3d,'std',std3d);

% Sort 2d pts by image first, then id.
[~,i]=sortrows(s.markPts.all,[1,2]);
s.markPts.all=s.markPts.all(i,:);

% Only keep 2d measurements that corresponds to 3D points.
keep=ismember(s.markPts.all(:,2)+1,pts3d.id);

id2d=s.markPts.all(keep,2)'+1;
imNo=s.markPts.all(keep,1)'+1;
pos2d=s.markPts.all(keep,3:4)';
res2d=nan(size(pos2d));

pts2d=struct('id',id2d,'imNo',imNo,'pos',pos2d,'res',res2d,'kept',keep);
