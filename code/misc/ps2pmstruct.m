function [prob,pmReport,pts3d,pts2d]=ps2pmstruct(s)
%PS2PMSTRUCT Convert PhotoScan structure to PhotoModeler structure.
%
%   PROB=PS2PMSTRUCT(S) converts the PhotoScan struct S, as returned
%   by LOADPSZ, to a PhotoModeler struct PROB, as loaded by LOADPM.
%
%   [PROB,PTS3D,PTS2D]=PS2PMSTRUCT(S) furthermore returns 3D and 2D
%   points.

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

% Copy mark points. Set uncertainty to 1 pixel.
markPts=[s.markPts.all,[0.1*ones(size(s.markPts.ctrl,1),2);1*ones(size(s.markPts.obj,1),2)]];

% Construct PM structure.
prob=struct('job',job,'images',images,'ctrlPts',ctrlPts,'objPts',objPts,...
             'markPts',markPts);

EO=[CC;ang;zeros(1,size(CC,2))];
EOstd=nan(size(EO));
    
pmReport=struct('EO',EO,'EOstd',EOstd);

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
% Create a 'raw' visibility map, i.e. for all object points.
rawVis=sparse(s.markPts.all(:,2),s.markPts.all(:,1)+1,1);
% Keep only the object points that were actually computed.
vis=rawVis(id3d,:);
pts3d=struct('id',id3d,'name',{name},'vis',vis,'pos',pos3d,'std',std3d);

% Only keep measurements that correspond to 3D points.
keep=ismember(s.markPts.all(:,2),pts3d.id);

id2d=s.markPts.all(keep,2)';
imNo=s.markPts.all(keep,1)'+1;
pos2d=s.markPts.all(keep,3:4)';
res2d=nan(size(pos2d));

pts2d=struct('id',id2d,'imNo',imNo,'pos',pos2d,'res',res2d);
