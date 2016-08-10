function pmtops(prob,psFile,newImDir)
%PMTOPS Convert Photomodeler prob structure to Photoscan file structure.
%
%   PMTOPS(PROB,PSFILE) converts the structure PROB loaded from a
%   Photomodeler file to Photoscan. The PSFILE file name should be
%   the path to the main doc.xml file.

% $Id$

psDir=fileparts(psFile);

psUnpackedDir=fullfile(psDir,'unpacked');

if ~exist(psUnpackedDir)
    mkdir(psUnpackedDir)
end

delete(fullfile(psUnpackedDir,'*.ply'));

noText=char(zeros(1,0));

camName='Photomodeler imported project camera';

eye3x3='1 0 0 0 1 0 0 0 1';

minPos=min(prob.objPts(:,2:4));
maxPos=max(prob.objPts(:,2:4));
center=struct('Text',sprintf('%.16e ',(minPos+maxPos)/2));
sz=struct('Text',sprintf('%.16e ',max(maxPos-minPos,0.01)));
R=struct('Text',eye3x3);
region=struct('center',center,'size',sz,'R',R);

rotation=struct('Text',eye3x3);
translation=struct('Text','0 0 0');
scale=struct('Text','1');
transform=struct('rotation',rotation,'translation',translation,'scale',scale);

tracksFileName='tracks0.ply';

tracks=struct('Text',noText,'Attributes',struct('path',tracksFileName));

pointsFileName='points0.ply';

points=struct('Text',noText,'Attributes',struct('path',pointsFileName));

projections=cell(1,length(prob.images));
projFileNames=cell(1,length(prob.images));
for i=1:length(projections)
    camId=i-1;
    projFileNames{i}=sprintf('projections%d.ply',camId);
    projections{i}=struct('Text',noText,...
                          'Attributes',...
                          struct('camera_id',sprintf('%d',camId),...
                                 'path',projFileNames{i}));
end
projFullFileNames=fullfile(psUnpackedDir,projFileNames);

point_cloud=struct('tracks',tracks,'points',points,...
                   'projections',{projections});

% Write all non-ctrl object points to the point0.ply file.
pointsFullFileName=fullfile(psUnpackedDir,pointsFileName);
isObj=~ismember(prob.objPts(:,1),prob.ctrlPts(:,1));
objPts=prob.objPts(isObj,:);
[~,i]=sort(objPts(:,1));
objPts=objPts(i,:);
if ~isempty(objPts)
    vertex=struct('x',objPts(:,2),'y',objPts(:,3),'z',objPts(:,4),'id',objPts(:,1));
    d=struct('vertex',vertex);
    ply_write(d,pointsFullFileName,'binary_little_endian');
end

% Write all mark points except those of control points to
% projectionsXX.ply files.
for i=1:length(projections)
    camId=i-1;
    j=prob.markPts(:,1)==camId & ~ismember(prob.markPts(:,2),prob.ctrlPts(:,1));
    markPts=prob.markPts(j,:);
    [~,k]=sort(markPts(:,2));
    markPts=markPts(k,:);
    if ~isempty(markPts)
        vertex=struct('x',markPts(:,3),'y',markPts(:,4),...
                      'size',4*ones(size(markPts,1),1),'id',markPts(:,2));
        d=struct('vertex',vertex);
        ply_write(d,projFullFileNames{i},'binary_little_endian');
    end
end

% Write dummy tracks0.ply file.
tracksFullFileName=fullfile(psUnpackedDir,tracksFileName);
maxId=max(setdiff(prob.markPts(:,2),prob.ctrlPts(:,1)));
if ~isempty(maxId)
    dummy=zeros(maxId+1,'uint8');
    d=struct('vertex',struct('red',dummy,'green',dummy,'blue',dummy));
    ply_write(d,tracksFullFileName,'binary_little_endian');
end

ctrlPts=prob.ctrlPts;
[~,i]=sort(ctrlPts(:,1));
ctrlPts=ctrlPts(i,:);

marker=cell(1,size(ctrlPts,1));

for i=1:length(marker)
    id=ctrlPts(i,1);
    measured=find(prob.markPts(:,2)==id);
    location=cell(1,length(measured));
    for j=1:length(location)
        r=measured(j);
        location{j}=struct('Text',noText,...
                           'Attributes',...
                           struct('camera_id',sprintf('%d',...
                                                      prob.markPts(r,1)),...
                                  'pinned','true',...
                                  'x',sprintf('%.8g',prob.markPts(r,3)),...
                                  'y',sprintf('%.8g',prob.markPts(r,4))));
    end
    mAttr=struct('marker_id',sprintf('%d',id-min(ctrlPts(:,1))));
    marker{i}=struct('location',{location},'Attributes',mAttr);
end

markers=struct('marker',{marker});

camera=cell(1,length(prob.images));

for i=1:length(prob.images)
    thumbnail=struct('Text',noText','Attributes',...
                     struct('path',sprintf('thumb%d.jpg',i-1)));
    propFocal=struct('Text',noText,'Attributes',...
                     struct('name','Exif/FocalLength',...
                            'value',sprintf('%.4f',prob.job.defCam(1))));
    propModel=struct('Text',noText,'Attributes',...
                     struct('name','Exif/Model',...
                            'value',camName));
    propHeight=struct('Text',noText,'Attributes',...
                     struct('name','File/ImageHeight',...
                            'value',sprintf('%d',prob.job.imSz(2))));
    propWidth=struct('Text',noText,'Attributes',...
                     struct('name','File/ImageWidth',...
                            'value',sprintf('%d',prob.job.imSz(1))));
    property={propFocal,propModel,propHeight,propWidth};
    
    meta=struct('property',{property});
    photoAttr=struct('path',strrep(prob.images(i).imName,'\','/'));
    photo=struct('meta',{meta},'Attributes',photoAttr);
    cameraAttr=struct('camera_id',sprintf('%d',i-1));
    camera{i}=struct('photo',photo,'thumbnail',thumbnail,'Attributes',cameraAttr);
end

cameras=struct('camera',{camera});

frameAttr=struct('id','0');
frame=struct('cameras',cameras,'markers',markers,'point_cloud',point_cloud,...
             'Attributes',frameAttr);

framesAttr=struct('next_id','1');
frames=struct('frame',frame,'Attributes',framesAttr);

ctrlPts=prob.ctrlPts;
[~,i]=sort(ctrlPts(:,1));
ctrlPts=ctrlPts(i,:);

marker=cell(1,size(ctrlPts,1));

for i=1:size(ctrlPts,1)
    cpAttrib=struct('id',sprintf('%d',ctrlPts(i,1)-min(ctrlPts(:,1))),...
                    'label',sprintf('CP%d',ctrlPts(i,1)));
    refText=char(zeros(1,0));
    refAttr=struct('enabled','true',...
                   'x',sprintf('%.16g',ctrlPts(i,2)),...
                   'y',sprintf('%.16g',ctrlPts(i,3)),...
                   'z',sprintf('%.16g',ctrlPts(i,4)),...
                   'sx',sprintf('%.16g',max(ctrlPts(i,5),0)),...
                   'sy',sprintf('%.16g',max(ctrlPts(i,6),0)),...
                   'sz',sprintf('%.16g',max(ctrlPts(i,7),0)));
    m=struct('reference',struct('Text',refText,'Attributes',refAttr),...
             'Attributes',cpAttrib);
    marker{i}=m;
end

markersAttr=struct('next_id',sprintf('%d',max(ctrlPts(:,1))+1));
markers=struct('marker',{marker},'Attributes',markersAttr);

camera=cell(1,length(prob.images));

for i=1:length(prob.images)
    [~,n,e]=fileparts(strrep(prob.images(i).imName,'\','/'));
    attrib=struct('enabled','true','id',sprintf('%d',i-1),...
                  'label',[n,e],'sensor_id','0');

    orient=struct('Text','1');
    
    % Construct camera matrix.
    outer=prob.images(i).outer;
    EO=nan(6,1);
    
    EO(1:3)=outer(1:3);
    % Euler angles, stored as kappa, phi, omega in PM file.
    EO(4:6)=outer([6,5,4])/180*pi;

    RR=pm_eulerrotmat(EO(4:6));
    CC=EO(1:3);
    P=RR*[eye(3),-CC];
    invP=inv([P;0,0,0,1])*diag([1,-1,-1,1]);
    xForm=struct('Text',sprintf('%.16e ',invP'));
    
    c=struct('transform',xForm,'orientation',orient,'Attributes',attrib);
    camera{i}=c;
end
camerasAttr=struct('next_group_id','0',...
                   'next_id',sprintf('%d',length(prob.images)));
cameras=struct('camera',{camera},'Attributes',camerasAttr);

resolutionAttr=struct('height',sprintf('%d',prob.job.imSz(2)),...
                      'width',sprintf('%d',prob.job.imSz(1)));
resolution=struct('Text',noText,'Attributes',resolutionAttr);
fx=struct('Text',sprintf('%.16f',prob.job.defCam(1)*prob.job.imSz(1)/prob.job.defCam(4)));
fy=struct('Text',sprintf('%.16f',prob.job.defCam(1)*prob.job.imSz(2)/prob.job.defCam(5)));
cx=struct('Text',sprintf('%.16f',prob.job.defCam(2)*prob.job.imSz(1)/prob.job.defCam(4)));
cy=struct('Text',sprintf('%.16f',prob.job.defCam(3)*prob.job.imSz(2)/prob.job.defCam(5)));
calibrationAttr=struct('class','initial','type','frame');
calibration=struct('resolution',resolution,'fx',fx,'fy',fy,'cx',cx,'cy',cy,...
                   'Attributes',calibrationAttr);

colors={'Red','Green','Blue'};

band=cell(size(colors));
for i=1:length(band)
    band{i}=struct('Text',noText,'Attributes',struct('label',colors{i}));
end

bands=struct('band',{band});

property=cell(1,4);
property{1}=struct('Text',noText,'Attributes',...
                   struct('name','pixel_width',...
                          'value',...
                          sprintf('%.16g',...
                                  prob.job.defCam(4)/prob.job.imSz(1))));
property{2}=struct('Text',noText,'Attributes',...
                   struct('name','pixel_height',...
                          'value',...
                          sprintf('%.16g',...
                                  prob.job.defCam(5)/prob.job.imSz(2))));
property{3}=struct('Text',noText,'Attributes',...
                   struct('name','focal_length',...
                          'value',...
                          sprintf('%.16g',...
                                  prob.job.defCam(1))));
property{4}=struct('Text',noText,'Attributes',...
                   struct('name','fixed',...
                          'value','true'));

resolutionAttr=struct('height',sprintf('%d',prob.job.imSz(2)),...
                      'width',sprintf('%d',prob.job.imSz(1)));
resolution=struct('Text',noText,'Attributes',resolutionAttr);

sensorAttr=struct('id','0','label',camName,...
                  'type','frame');
sensor=struct('resolution',resolution,'property',{property},'bands',bands,...
              'calibration',calibration,'Attributes',sensorAttr);

sensorsAttr=struct('next_id','1');
sensors=struct('sensor',sensor,'Attributes',sensorsAttr);

chunkAttr=struct('enabled','true','id','0','label','Chunk 1 (imported)');
chunk=struct('sensors',sensors,'cameras',cameras,'markers',markers,...
             'frames',frames,'transform',transform,'region',region,...
             'Attributes',chunkAttr);

chunksAttr=struct('next_id','1');
chunks=struct('chunk',chunk,'Attributes',chunksAttr);

docAttr=struct('version','1.2.0');
doc=struct('chunks',chunks,'Attributes',docAttr);

s=struct('document',doc);

docFileName='doc.xml';
docFullFile=fullfile(psUnpackedDir,docFileName);

struct2xml(s,docFullFile);

tmpFile=[tempname(psDir),'.zip'];

filesToZip={docFileName,tracksFileName,pointsFileName,projFileNames{:}};

for i=1:length(filesToZip)
    [~,~]=system(sprintf(['touch --date="1979-12-31 00:00:00 +0100" "%s"'],...
                         fullfile(psUnpackedDir,filesToZip{i})));
end

zip(tmpFile,filesToZip,psUnpackedDir);

movefile(tmpFile,psFile);
