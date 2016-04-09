function pmtops(prob,psFile)
%PMTOPS Convert Photomodeler prob structure to Photoscan file structure.
%
%   PMTOPS(PROB,PSFILE) converts the structure PROB loaded from a
%   Photomodeler file to Photoscan. The PSFILE file name should point to
%   the main doc.xml file. Both structures must exist and have the same
%   number of images.

% $Id$

if 0
% Load existing doc.xml.
s=xml2struct2(psFile);

% Sensor info
s.document.chunks.chunk.sensors.sensor.Attributes.label='Photomodeler camera';
s.document.chunks.chunk.sensors.sensor.resolution.Attributes.width=...
    sprintf('%d',prob.job.imSz(1));
s.document.chunks.chunk.sensors.sensor.resolution.Attributes.height=...
    sprintf('%d',prob.job.imSz(2));

strs={'pixel_width','pixel_height','focal_length'};
vals=[prob.job.defCam(4:5)./prob.job.imSz(1:2);prob.job.defCam(1)];

psStrs=cellfun(@(x)x.Attributes.name,s.document.chunks.chunk.sensors.sensor.property,'uniformoutput',false);

[~,ia,ib]=intersect(psStrs,strs);

for i=1:length(ia)
    s.document.chunks.chunk.sensors.sensor.property{ia(i)}.Attributes.value=...
        sprintf('%.16f',vals(ib(i)));
end

j=find(ismember('fixed',psStrs));
s.document.chunks.chunk.sensors.sensor.property{j}.Attributes.value='true';

s.document.chunks.chunk.sensors.sensor.calibration.Attributes.class='initial';


end

noText=char(zeros(1,0));

camName='Photomodeler imported project camera';

%% TODO: Continue backwards in this file.

eye3x3='1 0 0 0 1 0 0 0 1';

minPos=min(prob.objPts(:,2:4));
maxPos=max(prob.objPts(:,2:4));
center=struct('Text',sprintf('%.16e ',(minPos+maxPos)/2));
sz=struct('Text',sprintf('%.16e ',maxPos-minPos));
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
for i=1:length(projections)
    camId=i-1;
    projections{i}=struct('Text',noText,...
                          'Attributes',...
                          struct('camera_id',sprintf('%d',camId),...
                                 'path',sprintf('projections%d.ply',camId)));
end

point_cloud=struct('tracks',tracks,'points',points,...
                   'projections',{projections});

%% TODO: Write tracks0.ply file.
%% TODO: Write points0.ply file.
%% TODO: Write projectionsXX.ply file.

marker=cell(1,size(prob.ctrlPts,1));

for i=1:length(marker)
    id=prob.ctrlPts(i,1);
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
    mAttr=struct('marker_id',sprintf('%d',id));
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

marker=cell(1,size(prob.ctrlPts,1));

for i=1:size(prob.ctrlPts,1)
    cpAttrib=struct('id',sprintf('%d',prob.ctrlPts(i,1)),...
                    'label',sprintf('CP%d',prob.ctrlPts(i,1)));
    refText=char(zeros(1,0));
    refAttr=struct('enabled','true',...
                   'x',sprintf('%.16g',prob.ctrlPts(i,2)),...
                   'y',sprintf('%.16g',prob.ctrlPts(i,3)),...
                   'z',sprintf('%.16g',prob.ctrlPts(i,4)),...
                   'sx',sprintf('%.16g',max(prob.ctrlPts(i,5),0)),...
                   'sy',sprintf('%.16g',max(prob.ctrlPts(i,6),0)),...
                   'sz',sprintf('%.16g',max(prob.ctrlPts(i,7),0)));
    m=struct('reference',struct('Text',refText,'Attributes',refAttr),...
             'Attributes',cpAttrib);
    marker{i}=m;
end

markersAttr=struct('next_id',sprintf('%d',max(prob.ctrlPts(:,1))+1));
markers=struct('marker',{marker},'Attributes',markersAttr);

camera=cell(1,length(prob.images));

for i=1:length(prob.images)
    [~,n,e]=fileparts(strrep(prob.images(i).imName,'\','/'));
    attrib=struct('enable','true','id',sprintf('%d',i-1),...
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
    invP=inv([P;0,0,0,1]);
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

struct2xml(s,psFile);

