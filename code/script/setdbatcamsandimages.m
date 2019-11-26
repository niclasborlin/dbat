function s=setdbatcamsandimages(s,cams,ims)
%SETDBATCAMSANDIMAGES Set camera and image info in a DBAT structure
%
%   S=SETDBATCAMSANDIMAGES(S,CAMS,IMS) sets the camera and image info
%   in the DBAT structure S from the DBATCamera cell array CAMS and
%   the DBAT image struct IMS. S should have been initialized to its
%   proper size.

% Camera that was used for each image
camNo=ims.cam;

s.IO.model.camUnit=cams{1}.Unit;

s.IO.cam=cams;

s.IO.sensor.ssSize(:)=cell2mat(cellfun(@(x)x.SensorSize', ...
                                       cams(ims.cam),'UniformOutput',false));
s.IO.sensor.imSize(:)=cell2mat(cellfun(@(x)x.ImageSize', ...
                                       cams(ims.cam),'UniformOutput',false));
s.IO.sensor.pxSize(:)=cell2mat(cellfun(@(x)PixelSize(x)', ...
                                       cams(ims.cam),'UniformOutput',false));
s.IO.sensor.samePxSize=isscalar(unique(s.IO.sensor.pxSize));

s=setcamlenscoeff(s,nK(cams{1}),nP(cams{1}));
s=setcammodel(s,cams{1}.Model);

s.IO.struct.block=repmat(camNo,size(s.IO.val,1),1);

s.EO.name=ims.path;
s.EO.label=ims.name;
s.EO.cam=ims.id;
s.EO.id=ims.id;

% Extract camera values from structure.

cc=cellfun(@(x)x.CameraConstant,cams);
pp=cell2mat(cellfun(@(x)x.PrincipalPoint',cams,'UniformOutput',false));
skew=cellfun(@(x)x.Skew,cams);
aspectDiff=cellfun(@(x)AspectDiff(x),cams);
K=cell2mat(cellfun(@(x)x.K',cams,'UniformOutput',false));
P=cell2mat(cellfun(@(x)x.P',cams,'UniformOutput',false));

s=setcamvals(s,'prior','cc',cc(:,camNo),'pp',pp(:,camNo), ...
               'as', aspectDiff(:,camNo),'sk',skew(:,camNo), ...
               'K',K(:,camNo),'P',P(:,camNo));

s=parseblockvariant(s);
s=buildparamtypes(s);
