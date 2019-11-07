function s=setdbatcamsandimages(s,cams,ims)
%SETDBATCAMSANDIMAGES Set camera and image info in a DBAT structure
%
%   S=SETDBATCAMSANDIMAGES(S,CAMS,IMS) sets the camera and image info
%   in the DBAT structure S from the DBAT camera structure array CAMS
%   and the DBAT image struct IMS. S should have been initialized to
%   its proper size.

nIms=length(ims.id);
nK=cams(1).nK;
nP=cams(1).nP;
% Camera that was used for each image
camNo=ims.cam;

s.IO.model.camUnit=cams(1).unit;

s.IO.sensor.ssSize(:)=cat(1,cams(ims.cam).sensor)';
s.IO.sensor.imSize(:)=cat(1,cams(ims.cam).image)';
s.IO.sensor.pxSize(:)=s.IO.sensor.ssSize./s.IO.sensor.imSize;
s.IO.sensor.samePxSize=isscalar(unique(s.IO.sensor.pxSize));

s=setcamlenscoeff(s,cams(1).nK,cams(1).nP);
s=setcammodel(s,cams(1).model);

s.IO.struct.block=repmat(camNo,size(s.IO.val,1),1);

s.EO.name=ims.path;
s.EO.label=ims.name;
s.EO.cam=ims.id;
s.EO.id=ims.id;

% Extract camera values from structure.

cc=cat(1,cams.cc)';
pp=cat(1,cams.pp)';
skew=cat(1,cams.skew)';
aspectDiff=cat(1,cams.aspectDiff)';
K=cat(1,cams.K)';
P=cat(1,cams.P)';

s=setcamvals(s,'prior','cc',cc(:,camNo),'pp',pp(:,camNo), ...
               'as', aspectDiff(:,camNo),'sk',skew(:,camNo), ...
               'K',K(:,camNo),'P',P(:,camNo));

s=parseblockvariant(s);
s=buildparamtypes(s);
