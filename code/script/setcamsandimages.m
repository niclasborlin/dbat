function s=setcamsandimages(s,cams,ims)
%SETCAMSANDIMAGES Set camera and image info in a DBAT structure
%
%   S=SETIMAGES(S,CAMS,IMS) sets the camera and image info in the DBAT
%   structure S from the DBAT camera structure array CAMS and the DBAT
%   image struct IMS. S should have been initialized to its proper
%   size.

nIms=length(ims.id);
nK=cams(1).nK;
nP=cams(1).nP;
% Camera that was used for each image
camNo=ims.cam;

s.IO.model.camUnit=cams(1).unit;

s=setcamlenscoeff(s,cams(1).nK,cams(1).nP);
s=setcammodel(s,cams(1).model);

s.IO.sensor.ssSize(:)=cat(1,cams(ims.cam).sensor)';
s.IO.sensor.imSize(:)=cat(1,cams(ims.cam).image)';
s.IO.sensor.pxSize(:)=s.IO.sensor.ssSize./s.IO.sensor.imSize;
s.IO.sensor.samePxSize=isscalar(unique(s.IO.sensor.pxSize));

s.EO.name=ims.path;
s.EO.label=ims.name;
s.EO.cam=ims.id;
s.EO.id=ims.id;
