clc
fName=['/home/niclas/dbat/code/demo/data/weighted/ps/sxb/test2/unpacked/doc.xml'];
s=xml2struct2(fName);
nCams=length(s.document.chunks.chunk.cameras.camera);
cams=nan(3,4,nCams);
CC=nan(3,nCams);
for i=1:nCams
    t=s.document.chunks.chunk.cameras.camera{i}.transform.Text;
    P=inv(reshape(sscanf(t,'%g '),4,[])');
    cams(:,:,i)=P(1:3,:);
    CC(:,i)=euclidean(null(cams(:,:,i)));
end
tv=sscanf(s.document.chunks.chunk.transform.translation.Text,'%g ');
T0=[eye(3),tv(:);0,0,0,1];
rv=sscanf(s.document.chunks.chunk.transform.rotation.Text,'%g ');
R0=blkdiag(reshape(rv,3,[])',1);
sv=sscanf(s.document.chunks.chunk.transform.scale.Text,'%g ');
S0=diag([repmat(sv,1,3),1]);

cv=sscanf(s.document.chunks.chunk.region.center.Text,'%g ');
T1=[eye(3),cv(:);0,0,0,1];
rv=sscanf(s.document.chunks.chunk.region.R.Text,'%g ');
R1=blkdiag(reshape(rv,3,[])',1);

cal=s.document.chunks.chunk.sensors.sensor.calibration;
fx=sscanf(cal.fx.Text,'%g');
fy=sscanf(cal.fy.Text,'%g');
cx=sscanf(cal.cx.Text,'%g');
cy=sscanf(cal.cy.Text,'%g');
K=[fx,0,cx;0,fy,cy;0,0,1];

CP=nan(3,length(s.document.chunks.chunk.markers.marker));
for i=1:size(CP,2);
    m=s.document.chunks.chunk.markers.marker{i};
    x=sscanf(m.reference.Attributes.x,'%g');
    y=sscanf(m.reference.Attributes.y,'%g');
    z=sscanf(m.reference.Attributes.z,'%g');
    CP(:,i)=[x,y,z]';
end

TSR=T0*S0*R0;

euclidean(TSR*homogeneous(CC))

xy1=euclidean(K*cams(:,:,1)*TSR*homogeneous(CP(:,1)))

m1=s.document.chunks.chunk.frames.frame.markers.marker{1}.location.Attributes;
mark1=[sscanf(m1.x,'%g');sscanf(m1.y,'%g')]

res=xy1-mark1
