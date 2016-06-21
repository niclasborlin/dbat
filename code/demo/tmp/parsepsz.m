psFile=fullfile(getenv('HOME'),...
                'dbat/code/demo/data/weighted/ps/camcal7dout/camcal7dout.psz');
unpackpsz(psFile);
psUnpackedDir=fullfile(fileparts(psFile),'unpacked');
fName=fullfile(psUnpackedDir,'doc.xml');
s=xml2struct2(fName);
camera=s.document.chunks.chunk.cameras.camera;
nCams=length(camera);
cams=nan(3,4,nCams);
CC=nan(3,nCams);
PP=nan(4,4,nCams);
for i=1:nCams
    t=camera{i}.transform.Text;
    P=reshape(sscanf(t,'%g '),4,[])';
    PP(:,:,i)=P;
    cams(:,:,i)=eye(3,4)*inv(P);
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
    id=sscanf(m.Attributes.id,'%d');
    x=sscanf(m.reference.Attributes.x,'%g');
    y=sscanf(m.reference.Attributes.y,'%g');
    z=sscanf(m.reference.Attributes.z,'%g');
    CP(:,id+1)=[x,y,z]';
end

TSR=T0*S0*R0;

CCglobal=euclidean(TSR*homogeneous(CC))

imNo=nCams;
mNo=1;

measuredMarkers=s.document.chunks.chunk.frames.frame.markers.marker;
markerId=cellfun(@(x)sscanf(x.Attributes.marker_id,'%d'),measuredMarkers);

xy=zeros(2,size(CP,2));

sortedMarkerId=sort(markerId);

for mNo=1:4
    mIx=find(markerId==sortedMarkerId(mNo));
    camId=cellfun(@(x)sscanf(x.Attributes.camera_id,'%d')+1,measuredMarkers{mIx}.location);
    lIx=find(camId==imNo);

    loc=measuredMarkers{mIx}.location{lIx};
    xy(:,mIx)=[sscanf(loc.Attributes.x,'%g');
               sscanf(loc.Attributes.y,'%g')];
end

T=P(:,:,end);
proj1=euclidean(K*eye(3,4)*inv(T)*inv(TSR)*homogeneous(CP))

xy

s=loadpsz2(psFile);

fig=tagfigure('image');
figure(fig);
imshow(s.imNames{imNo},'parent',gca(fig))
line(proj1(1,:),proj1(2,:),'color','r','marker','o','linestyle','none');
line(xy(1,:),xy(2,:),'color','g','marker','x','linestyle','none');

fig=tagfigure('3d');
figure(fig)
XYZ=euclidean(TSR*homogeneous(s.objPts(:,2:4)'));
plot3(XYZ(1,:),XYZ(2,:),XYZ(3,:),'x')
line(CCglobal(1,:),CCglobal(2,:),CCglobal(3,:),'marker','o','linestyle','none','color','r')
line(CP(1,:),CP(2,:),CP(3,:),'marker','^','linestyle','none','color','y')
axis equal

asfds
xy1=euclidean(K*cams(:,:,1)*TSR*homogeneous(CP(:,1)))
xy2=euclidean(K*cams(:,:,2)*TSR*homogeneous(CP(:,1)))

m1=measuredMarkers{1}.location{end}.Attributes;
mark1=[sscanf(m1.x,'%g');sscanf(m1.y,'%g')]

m2=s.document.chunks.chunk.frames.frame.markers.marker{1}.location{2}.Attributes;
mark2=[sscanf(m2.x,'%g');sscanf(m2.y,'%g')]

res1=euclidean(K*cams(:,:,1)*TSR*homogeneous(X))-mark1
res2=euclidean(K*cams(:,:,2)*TSR*homogeneous(X))-mark2
res=[res1,res2];
resVec=sqrt(sum(res.^2,1));
errorPix=sqrt(mean(resVec.^2))

pt3dName=fullfile(fDir,s.document.chunks.chunk.frames.frame.point_cloud.points.Attributes.path);

[~,~,d,~]=ply_read(pt3dName,'tri');

pts3d=[d.vertex.id,d.vertex.x,d.vertex.y,d.vertex.z];

tracksName=fullfile(fDir,s.document.chunks.chunk.frames.frame.point_cloud.tracks.Attributes.path);

[~,~,d,~]=ply_read(tracksName,'tri');

tracks=[d.vertex.red,d.vertex.green,d.vertex.blue];

f=s.document.chunks.chunk.frames.frame;
pts2d=zeros(0,6);

for i=1:length(f.point_cloud.projections)
    imNo=sscanf(f.point_cloud.projections{i}.Attributes.camera_id,'%d');
    pName=fullfile(fDir,f.point_cloud.projections{i}.Attributes.path);
    [~,~,d,~]=ply_read(pName,'tri');
    p=[repmat(imNo,size(d.vertex.id)),d.vertex.id,d.vertex.x,d.vertex.y....
       d.vertex.size,d.vertex.size];
    pts2d=[pts2d;[p,zeros(size(p,1),2)]];
end

imNames=cellfun(@(x)x.photo.Attributes.path,f.cameras.camera,...
                'uniformoutput',false)

vis=sparse(pts2d(:,2)+1,pts2d(:,1)+1,1);

imNo=1;
imshow(fullfile(fDir,'..',imNames{imNo}));

i=pts2d(:,1)==imNo-1;
imPts=pts2d(i,:);
id=imPts(:,2)+1;
xy=imPts(:,3:4);
sz=imPts(:,5);
%p=PTCircle2D(PTGaussian(xy'),sz',id');

n=sum(vis(id,:),2);
%
%
j=ismember(id,pts3d(:,1)+1);
line(imPts(j,3),imPts(j,4),'marker','o','linestyle','none','color','b');
line(imPts(~j,3),imPts(~j,4),'marker','o','linestyle','none','color','y');
%line(imPts(n>2,3),imPts(n>2,4),'marker','o','linestyle','none','color','b');
%hold on
%plot(p)
%hold off
