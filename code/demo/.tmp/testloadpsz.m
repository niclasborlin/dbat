psStub1='camcal7dout2';
psFile1=fullfile(fileparts(mfilename('fullpath')),'data', ...
                 'weighted','ps',psStub1,[psStub1,'.psz']);

s1=loadpsz2(psFile1);

imNo=length(s1.imNames);

imshow(s1.imNames{imNo})

xo=euclidean(s1.K*s1.global.P(:,:,end)*homogeneous(s1.global.objPts(:,2:4)'));

line(xo(1,:),xo(2,:),'marker','x','linestyle','none','color','g')

xc=euclidean(s1.K*s1.global.P(:,:,end)*homogeneous(s1.global.ctrlPts(:,2:4)'));

line(xc(1,:),xc(2,:),'marker','o','linestyle','none','color','g')

ixc=s1.markPts.ctrl(:,1)==imNo-1;

line(s1.markPts.ctrl(ixc,3),s1.markPts.ctrl(ixc,4),'marker','o','color','r',...
     'linestyle','none');

ixo=s1.markPts.obj(:,1)==imNo-1;

line(s1.markPts.obj(ixo,3),s1.markPts.obj(ixo,4),'marker','x','color','r',...
     'linestyle','none');

% Project control points into each image.
proj=nan(size(s1.markPts.ctrl,1),2);
for i=1:size(proj,1)
    % For each measured marker.
    imNo=find(s1.markPts.ctrl(i,1)+1==s1.cameraIds);
    ctrlId=s1.markPts.ctrl(i,2);
    cIx=find(s1.global.ctrlPts(:,1)==ctrlId);
    proj(i,:)=euclidean(s1.K*s1.global.P(:,:,imNo)*homogeneous(s1.global.ctrlPts(cIx,2:4)'))';
end

% Compute residuals.
res=proj-s1.markPts.ctrl(:,3:4);

rr=cell(1,size(s1.global.ctrlPts,1));
% Group residuals per point.
for i=1:size(s1.global.ctrlPts,1)
    rr{i}=res(s1.markPts.ctrl(:,2)==s1.global.ctrlPts(i,1),:);
end

asdf
fName='/Home/staff/niclas/photoscan/cptest2.psz';


[p,doc,pts]=loadpsz(fName);

TT=zeros(4,4,length(doc.document.chunks.chunk.cameras.camera));

for i=1:size(T,3)
    TT(:,:,i)=reshape(sscanf(doc.document.chunks.chunk.cameras.camera{i}.transform.Text,'%g'),4,4)';
end

% Control points in 
% doc.document.chunks.chunk.markers.marker{i}.reference

t=sscanf(doc.document.chunks.chunk.transform.translation.Text,'%g');
s=sscanf(doc.document.chunks.chunk.transform.scale.Text,'%g');
R=reshape(sscanf(doc.document.chunks.chunk.transform.rotation.Text,'%g'),3,3)';

T0=[eye(3),t;0,0,0,1];
S0=diag([s,s,s,1]);
R0=blkdiag(R,1);

X0=T0*S0*R0;

I34=eye(3,4);

% Reconstruct EO
rEO=zeros(size(EO));
for i=1:size(TT,3)
    P=I34*inv(TT(:,:,i))*inv(X0);
    rEO(1:3,i)=euclidean(null(P));
    M=P(1:3,1:3);
    [U,S,V]=svd(M);
    Mh=diag([1,-1,-1])*U*diag([1,1,det(U*V')])*V';
    rEO(4:6,i)=derotmat3d(Mh)';
end
