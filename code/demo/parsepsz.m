clc
fName=['/home/niclas/dbat/code/demo/data/weighted/ps/sxb/weighted/unpacked/doc.xml'];
s=xml2struct2(fName);
nCams=length(s.document.chunks.chunk.cameras.camera);
cams=nan(3,4,nCams);
CC=nan(3,nCams);
for i=1:nCams
    t=s.document.chunks.chunk.cameras.camera{i}.transform.Text;
    P=reshape(sscanf(t,'%g '),4,[])';
    cams(:,:,i)=P(1:3,:);
    CC(:,i)=euclidean(null(cams(:,:,i)));
end
tv=sscanf(s.document.chunks.chunk.transform.translation.Text,'%g ');
T0=[eye(3),tv(:);0,0,0,1];
rv=sscanf(s.document.chunks.chunk.transform.rotation.Text,'%g ');
R0=blkdiag(reshape(rv,3,[])',1);
sv=sscanf(s.document.chunks.chunk.transform.scale.Text,'%g ');
S0=diag([repmat(sv,1,3),1]);

TSR=T0*S0*R0;
TRS=T0*R0*S0;
STR=S0*T0*R0;
SRT=S0*R0*T0;
RST=R0*S0*T0;
RTS=R0*T0*S0;

disp('TSR')
euclidean(TSR*homogeneous(CC))
disp('TRS')
euclidean(TRS*homogeneous(CC))
disp('STR')
euclidean(STR*homogeneous(CC))
disp('SRT')
euclidean(SRT*homogeneous(CC))
disp('RST')
euclidean(RST*homogeneous(CC))
disp('RTS')
euclidean(RTS*homogeneous(CC))

pts3d=ply_read