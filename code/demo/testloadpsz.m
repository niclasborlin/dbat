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
