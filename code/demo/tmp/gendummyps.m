ctrlPts=[1001 0 0 0
         1002 1 0 0
         1003 0 1 0
         1004 1 1 0
         1005 0 0 1
         1006 1 0 1
         1007 0 1 1
         1008 1 1 1
         1009 0.5 0.5 0.5
         1010 1.25 1.25 1.25];
ctrlPts=[ctrlPts,1e-6*ones(size(ctrlPts,1),3)];

nOtherPts=100;
otherPts=[1:nOtherPts;randn(3,nOtherPts);zeros(3,nOtherPts)]';

CC=[kron(ones(1,3),[-0.5,0,0.5])
    kron([-0.5,0,0.5],ones(1,3))
    5*ones(1,9)];
nCams=size(CC,2);

isCtrl=ismember(ctrlPts(:,1),1001:1010);

f=10;
imSz=[1000,800];
wh=[10,8];
defCam=[20,wh/2,wh,zeros(1,5)]';
p2=struct('job',struct('imSz',imSz,'defCam',defCam));

imName=fullfile(getenv('HOME'),'winshare','pm','square','white1000x800.png');
p2.images=repmat(struct('imName',imName,...
                        'outer',[CC(:,1);zeros(3,1)]'),1,nCams);
for i=1:length(p2.images)
    p2.images(i).outer(1:3)=CC(:,i)';
end

p2.ctrlPts=ctrlPts(isCtrl,:);
p2.objPts=[ctrlPts;otherPts];

K=[f,0,wh(1)/2;0,f,wh(2)/2;0,0,1];

K=diag([imSz./wh,1])*K;

P=zeros(3,4,nCams);
for i=1:nCams
    R=pm_eulerrotmat(zeros(1,3));
    P(:,:,i)=R*[eye(3),-CC(:,i)];
end

markPts=[];
for i=1:nCams
    newPts=euclidean(K*P(:,:,i)*homogeneous(p2.objPts(:,2:4)'));
    markPts=[markPts;repmat(i-1,size(p2.objPts(:,1))),p2.objPts(:,1),newPts'];
end

p2.markPts=markPts;

pmtops(p2,'/tmp/dummyps.psz')