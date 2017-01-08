function analyzepsz(s)

ims=unique(s.markPts.obj(:,1));

OPid=[s.global.ctrlPts(:,1)',s.global.objPts(:,1)'];
OP=[s.global.ctrlPts(:,2:4)',s.global.objPts(:,2:4)'];
isCtrl=ismember(OPid,s.global.ctrlPts(:,1));

% Project 3D points.
pp=cell(size(s.raw.projections));
for i=1:length(ims)
    pp{i}=euclidean(s.K*s.global.P(:,:,i)*homogeneous(OP(:,s.vis(:,i) ...
                                                      & ~isCtrl')));
end

% Measured points.
mm=cell(size(ims));
mmSize=cell(size(ims));
for i=1:length(ims)
    j=s.markPts.obj(:,1)==i;
    mm{i}=s.markPts.obj(j,3:4)';
    mmSize{i}=s.markPts.sz(j);
end

objRes=cat(2,pp{:})-cat(2,mm{:});
keyPtSize=cat(1,mmSize{:})';
objResNorm=sqrt(sum(objRes.^2,1));
RMSreprojScaled=sqrt(mean((objResNorm./keyPtSize).^2))
RMSreprojPx=sqrt(mean(objResNorm.^2))

maxReprojScaled=max(objResNorm./keyPtSize)
maxReprojPx=max(objResNorm)

meanKeyPtSize=mean(keyPtSize)