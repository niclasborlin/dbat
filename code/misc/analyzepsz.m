function analyzepsz(s)
%ANALYZEPSZ Analyze photoscan project
%
%   ANALYZEPSZ(S), where S is a Photoscan project, presents the
%   same numbers as Generate Report function of Photoscan.

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

fprintf('Survey data:\n')
fprintf('  Number of images         : %d\n',nnz(ims));
fprintf('  Camera stations          : %d\n',nnz(ims));
flyAlt=mean(s.global.CC(3,:))-mean(s.global.objPts(:,4));
nTiePts=size(s.global.objPts,2);
fprintf('  Flying altitude (height) : %d m\n',round(flyAlt));
fprintf('  Tie points               : %d\n',nTiePts);

GSD=s.camera.pixelSz(1)*flyAlt/s.camera.focal;
fprintf('  Ground resolution (GSD)  : %.2f cm/pix\n',GSD*100);

nProjs=size(s.markPts.obj,1);
fprintf('  Projections              : %d\n',nProjs);

pp2=cat(2,pp{:});
mm2=cat(2,mm{:});
objRes=pp2-mm2;
keyPtSize=cat(1,mmSize{:})';
objResNorm=sqrt(sum(objRes.^2,1));
RMSreprojPx=sqrt(mean(objResNorm.^2));


fprintf('  Coverage area            : %.2g sq m\n',nan);
fprintf('  Reprojection error       : %.3f pix\n',RMSreprojPx);

fprintf('\nCamera:\n');
fprintf('  Model         : %s\n',s.camera.name);
fprintf('  Type          : %s\n',s.camera.type);
fprintf('  Resolution    : %d x %d\n',s.camera.imSz);
fprintf('  Focal length  : %.3f mm\n',s.camera.focal);
fprintf('  Pixel size    : %g x %g um\n',s.camera.pixelSz*1000);
if s.camera.isAdjusted
    fprintf('  Precalibrated : No\n');
else
    fprintf('  Precalibrated : Yes\n');
end
fprintf('  F             : %.2f\n',s.camera.focal/s.camera.pixelSz(1));
cxcy=s.camera.pp./s.camera.pixelSz-s.camera.imSz/2;
fprintf('  Cx            : %.2f\n',cxcy(1));
fprintf('  Cy            : %.2f\n',cxcy(2));
k=s.camera.k;
k(end+1:4)=0;
for i=1:4
    fprintf('  K%d            : %.2g\n',i,k(i));
end
b=zeros(1,2);
for i=1:2
    fprintf('  B%d            : %.2g\n',i,b(i));
end
p=s.camera.p;
p(end+1:4)=0;
for i=1:4
    fprintf('  P%d            : %.2g\n',i,p(i));
end

fprintf('\nGround Control Points\n');

fprintf(['  Label  X error (um)  Y error (um)  Z error (um)  Total ' ...
         '(um)  Image (pix)\n']);
cpErr=nan(1,3);
for i=find(s.raw.ctrlPtsEnabled)'
    fprintf('%6s  %13g %13g %13g %11g %6.3f (%d)\n',s.raw.ctrlPtsLabels{i},cpErr,...
            norm(cpErr),norm(cpErr)/s.camera.pixelSz(1),...
            nnz(s.markPts.ctrl(:,2)==i));
end
RMSreprojScaled=sqrt(mean((objResNorm./keyPtSize).^2));
maxReprojScaled=max(objResNorm./keyPtSize);
maxReprojPx=max(objResNorm);
meanKeyPtSize=mean(keyPtSize);

fprintf('\nProcessing parameters\n');
fprintf('  General\n');
fprintf('    Cameras                   : %d\n',nnz(ims));
fprintf('    Aligned cameras           : %d\n',nnz(ims));
fprintf('    Markers                   : %d\n',nnz(s.raw.ctrlPtsEnabled));
fprintf('  Point cloud\n');
fprintf('    Points                    : %d of %d\n',size(s.raw.objPts,1),...
        size(s.raw.objPts,1));
fprintf('    RMS reprojection error    : %.6f (%.6f pix)\n',RMSreprojScaled,...
        RMSreprojPx);
fprintf('    Max reprojection error    : %.6f (%.6f pix)\n',maxReprojScaled,...
        maxReprojPx);
fprintf('    Mean key point size       : %.6f pix\n',meanKeyPtSize);
fprintf('    Effective overlap         : %.1f pix\n',full(mean(sum(s.vis,2))));
