function s=prob2dbatstruct(prob,individualCameras)
%PROB2DBATSTRUCT Convert Photomodeler data to DBAT struct.
%
%   S=PROB2DBATSTRUCT(PROB) converts the structure PROB returned by LOADPM
%   to the structure S used by the Damped Bundle Adjustment Toolbox. The
%   default camera is assumed to be used for all images. A warning
%   message is printed if the actual camera parameters differ more than
%   1 standard deviation from the default.
%
%   S=PROB2DBATSTRUCT(PROB,TRUE) forces each image to have its own camera.
%
%   The struct S has fields:
%       IO       - 16-by-nCams array with estimates of the internal
%                  orientation for each camera.
%       EO       - 7-by-nImages array with the external orientation for
%                  each image.
%       cams     - 1-by-nImages numerical array indicating which IO column
%                  correspond to which image. 
%       OP       - 3-by-nOP array with object and control points.
%       OPid     - 1-by-nOP array with object points ids.
%       isCtrl   - 1-by-nOP logical vector indicating which OP are control
%                  points.
%       markPts  - 2-by-nMarkPts array with measured image coordinates in
%                  pixels, stored in image-major order.
%       ptCams   - 1-by-nMarkPts array indicating which IO column
%                  correspond to which measured point.
%       markStd  - 2-by-nMarkPts array with standard deviations for the
%                  markPts coordinates.
%       vis      - nOP-by-nImage sparse logical array indicating in which
%                  image(s) each OP is visible. vis(I,J)==true if object
%                  point I has a measured coordinate in image J.
%       colPos   - nOP-by-nImage numerical array indicating which column
%                  in markPts the corresponding measurement is
%                  stored. colPos(I,J)==K indicates the the measurement of
%                  object point I in image J is stored in column K of markPts.
%       prior    - struct with prior observations
%                  IO     - 16-by-nCams array with prior observations of
%                           the IO parameters, or NaN if no observation.
%                  IOstd  - 16-by-nCams array with prior standard
%                           deviations for the IO parameters, 0 if
%                           exact, NaN if none.
%                  IOcov  - 16-by-16-by-nCams array with prior
%                           covariance matrices for the IO
%                           parameters, or empty if none.
%                  EO     - 7-by-nImages array with prior observations of
%                           the EO parameters, or NaN if none.
%                  EOstd  - 7-by-nImages array with prior standard
%                           deviations for the EO parameters, 0 if exact,
%                           NaN if none.
%                  EOcov  - 6-by-6-by-nImages array with prior covariance
%                           matrices for the EO parameters, or empty if none.
%                  OP     - 3-by-nOP array with prior observations of
%                           control points, NaN if none.
%                  OPstd  - 3-by-nOP array with prior OP standard
%                           deviations, 0 if exact, NaN if none.
%                  OPcov  - 3-by-3-by-nOP array with prior covariance
%                           matrices for the OP parameters, or empty if none.
%                  sigmas - single or multiple sigmas for different
%                           measurement types.
%       residuals - posterior residuals after the bundle
%                  markPt - 2-by-nMarkPts array with mark point residuals
%                           in pixels. Filled in by the bundle. 
%                  IO     - 16-by-nCams array with IO residuals if prior IO
%                           observations were used in the bundle.
%                  EO     - 7-by-nImages array with EO residuals if prior EO
%                           observations were used in the bundle.
%                  OP     - 3-by-nOP array with OP and CP residuals if prior 
%                           OP/CP observations were used in the bundle.
%       sigmas   - vector with a posteriori standard deviations
%                  (prior.sigmas scaled by estimated sigma0).
%       estIO    - 12-by-nCams logical array indicating which internal
%                  parameters should be estimated by the bundle. Defaults
%                  to all false.
%       estEO    - 7-by-nImages logical array indicating which external
%                  parameters should be estimated by the bundle. Defaults
%                  to true for all real camera parameters (first 6 rows).
%       estOP    - 3-by-nOP logical array indicating which OP parameters
%                  are considered free and should be estimated by the
%                  bundle. Defaults to true for all but fixed control
%                  points.
%       useIOobs - 12-by-nCams logical array indicating which prior IO
%                  observations should be used by the bundle. Defaults 
%                  to all false.
%       useEOobs - 7-by-nCams logical array indicating which prior EO
%                  observations should be used by the bundle. Defaults 
%                  to all false.
%       useOPobs - 3-by-nOP logical array indicating which prior OP
%                  observations should be used by the bundle. Defaults 
%                  to true for non-fixed control points.
%       nK       - scalar indicating how many (potentially zero) K values
%                  are used in the model. nK=3.
%       nP       - scalar indicating how many (potentially zero) P values
%                  are used in the model. nK=2.
%       camUnit  - string with the unit used internally by the camera
%                  mm     - nominal mm,
%                  35mm   - '35 mm equivalent' units, i.e. sensor height=24mm,
%                  pixels - pixels,
%                  unity  - sensor height=1.
%       objUnit  - string with the object space unit.
%       x0desc   - comment string on the initial values used by bundle.
%       title    - title string.
%       imNames  - nEO-cell array with image names.
%       imDir    - string with the image directory.
%
%   Each IO column stores the parameters below. Currently, only the first
%   8 may be estimated by the bundle.
%       px,
%       py      - principal point in camera units (typically mm).
%       c       - camera constant in camera units.
%       K1,
%       K2,
%       K3      - radial distortion parameters of Brown (1971).
%       P1,
%       P2      - tangential distortion parameters of Brown (1971).
%       a1,
%       a2      - skew parameters (not used).
%       sw,
%       sh      - sensor width and height in camera units.
%       iw,
%       ih      - image width in pixels.
%       rx,
%       ry      - image resolution.
%
%   Each EO column stores the parameters below. The first 6 parameters
%   may be estimated by the bundle.
%       X,
%       Y,
%       Z       - external coordinates of the camera center in project units.
%       omega,
%       phi,
%       kappa   - Euler angles for the camera orientation in radians.
%       t       - parameter indicating which Euler convention is
%                 used. Currently only t=0 (omega, phi,kappa) is supported.
%
%See also: LOADPM.


if nargin<2, individualCameras=false; end

% Determine number of each type of object.
if individualCameras
    nCams=length(prob.images);
else
    nCams=1;
end

nImages=length(prob.images);
nOP=length(unique([prob.ctrlPts(:,1);prob.objPts(:,1)]));

% Internal orientation.
IO=nan(16,nCams);
IOstd=nan(16,nCams);
IOcov=[];

if individualCameras
    inner=cat(1,prob.images.inner)';
    innerStd=cat(1,prob.images.innerStd)';
    imSz=reshape(cat(1,prob.images.imSz),2,[]);
else
    inner=prob.job.defCam;
    innerStd=prob.job.defCamStd;
    imSz=prob.job.imSz(:);
end

% Principal point. Flip y coordinate.
IO(1:2,:)=diag([1,-1])*inner(2:3,:);
IOstd(1:2,:)=innerStd(2:3,:);
% Camera constant.
IO(3,:)=inner(1,:);
IOstd(3,:)=innerStd(1,:);
% Radial distortion coefficients K1, K2, K3.
nK=3;
IO(3+(1:nK),:)=-inner(5+(1:nK),:);
IOstd(3+(1:nK),:)=innerStd(5+(1:nK),:);
% Tangential distortion coefficients P1, P2.
nP=2;
IO(3+nK+(1:nP),:)=-inner(5+nK+(1:nP),:);
IOstd(3+nK+(1:nP),:)=innerStd(5+nK+(1:nP),:);
% Skew (not used).
IO(3+nK+nP+(1:2),:)=0;
IOstd(3+nK+nP+(1:2),:)=0;
% Sensor size in camera units.
IO(3+nK+nP+2+(1:2),:)=inner(4:5,:);
IOstd(3+nK+nP+2+(1:2),:)=innerStd(4:5,:);
% Image size in pixels.
IO(3+nK+nP+4+(1:2),:)=imSz;
IOstd(3+nK+nP+4+(1:2),:)=0;
% Sensor resolution.
IO(3+nK+nP+6+(1:2),:)=IO(3+nK+nP+4+(1:2),:)./IO(3+nK+nP+2+(1:2),:);
% Fix to force square pixels.
IO(3+nK+nP+6+(1:2),:)=mean(IO(3+nK+nP+6+(1:2),:),1);

% First-order error propagation.
% C=A/B; std(C) = abs(A/B^2)*std(B).
IOstd(3+nK+nP+6+(1:2),:)=...
    abs(IO(3+nK+nP+4+(1:2),:)./...
        IO(3+nK+nP+2+(1:2),:).^2).*IOstd(3+nK+nP+2+(1:2),:);


% External orientation.    
EO=nan(7,nImages);
EOcov=[];

outer=cat(1,prob.images.outer)';
outerStd=cat(1,prob.images.outerStd)';
% Camera centers.
EO(1:3,:)=outer(1:3,:);
EOstd(1:3,:)=outerStd(1:3,:);
% Euler angles, stored as kappa, phi, omega in PM file.
EO(4:6,:)=outer([6,5,4],:)/180*pi;
EOstd(4:6,:)=outerStd([6,5,4],:)/180*pi;
EO(7,:)=0;
EOstd(7,:)=0;

if individualCameras
    cams=1:nImages;
else
    cams=ones(1,nImages);
end

imNames=cellfun(@(x)strrep(x,'\','/'),{prob.images.imName},...
                'uniformoutput',false);

% Find shortest common dir prefix.
imDirs=unique(cellfun(@fileparts,imNames,'uniformoutput',false));

while length(imDirs)>1
    % More than one, check if shortest is prefix to others.
    [~,i]=min(cellfun(@length,imDirs));
    testDir=fullfile(imDirs{i},filesep);
    isPrefixed=strncmp(testDir,imNames,length(testDir));
    if all(isPrefixed)
        % Yes.
        imDirs=imDirs{i};
    else
        % No, trim again.
        imDirs=unique(cellfun(@fileparts,imDirs,'uniformoutput',false));
    end
end

if isempty(imDirs)
    imDir='';
else
    % Pick remaining dir and append / safely.
    imDir=fullfile(imDirs{:},filesep);
end

% Remove image directory from image names.
imNames=cellfun(@(x)x(length(imDir)+1:end),imNames,'uniformoutput',false);

% Object and control points.
OP=nan(3,nOP);
OPstd=nan(3,nOP);
CP=nan(3,nOP);
CPstd=nan(3,nOP);
CPcov=[];
OPid=unique([prob.ctrlPts(:,1);prob.objPts(:,1)]);
isCtrl=ismember(OPid,prob.ctrlPts(:,1));

% Copy object point coordinates...
[~,ia,ib]=intersect(OPid,prob.objPts(:,1));
OP(:,ia)=prob.objPts(ib,2:4)';
OPstd(:,ia)=prob.objPts(ib,5:7)';

% ...and control point coordinates.
[~,ia,ib]=intersect(OPid,prob.ctrlPts(:,1));
CP(:,ia)=prob.ctrlPts(ib,2:4)';
CPstd(:,ia)=prob.ctrlPts(ib,5:7)';

% Find out how many mark points have corresponding object/control points.
imId=unique(prob.markPts(:,1:2),'rows');
nMarkPts=nnz(ismember(imId(:,2),OPid));

% Measured coordinates and visibility matrix.
markPts=nan(2,nMarkPts);
markStd=nan(2,nMarkPts);
% Visibility matrix.
vis=logical(spalloc(nOP,nImages,nMarkPts));
% Column "pointer".
colPos=spalloc(nOP,nImages,nMarkPts);

% Where to put next set of measured points.
ii=0;
for i=1:nImages
    % Extract measured points from this image.
    j=prob.markPts(:,1)==i-1;
    measured=prob.markPts(j,:);
    % Sort by id.
    [~,k]=sort(measured(:,2));
    measured=measured(k,:);
    % Find out which measured points correspond to object/control points.
    valid=ismember(measured(:,2),OPid);
    markPts(:,ii+(1:nnz(valid)))=measured(valid,3:4)';
    markStd(:,ii+(1:nnz(valid)))=measured(valid,5:6)';
    % Update visibility column and pointers.
    vis(:,i)=ismember(OPid,measured(:,2));
    colPos(vis(:,i),i)=ii+(1:nnz(valid)); % Consider replacing indexing for speed
    ii=ii+nnz(valid);
end

priorSigmas=unique(markStd(:));
if isscalar(priorSigmas)
    priorSigma=priorSigmas;
else
    warning('Multiple prior sigmas detected. Using prior sigma==1.');
    priorSigmas %#ok<NOPRT>
    priorSigma=1;
end
    
% Pre-calculate which camera corresponds to each point.
[i,j]=find(vis); %#ok<ASGLU>
ptCams=cams(j);

prior=struct('IO',IO,'IOstd',IOstd,'IOcov',IOcov,...
             'EO',EO,'EOstd',EOstd,'EOcov',EOcov,...
             'OP',CP,'OPstd',CPstd,'OPcov',CPcov,'sigmas',priorSigma);

residuals=struct('markPt',nan(2,nMarkPts),'IO',nan(size(IO)),...
                 'EO',nan(size(EO)),'OP',nan(size(OP)));

% Treat IO as fixed.
estIO=false(size(IO));
useIOobs=false(size(IO));
% Treat EO as free.
estEO=true(size(EO));
estEO(end,:)=false;
useEOobs=false(size(EO));
% Estimate all non-fixed OP.
estOP=~(prior.OPstd==0);
% Use all non-fixed CP observations. For fixed CP, we only need the
% current "observation".
useOPobs=~isnan(prior.OP) & ~(prior.OPstd==0);

% Default camera and object space units.
camUnit='mm';
objUnit='m';

s=struct('title',prob.job.title,'imDir',imDir,'imNames',{imNames}, ...
         'IO',IO,'IOstd',IOstd, ...
         'EO',EO,'EOstd',EOstd, ...
         'cams',cams,...
         'OP',OP,'OPstd',OPstd,'OPid',OPid,...
         'isCtrl',isCtrl,...
         'markPts',markPts,'markStd',markStd,'ptCams',ptCams,...
         'vis',vis,'colPos',colPos,...
         'prior',prior,...
         'residuals',residuals,...
         'estIO',estIO,'estEO',estEO,'estOP',estOP,...
         'useIOobs',useIOobs,'useEOobs',useEOobs,'useOPobs',useOPobs,...
         'nK',nK,'nP',nP,'camUnit',camUnit,...
         'objUnit',objUnit,'x0desc','');
