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
%   The struct S has the following fields:
%       IO       - 16-by-nImages array with estimates of the internal
%                  orientation for each camera.
%       IOblock  - 16-by-nImages array with numbering indicating
%                  what IO values are distinct. Block-variant
%                  projects have only one unique value.
%                  Image-variant projects have all values distinct.
%       EO       - 7-by-nImages array with the external orientation for
%                  each image.
%       EOblock  - 7-by-nImages array with numbering indicating
%                  what EO values are distinct.
%       OP       - 3-by-nOP array with object and control points.
%       OPid     - 1-by-nOP array with object points ids.
%       OPrawId  - 1-by-nOP array with original object point ids.
%       OPlabels - 1-by-nOP cell array with labels of the original ctrl pts.
%       isCtrl   - 1-by-nOP logical vector indicating which OP are control
%                  points.
%       isCheck  - 1-by-nOP logical vector indicating which OP are control
%                  points. Currently set to 
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
%                  IO     - 16-by-nImages array with prior observations of
%                           the IO parameters, or NaN if no observation.
%                  IOstd  - 16-by-nImages array with prior standard
%                           deviations for the IO parameters, 0 if
%                           exact, NaN if none.
%                  IOcov  - 16-by-16-by-nImages array with prior
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
%                  IO     - 16-by-nImages array with IO residuals if prior IO
%                           observations were used in the bundle.
%                  EO     - 7-by-nImages array with EO residuals if prior EO
%                           observations were used in the bundle.
%                  OP     - 3-by-nOP array with OP and CP residuals if prior 
%                           OP/CP observations were used in the bundle.
%       paramTypes - struct with fields IO, EO, OP that indicate
%                  what type of parameter is stored at the
%                  respective position. See PARAMETER TYPES below.
%       sigmas   - vector with a posteriori standard deviations
%                  (prior.sigmas scaled by estimated sigma0).
%       estIO    - 16-by-nImages logical array indicating which internal
%                  parameters should be estimated by the bundle. Defaults
%                  to all false.
%       estEO    - 7-by-nImages logical array indicating which external
%                  parameters should be estimated by the bundle. Defaults
%                  to true for all real camera parameters (first 6 rows).
%       estOP    - 3-by-nOP logical array indicating which OP parameters
%                  are considered free and should be estimated by the
%                  bundle. Defaults to true for all but fixed control
%                  points.
%       serial   - struct with serialisation indices used when
%                  constructing the vector x of unknowns
%                  IOIO - where from in IO should the values be copied?
%                  IOx  - where in x should the values end up?
%                  EOEO - where from in EO should the values be copied?
%                  EOx  - where in x should the values end up?
%                  OPOP - where from in OP should the values be copied?
%                  OPx  - where in x should the values end up?
%       deserial - struct with deserialisation indices used when
%                  deconstructing the vector of unknowns
%                  IOx  - where from in x should the IO values be copied?
%                  IOIO - where in IO should the elements end up?
%                  EOx  - where from in x should the EO values be copied?
%                  EOEO - where in EO should the elements end up?
%                  OPx  - where from in x should the OP values be copied?
%                  OPOP - where in OP should the elements end up?
%       useIOobs - 16-by-nImages logical array indicating which prior IO
%                  observations should be used by the bundle. Defaults 
%                  to all false.
%       useEOobs - 7-by-nImages logical array indicating which prior EO
%                  observations should be used by the bundle. Defaults 
%                  to all false.
%       useOPobs - 3-by-nOP logical array indicating which prior OP
%                  observations should be used by the bundle. Defaults 
%                  to true for non-fixed control points.
%       nK       - scalar indicating how many (potentially zero) K values
%                  are used in the model. Default: nK=3.
%       nP       - scalar indicating how many (potentially zero) P values
%                  are used in the model. Default: nP=2.
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
%       fileName - name of the original project file.
%       imLabels - nEO-cell array with image labels.
%       camId    - nEO-vector with camera ids.
%
%
%   PARAMETER TYPES:
%
%   Each IO column stores the parameters below. Only the first 10
%   may be estimated by the bundle.
%       px,
%       py      - principal point in camera units (typically mm).
%       cc      - camera constant in camera units.
%       K1,
%       K2,
%       K3      - radial distortion parameters of Brown (1971).
%       P1,
%       P2      - tangential distortion parameters of Brown (1971).
%       fa      - affine parameter. Aspect will be (1+fa):1.
%       fs      - skew parameters.
%       sw,
%       sh      - sensor width and height in camera units.
%       iw,
%       ih      - image width in pixels.
%       rx,
%       ry      - image resolution.
%
%   The names above are stored in the paramTypes.IO field. If
%   multiple IO columns are present, the column number is appended.
%
%   Each EO column stores the parameters below. The first 6 parameters
%   may be estimated by the bundle.
%       EX,
%       EY,
%       EZ       - external coordinates of the camera center in project units.
%       omega,
%       phi,
%       kappa   - Euler angles for the camera orientation in radians.
%       tt      - parameter indicating which Euler convention is
%                 used. Currently only t=0 (omega, phi,kappa) is supported.
%
%   The first two letters of the names above are stored in the
%   paramTypes.EO field. If multiple EO columns are present, a camera
%   identifier is appended to each parameter. The camera parameter is
%   consists of the camera sequence number and camera id.
%
%   Each OP column stores the X, Y, Z coordinates.
%
%   The OP names are stored in the paramTypes.OP field. Object points
%   are prefixed with 'O', i.e. 'OX', 'OY', 'OZ'. Control points are
%   prefixed with 'C'. Check points are prefixed with 'H'.
%   Furthermore, a point identifier is appended. The point
%   identifier consists of the sequence number and if necessary,
%   the OP id, the OP raw id, and the OP label.
%
%See also: LOADPM.


if nargin<2, individualCameras=false; end

if ~isstruct(prob)
    error(['PROB2DBATSTRUCT: parameter PROB is not a structure. Was ' ...
           'loading it OK?']);
end

% Determine number of each type of object.
if individualCameras
    nCams=length(prob.images);
else
    nCams=1;
end

nImages=length(prob.images);
nOP=length(unique([prob.ctrlPts(:,1);prob.objPts(:,1)]));

% Internal orientation.
IO=nan(16,nImages);
IOstd=nan(size(IO));
IOcov=[];

if individualCameras
    % Image-invariant
    inner=cat(1,prob.images.inner)';
    innerStd=cat(1,prob.images.innerStd)';
    imSz=reshape(cat(1,prob.images.imSz),2,[]);
    IOblock=repmat(1:nImages,1,16);
else
    % Block-invariant
    inner=repmat(prob.job.defCam,1,nImages);
    innerStd=repmat(prob.job.defCamStd,1,nImages);
    imSz=repmat(prob.job.imSz(:),1,nImages);
    IOblock=ones(16,nImages);
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
% Sensor size in camera units.
IO(3+nK+nP+2+(1:2),:)=inner(4:5,:);
IOstd(3+nK+nP+2+(1:2),:)=innerStd(4:5,:);
% Image size in pixels.
IO(3+nK+nP+4+(1:2),:)=imSz;
IOstd(3+nK+nP+4+(1:2),:)=0;
% Sensor resolution.
IO(3+nK+nP+6+(1:2),:)=IO(3+nK+nP+4+(1:2),:)./IO(3+nK+nP+2+(1:2),:);
% Use y as the pixel size. Store diff as aspect.
pixelSize=IO(3+nK+nP+6+(1:2),:);
IO(3+nK+nP+6+(1:2),:)=pixelSize([2,2],:);
% Aspect
IO(3+nK+nP+1,:)=1-pixelSize(1,:)./pixelSize(2,:);
IOstd(3+nK+nP+1,:)=0; % TODO: Fix this estimate.
% Skew.
IO(3+nK+nP+2,:)=0;
IOstd(3+nK+nP+2,:)=0;

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

imLabels=cellfun(@(x)strrep(x,'\','/'),{prob.images.label},...
                 'uniformoutput',false);

camIds=cell2mat({prob.images.id});

% Default to no common cam stations.
EOblock=repmat(1:nImages,7,1);

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
CCP=nan(3,nOP);
CCPstd=nan(3,nOP);
CCPcov=[];
[OPid,i]=sort(prob.objPts(:,1),'ascend');
OPrawId=prob.rawOPids(i);
OPlabels=prob.OPlabels(i);

isCtrl=ismember(OPid,prob.ctrlPts(:,1));
isCheck=ismember(OPid,prob.checkPts(:,1));

% Copy object point coordinates...
[~,ia,ib]=intersect(OPid,prob.objPts(:,1));
OP(:,ia)=prob.objPts(ib,2:4)';
OPstd(:,ia)=prob.objPts(ib,5:7)';

% ...and control point coordinates...
[~,ia,ib]=intersect(OPid,prob.ctrlPts(:,1));
CP(:,ia)=prob.ctrlPts(ib,2:4)';
CPstd(:,ia)=prob.ctrlPts(ib,5:7)';

% ...and check point coordinates...
[~,ia,ib]=intersect(OPid,prob.checkPts(:,1));
CCP(:,ia)=prob.checkPts(ib,2:4)';
CCPstd(:,ia)=prob.checkPts(ib,5:7)';

OPtypes={'OX','OY','OZ'}';
if size(OP,2)>1
    OPtypes=repmat(OPtypes,1,size(OP,2));
    if any(isCtrl)
        OPtypes(:,isCtrl)=repmat({'CX','CY','CZ'}',1,nnz(isCtrl));
    end
    if any(isCheck)
        OPtypes(:,isCheck)=repmat({'HX','HY','HZ'}',1,nnz(isCheck));
    end

    for i=1:size(OP,2)
        % Id for this OP.
        OPstr=sprintf('-%d',i);
        if OPid(i)~=i
            OPstr=[OPstr,sprintf('/%d',OPid(i))];
        end
        if OPrawId(i)~=OPid(i)
            OPstr=[OPstr,sprintf('/%d',OPrawId(i))];
        end
        if ~isempty(OPlabels{i})
            OPstr=[OPstr,'-',OPlabels{i}];
        end
        OPtypes(:,i)=cellfun(@(x)[x,OPstr],OPtypes(:,i),'uniformoutput',false);
    end
end

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
[~,j]=find(vis);
ptCams=j;

prior=struct('IO',IO,'IOstd',IOstd,'IOcov',IOcov,...
             'EO',EO,'EOstd',EOstd,'EOcov',EOcov,...
             'OP',CP,'OPstd',CPstd,'OPcov',CPcov,...
             'CCP',CCP,'CCPstd',CCPstd,'CCPcov',CCPcov,...
             'sigmas',priorSigma);

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

% Default to model 1 - backward Brown.
IOdistModel=ones(1,size(IO,2));

s=struct('fileName',prob.job.fileName,'title',prob.job.title,'imDir',imDir,...
         'imNames',{imNames},'imLabels',{imLabels},'camIds',camIds,...
         'IO',IO,'IOstd',IOstd,'IOdistModel',IOdistModel,'IOblock',IOblock,...
         'EO',EO,'EOstd',EOstd,'EOblock',EOblock,...
         'OP',OP,'OPstd',OPstd,'OPid',OPid,...
         'OPrawId',OPrawId,'OPlabels',{OPlabels},...
         'isCtrl',isCtrl,'isCheck',isCheck,...
         'markPts',markPts,'markStd',markStd,'ptCams',ptCams,...
         'vis',vis,'colPos',colPos,...
         'prior',prior,...
         'residuals',residuals,...
         'estIO',estIO,'estEO',estEO,'estOP',estOP,...
         'useIOobs',useIOobs,'useEOobs',useEOobs,'useOPobs',useOPobs,...
         'paramTypes',[],...
         'nK',nK,'nP',nP,'camUnit',camUnit,...
         'objUnit',objUnit,'x0desc','');

[IOtypes,EOtypes,OPtypes]=buildparamtypes(s);
s.paramTypes=struct('IO',{IOtypes},'EO',{EOtypes},'OP',{OPtypes});
