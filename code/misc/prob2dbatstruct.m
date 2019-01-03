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
%   - IO     - struct with internal orientation (IO) data.
%   - EO     - struct with external orientation (EO) data.
%   - OP     - struct with object points (OP) data.
%   - IP     - struct with image points (IP) data.
%   - prior  - struct with prior observations of IO, EO, and OP parameters.
%   - bundle - struct with data related to the bundle process.
%   - post   - struct with post-bundle data.
%   - proj   - struct with global project information.
%
%   The IO field contains the following subfields:
%   - model 
%     - distModel  - 1-by-nImages with the used lens distortion model.
%     - nK         - number of radial distortion coefficients.
%     - nP         - number of tangential distortion coefficients.
%     - camUnit    - string with the unit of the physical camera parameters.
%   - val          - NC-by-nImages array with estimates of the internal
%                    orientation for each camera, where NC=5+nK+nP.
%   - type         - NC-by-nImages cell array with strings indicating the
%                    type of parameter, e.g. cc, px, K1, etc.
%   - struct       - struct indicating the block structure of the IO
%                    parameters. See below.
%   - sensor       - sensor/image size information
%     - imSize     - 2-by-nImages array with image [w;h] size in pixels.
%     - ssSize     - 2-by-nImages array with sensor [w;h] in physical units.
%     - pxSize     - 2-by-nImages array with pixel [w;h] in physical units.
%     - samePxSize - logical, true if all pixel sizes are equal.
%
%   The EO field contains the following subfields:
%   - cam    - 1-by-nImages array with physical camera number in IO.
%   - model  - 1-by-nImages array with rotation matrix model. 0=Euler x-y-z.
%   - val    - 6-by-nImages array with EO parameters [X;Y;Z;omega;phi;kappa].
%   - type   - 6-by-nImages cell array of strings indicating parameter type;
%              EX, EY, EZ, om, ph, or, ka.
%   - name   - 1-by-nImages cell array with image names. See also proj.imDir.
%   - id     - 1-by-nImages array with positive integer image ids.
%   - label  - 1-by-nImages cell array with image labels.
%   - struct - struct indicating the block structure of the EO parameters.
%              See below.
%
%   The OP structure contains the following fields:
%   - val   - 3-by-nOPs array with object point coordinates.
%   - type  - 3-by-nOPs cell array with strings indicating
%             parameter type, e.g. OX, OY, OZ for object points,
%             CX, CY, CZ for control points, HX, HY, HZ for check points.
%   - id    - 1-by-nOPs array with positive integer point ids.
%   - rawId - 1-by-nOPs array with native integer ids (depends on source).
%   - label - 1-by-nOPs cell array of strings with point labels.
%      
%   The IP structure contains the following fields:
%   - val    - 2-by-nIPs array with image measurents.
%   - std    - 2-by-nIPs array with assumed standard deviation of the
%              image measurements.
%   - cov    - empty or 2-by-2-by-nIPs with individual covariance
%              for each IP measurement.
%   - type   - 1-by-nIPs cell array of strings with marker types. (Unused).
%   - id     - 1-by-nIPs array with positive integer id for each IP.
%   - cams   - 1-by-nIPs array with camera number (index into IO) for each IP.
%   - vis    - sparse nOPs-by-nImages logical array indicating
%              which IPs are measurements of which OPs.
%              vis(i,j)=true if OP number i was measured in image j
%   - ix     - sparse nOPs-by-nImages array indicating which IP
%              corresponds to a certain measurement. If
%              vis(i,j)=true, then ix(i,j) is the IP number for the
%              measurement of OP i in image j.
%   - sigmas - vector with used sigma values for the IPs.
%
%   The IO.struct and EO.struct fields contain the IO/EO block
%   structure and has the following fields:
%   - block    - NS-by-nImages array, NS=NC (IO) or NS=6 (EO), with
%                numbering indicating what IO/EO parameter values are
%                distinct. Block-variant projects have only one unique
%                value. Image-variant projects have all values distinct.
%   - uniq     - 1-by-nImages logical array indicating if the camera/station
%                is unique.
%   - no       - 1-by-nImages array indicating camera/station number among
%                of the unique numbers.
%   - isSimple - 1-by-nImages logical array indicating if the camera/station
%                is simple, i.e. contain parameters from one block only.
%   - leading  - NS-by-nImages logical array indicating which parameters 
%                are leading a block, i.e. are the first parameter
%                in each row that are to be estimated.
%   The block field should be populated by the user. The uniq, no, and
%   isSimple fields are populated by PARSEBLOCKVARIANT from the block field.
%   The leading field is populated by BUILDSERIALINDICES based on the block
%   info and the corresponding bundle.est.XX and prior.XX.use fields.
%
%   The prior field contain the fields IO, EO, and OP fields with prior
%   observations information, each with subfields:
%   - use - NS-by-nObs logical array indicating whether the corresponding
%           parameter has a prior observation. NS=NC (IO), 6 (EO), or 3
%           (OP). nObs=nImages (IO or EO) or nOPs (OP).
%   - val - NS-by-nObs array with prior values, or NaN if no prior value.
%   - std - NS-by-nObs array with prior standard deviations, or NaN if no
%           prior value. Exact observations are indicated by std=0.
%   - cov - Empty or NS-by-NS-by-nObs array with prior covariance for
%           each columns of val. If empty, the std values are used instead.
%
%   The prior.OP furthermore has two fields:
%   - isCtrl  - 1-by-nOP logical array indicating whether the OP is a
%               control point.
%   - isCheck - 1-by-nOP logical array indicating whether the OP is a
%               check point.
%
%   The bundle field contains information related to the bundle
%   adjustment, with fields
%   - est       - struct with information about which parameters should be
%                 estimated by the bundle.
%     - IO      - NC-by-nImages logical array indicating which IO parameters
%                 should be estimated. Used together with IO.struct.block.
%     - EO      - 6-by-nImages logical array indicating which EO parameters
%                 should be estimated. Used together with EO.struct.block.
%     - OP      - 3-by-nOPs logical array indicating which OP parametes
%                 should be estimated.
%   - serial    - struct with indices describing how to serialize the
%                 bundle data, i.e. generate an x vector.
%     - IO.src  - where from in IO should the values be copied?
%     - IO.dest - where in x should the values end up?
%     - IO.obs  - what IO values should be used as observations?
%     - EO.src  - where from in EO should the values be copied?
%     - EO.dest - where in x should the values end up?
%     - EO.obs  - what EO values should be used as observations?
%     - OP.src  - where from in OP should the values be copied?
%     - OP.dest - where in x should the values end up?
%     - OP.obs  - what OP values should be used as observations?
%     - n       - total number of unknowns.
%   - deserial  - struct with indices describing how to deserialize the x
%                 vector, i.e. update the IO, EO, OP values from x.
%     - IO.src  - where from in x should the IO values be copied?
%     - IO.dest - where in IO should the elements end up?
%     - EO.src  - where from in x should the EO values be copied?
%     - EO.dest - where in EO should the elements end up?
%     - OP.src  - where from in x should the OP values be copied?
%     - OP.dest - where in OP should the elements end up?
%     - n       - total number of unknowns.
%   - resIx     - struct with indices into the residual vector, or 0 if unused.
%     - IP      - 2-by-nIPs with indices for IP residuals
%     - IO      - NC-by-nImages with residual indices for IO obs.
%     - EO      - 6-by-nImages with residual indices for EO obs.
%     - OP      - 3-by-nImages with residual indices for OP obs.
%
%   The post field contains information computed by the bundle.
%   - sigma0    - estimated standard deviation of unit weight.
%   - res, wres - computed unweighted/weighted residuals for each observation,
%                 or NaN if not used as observation
%     - IP      - 2-by-nIPs with IP residuals.
%     - IO      - NC-by-nImages with IO residuals.
%     - EO      - 6-by-nImages with EO residuals.
%     - OP      - 3-by-nImages with OP residuals.
%   - sigmas    - rescaled prior standard deviations (IP.sigmas*sigma0).
%   - std       - posterior standard deviations
%     - IO      - NC-by-nImages with IO posterior std.
%     - EO      - 6-by-nImages with EO posterior std.
%     - OP      - 3-by-nImages with OP posterior std.
%   - cov       - posterior covariances
%     - IO      - NC-by-NC-by-nImages with posterior IO covariance.
%     - EO      - 6-by-6-by-nImages with posterior EO covariance.
%     - OP      - 3-by-3-by-nOPs with posterior OP covariance.
%   - sensor    - posterior estimates of the sensor size
%     - imSize  - 2-by-nImages array with image [w;h] size in pixels.
%     - ssSize  - 2-by-nImages array with sensor [w;h] in physical units.
%     - pxSize  - 2-by-nImages array with pixel [w;h] in physical units.
%
%   The proj field contains global information about the project.
%   - objUnit  - string with the object space unit.
%   - x0desc   - comment string on the initial values used by bundle.
%   - title    - title string.
%   - imDir    - string with the image directory.
%   - fileName - name of the original project file.
%   - cptFile  - name of the control point file.
%   - EOfile   - name of prior EO observations.
%
%See also: LOADPM, BUILDPARAMTYPES, BUILDSERIALINDICES.


if nargin<2, individualCameras=false; end

if ~isstruct(prob)
    error(['PROB2DBATSTRUCT: parameter PROB is not a structure. Was ' ...
           'loading it OK?']);
end

nImages=length(prob.images);
nOP=length(unique([prob.ctrlPts(:,1);prob.objPts(:,1)]));

% Number of lens distortion coefficients.
nK=3;
nP=2;
ixK=5+(1:nK);
ixP=5+nK+(1:nP);

% Internal orientation.
IO=nan(5+nK+nP,nImages);
IOstd=nan(size(IO));
IOcov=[];

if individualCameras
    % Image-invariant
    inner=cat(1,prob.images.inner)';
    innerStd=cat(1,prob.images.innerStd)';
    imSz=reshape(cat(1,prob.images.imSz),2,[]);
    IOblock=repmat(1:nImages,5+nK+nP,1);
else
    % Block-invariant
    inner=repmat(prob.job.defCam,1,nImages);
    innerStd=repmat(prob.job.defCamStd,1,nImages);
    imSz=repmat(prob.job.imSz(:),1,nImages);
    IOblock=ones(5+nK+nP,nImages);
end

% Principal point. Flip y coordinate.
IO(2:3,:)=diag([1,-1])*inner(2:3,:);
IOstd(2:3,:)=innerStd(2:3,:);
% Camera constant.
IO(1,:)=inner(1,:);
IOstd(1,:)=innerStd(1,:);
% Radial distortion coefficients K1, K2, K3.
IO(ixK,:)=-inner(5+(1:nK),:);
IOstd(ixK,:)=innerStd(5+(1:nK),:);
% Tangential distortion coefficients P1, P2.
IO(ixP,:)=-inner(5+nK+(1:nP),:);
IOstd(ixP,:)=innerStd(5+nK+(1:nP),:);

% Sensor size in camera units.
sensorSize=inner(4:5,:);

% Pixel size
pixelSize=sensorSize./imSz;

% Use y as the pixel size. Store diff as aspect.
aspect=1-pixelSize(1,:)./pixelSize(2,:);
pixelSize=pixelSize([2,2],:);
samePxSize=isscalar(unique(pixelSize));

% No skew.
skew=0;

IO(4,:)=aspect;
IO(5,:)=skew;

% External orientation.    
EO=nan(6,nImages);
EOcov=[];

outer=cat(1,prob.images.outer)';
outerStd=cat(1,prob.images.outerStd)';
% Camera centers.
EO(1:3,:)=outer(1:3,:);
EOstd(1:3,:)=outerStd(1:3,:);
% Euler angles, stored as kappa, phi, omega in PM file.
EO(4:6,:)=outer([6,5,4],:)/180*pi;
EOstd(4:6,:)=outerStd([6,5,4],:)/180*pi;

imNames=cellfun(@(x)strrep(x,'\','/'),{prob.images.imName},...
                'uniformoutput',false);

imLabels=cellfun(@(x)strrep(x,'\','/'),{prob.images.label},...
                 'uniformoutput',false);

camIds=cell2mat({prob.images.id});

% Default to no common cam stations.
EOblock=repmat(1:nImages,6,1);

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
%OPstd=nan(3,nOP);
priorCP=nan(size(OP));
priorCPstd=nan(size(OP));
[OPid,i]=sort(prob.objPts(:,1)','ascend'); %#ok<UDIM>
OPrawId=prob.rawOPids(i)';
OPlabels=prob.OPlabels(i)';

isCtrl=ismember(OPid,prob.ctrlPts(:,1)');
isCheck=ismember(OPid,prob.checkPts(:,1)');

% Copy current object point coordinates...
[~,ia,ib]=intersect(OPid,prob.objPts(:,1));
OP(:,ia)=prob.objPts(ib,2:4)';
%OPstd(:,ia)=prob.objPts(ib,5:7)';

% ...and control point coordinates...
[~,ia,ib]=intersect(OPid,prob.ctrlPts(:,1));
priorCP(:,ia)=prob.ctrlPts(ib,2:4)';
priorCPstd(:,ia)=prob.ctrlPts(ib,5:7)';

% ...and check point coordinates...
[~,ia,ib]=intersect(OPid,prob.checkPts(:,1));
priorCP(:,ia)=prob.checkPts(ib,2:4)';
priorCPstd(:,ia)=prob.checkPts(ib,5:7)';

% Find out how many mark points have corresponding object/control points.
imId=unique(prob.markPts(:,1:2),'rows');
nMarkPts=nnz(ismember(imId(:,2),OPid));

% Measured coordinates and visibility matrix.
markPts=nan(2,nMarkPts);
markStd=nan(2,nMarkPts);
IPid=nan(1,nMarkPts);
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
    valid=ismember(measured(:,2)',OPid);
    IPid(ii+(1:nnz(valid)))=measured(valid,2)';
    markPts(:,ii+(1:nnz(valid)))=measured(valid,3:4)';
    markStd(:,ii+(1:nnz(valid)))=measured(valid,5:6)';
    % Update visibility column and pointers.
    vis(:,i)=ismember(OPid,measured(:,2)')';
    colPos(vis(:,i),i)=ii+(1:nnz(valid)); %#ok<SPRIX> % Consider replacing indexing for speed
    ii=ii+nnz(valid);
end

priorSigmas=unique(markStd(:));
if any(priorSigmas==0)
    priorSigmas %#ok<NOPRT>
    warning(['Zero and non-zero prior sigma detected. Using prior sigma==1 ' ...
             'for all image points.']);
    priorSigmas=1;
    markStd(:)=priorSigmas;
end
    
% Pre-calculate which camera corresponds to each point.
[~,j]=find(vis);
ptCams=j';

% Treat IO as fixed.
estIO=false(size(IO));
useIOobs=false(size(IO));
% Treat EO as free.
estEO=true(size(EO));
useEOobs=false(size(EO));
% Estimate all non-fixed OP.
estOP=~(priorCPstd==0);
% Use all non-fixed CP observations. For fixed CP, we only need the
% current "observation".
useOPobs=repmat(isCtrl,3,1);

% Default camera and object space units.
camUnit='mm';
objUnit='m';

% Default to model 1 - backward Brown.
IOdistModel=ones(1,size(IO,2));

% Create project structure.
proj=struct('objUnit',objUnit,...
            'x0desc','',...
            'title',prob.job.title,...
            'imDir',imDir,...
            'fileName',prob.job.fileName,...
            'cptFile','',...
            'EOfile','');

% Create structure for sensor data.
IOsensor=struct('ssSize',sensorSize,...
                'imSize',imSz,...
                'pxSize',pixelSize,...
                'samePxSize',samePxSize);

% Create structure for inner orientation model.
IOmodel=struct('distModel',IOdistModel,...
               'nK',nK,...
               'nP',nP,...
               'camUnit',camUnit);

% Create structure for prior IO data.
priorIO=struct('val',IO,...
               'std',IOstd,...
               'cov',IOcov,...
               'use',useIOobs);

% Create structure for IO data structure.
IOstruct=struct('block',IOblock,...
                'leading',[],...
                'uniq',[],...
                'no',[],...
                'isSimple',[]);

% Create struct with internal orientation data.
IOdata=struct('val',IO,...
             'model',IOmodel,...
             'sensor',IOsensor,...
             'type',[],...
             'struct',IOstruct);

% Create structure for prior EO data.
priorEO=struct('val',EO,...
               'std',EOstd,...
               'cov',EOcov,...
               'use',useEOobs);

% Create structure for EO data structure.
EOstruct=struct('block',EOblock,...
                'leading',[],...
                'uniq',[],...
                'no',[],...
                'isSimple',[]);

EOcam=1:size(EO,2);

% Create structure with external orientation data.
EOdata=struct('val',EO,...
             'model',zeros(1,size(EO,2)),...
             'type',[],...
             'cam',EOcam,...
             'name',{imNames},...
             'id',camIds,...
             'label',{imLabels},...
             'struct',EOstruct);

% Insert prior camera positions.
if ~isempty(prob.priorCamPos)
    [~,i,j]=intersect(EOdata.id,prob.priorCamPos(:,1)');
    priorEO.val(1:3,i)=prob.priorCamPos(j,2:4)';
    priorEO.std(1:3,i)=prob.priorCamPos(j,5:7)';
    priorEO.use(1:3,i)=true;
end

% Create structure for prior OP data.
priorOP=struct('val',priorCP,...
               'std',priorCPstd,...
               'cov',[],...
               'use',useOPobs,...
               'isCtrl',isCtrl,...
               'isCheck',isCheck);

% Create structure with object point data.
OPdata=struct('val',OP,...
             'type',[],...
             'id',OPid,...
             'rawId',OPrawId,...
             'label',{OPlabels});

% Create structure with image point data.
IPdata=struct('val',markPts,...
             'std',markStd,...
             'cov',[],...
             'type',[],...
             'id',IPid,...
             'cam',ptCams,...
             'vis',vis,...
             'ix',colPos,...
             'sigmas',priorSigmas);

% Create struct controlling prior observations.
prior=struct('IO',priorIO,...
             'EO',priorEO,...
             'OP',priorOP);

% Create struct with bundle estimation data.
bundle=struct(...
    'est',struct(...
        'IO',estIO,...
        'EO',estEO,...
        'OP',estOP),...
    'serial',[],...
    'deserial',[]);

% Create struct with posterior bundle data.
post=struct(...
    'res',struct(...
        'IP',[],...
        'IO',[],...
        'EO',[],...
        'OP',[],...
        'ix',[]),...
    'sigmas',[],...
    'std',struct(...
        'IO',[],...
        'EO',[],...
        'OP',[]),...
    'cov',struct(...
        'IO',[],...
        'EO',[],...
        'OP',[]),...
    'sensor',struct(...
        'ssSize',[],...
        'imSize',[],...
        'pxSize',[]));

s=struct('proj',proj,...
         'IO',IOdata,...
         'EO',EOdata,...
         'OP',OPdata,...
         'IP',IPdata,...
         'prior',prior,...
         'bundle',bundle,...
         'post',post);

s=parseblockvariant(s);
s=buildparamtypes(s);
