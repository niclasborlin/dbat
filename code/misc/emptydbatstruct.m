function s=emptydbatstruct(nImages,nOP,nIP)
%EMPTYDBATSTRUCT Create an empty DBAT struct.
%
%   S=EMPTYDBATSTRUCT(M,N,P) creates an empty DBAT structure prepared
%   to hold M images, N object points (including ctrl and check pts) and P
%   image points.
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
%     - CIO     - NC*nImages-by-NC*nImages with full IO covariance matrix.
%     - CEO     - 6*nImages-by-6*nImages with full EO covariance matrix.
%     - COP     - 3*nOPs-by-3*nOPs with full OP covariance matrix.
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
%   - UUID     - Random unique UUID of this computation.
%
%See also: PROB2DBATSTRUCT, RUNDBATSCRIPT, BUILDPARAMTYPES, BUILDSERIALINDICES.

% Default number of lens distortion coefficients.
nK=3;
nP=2;

% Mark the units as unknown.
objUnit='<unknown>';
camUnit='<unknown>';

% Create project structure.
proj=struct('objUnit',objUnit,...
            'x0desc','',...
            'title','',...
            'imDir','',...
            'fileName','',...
            'cptFile','',...
            'EOfile','',...
            'UUID',uuid);

% Create default structure for sensor data.
sensorSize=nan(2,nImages);
imSz=nan(2,nImages);
pixelSize=nan(2,nImages);
samePxSize=false;

IOsensor=struct('ssSize',sensorSize,...
                'imSize',imSz,...
                'pxSize',pixelSize,...
                'samePxSize',samePxSize);

% Create default structure for inner orientation model.
IOdistModel=nan(1,nImages);
camUnit='unknown';

IOmodel=struct('distModel',IOdistModel,...
               'nK',nK,...
               'nP',nP,...
               'camUnit',camUnit);

% Create default structure for prior IO data.
IO=nan(5+nK+nP,nImages);
IOstd=nan(size(IO));
IOcov=[];
useIOobs=false(size(IO));

priorIO=struct('val',IO,...
               'std',IOstd,...
               'cov',IOcov,...
               'use',useIOobs);

% Create structure for IO data structure.
IOblock=nan(size(IO));
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

% Create default structure for prior EO data.
EO=nan(6,nImages);
EOstd=nan(size(EO));
EOcov=[];
useEOobs=false(size(EO));

priorEO=struct('val',EO,...
               'std',EOstd,...
               'cov',EOcov,...
               'use',useEOobs);

% Create default structure for EO data.

% Default to no common cam stations.
EOblock=repmat(1:nImages,6,1);

EOstruct=struct('block',EOblock,...
                'leading',[],...
                'uniq',[],...
                'no',[],...
                'isSimple',[]);

EOcam=1:size(EO,2);
imNames=repmat({''},1,nImages);
imLabels=repmat({''},1,nImages);
camIds=nan(1,nImages);

% Create structure with external orientation data.
EOdata=struct('val',EO,...
             'model',zeros(1,size(EO,2)),...
             'type',[],...
             'cam',EOcam,...
             'name',{imNames},...
             'id',camIds,...
             'label',{imLabels},...
             'struct',EOstruct);

% Create default structure for prior OP data.
OP=nan(3,nOP);
priorCP=nan(size(OP));
priorCPstd=nan(size(OP));
isCtrl=false(1,nOP);
isCheck=false(1,nOP);
useOPobs=false(size(OP));
OPid=nan(1,nOP);
OPrawId=nan(1,nOP);
OPlabels=repmat({''},1,nOP);

priorOP=struct('val',priorCP,...
               'std',priorCPstd,...
               'cov',[],...
               'use',useOPobs,...
               'isCtrl',isCtrl,...
               'isCheck',isCheck);

% Create default structure with OP data.
OPdata=struct('val',OP,...
             'type',[],...
             'id',OPid,...
             'rawId',OPrawId,...
             'label',{OPlabels});

% Create default structure with image point data.
markPts=nan(2,nIP);
markStd=nan(size(markPts));
IPid=nan(1,nIP);
ptCams=nan(1,nIP);
vis=logical(sparse(nOP,nImages));
colPos=sparse(nOP,nImages);
priorSigmas=nan;
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

% Create default struct with bundle estimation data.
estIO=false(size(IO));
estEO=false(size(EO));
estOP=false(size(OP));
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
        'times',[],...
        'IO',[],...
        'EO',[],...
        'OP',[],...
        'CIO',[],...
        'CEO',[],...
        'COP',[]),...
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
