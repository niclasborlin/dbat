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
%       IO      - 16-by-nCams array with internal orientation for each camera.
%       IOstd   - 16-by-nCams array with standard deviations for the IO
%                 parameters.
%       EO      - 7-by-nImages array with external orientation for each image.
%       EOstd   - 7-by-nImages array with standard deviations for the EO
%                 parameters.
%       cams    - 1-by-nImages numerical array indicating which column in
%                 IO corresponds to each image.
%       OP      - 3-by-nOP array with object and control points.
%       OPstd   - 3-by-nOP array with standard deviations for the OP
%                 coordinates.
%       OPid    - 1-by-nOP array with object points ids.
%       isCtrl  - 1-by-nOP logical vector indicating which OP are control
%                 points.
%       markPts - 2-by-nMarkPos array with measured image coordinates in
%                 pixels, stored in image-major order.
%       markStd - 2-by-nMarkPos array with standard deviations for the
%                 markPts coordinates.
%       vis     - nOP-by-nImage sparse logical array indicating in which
%                 image(s) each OP is visible. vis(I,J)==true if object
%                 point I has a measured coordinate in image J.
%       colPos  - nOP-by-nImage numerical array indicating which column
%                 in markPts the corresponding measurement is
%                 stored. colPos(I,J)==K indicates the the measurement of
%                 object point I in image J is stored in column K of markPts.
%       cIO     - 12-by-nCams logical array indicating which internal
%                 parameters are considered free and should be estimated
%                 by the bundle. Defaults to all false.
%       cEO     - 7-by-nImages logical array indicating which external
%                 parameters are considered free and should be estimated
%                 by the bundle. Defaults to true for all real camera
%                 parameters (first 6 rows).
%       cOP     - 3-by-nOP logical array indicating which OP parameters
%                 are considered free and should be estimated by the
%                 bundle. Defaults to true for object points, false for
%                 control points.
%       nK      - scalar indicating how many (potentially zero) K values
%                 are used in the model. nK=3.
%       nP      - scalar indicating how many (potentially zero) P values
%                 are used in the model. nK=2.
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

% $Id$

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
% First-order error propagation.
% C=A/B; std(C) = abs(A/B^2)*std(B).
IOstd(3+nK+nP+6+(1:2),:)=...
    abs(IO(3+nK+nP+4+(1:2),:)./...
        IO(3+nK+nP+2+(1:2),:).^2).*IOstd(3+nK+nP+2+(1:2),:);


% External orientation.    
EO=nan(7,nImages);

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


% Object and control points.
OP=nan(3,nOP);
OPstd=nan(3,nOP);
OPid=unique([prob.ctrlPts(:,1);prob.objPts(:,1)]);
isCtrl=ismember(OPid,prob.ctrlPts(:,1));

% Copy object point coordinates...
[~,ia,ib]=intersect(OPid,prob.objPts(:,1));
OP(:,ia)=prob.objPts(ib,2:4)';
OPstd(:,ia)=prob.objPts(ib,5:7)';
% ...and control point coordinates (duplicates will overwrites OP coords).
[~,ia,ib]=intersect(OPid,prob.ctrlPts(:,1));
OP(:,ia)=prob.ctrlPts(ib,2:4)';
OPstd(:,ia)=prob.ctrlPts(ib,5:7)';

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
    colPos(vis(:,i),i)=ii+(1:nnz(valid));
    ii=ii+nnz(valid);
end

cIO=false(size(IO));
cEO=true(size(EO));
cEO(end,:)=false;
cOP=repmat(~isCtrl,3,1);

s=struct('IO',IO,'IOstd',IOstd,'EO',EO,'EOstd',EOstd,'cams',cams, ...
         'OP',OP,'OPstd',OPstd,'OPid',OPid,'isCtrl',isCtrl, 'markPts', ...
         markPts,'markStd',markStd,'vis',vis,'colPos',colPos,'cIO',cIO, ...
         'cEO',cEO,'cOP',cOP,'nK',nK,'nP',nP);
