function [rr,s0,prob]=camcaldemo_missing_obs(damping,doPause)
%CAMCALDEMO_MISSING_OBS Camera calibration demo for DBAT with missing obs.
%
%   CAMCALDEMO_1RAY runs a camera calibration bundle on a PhotoModeler
%   export file. The project is the same 21-image data set of the
%   single 2D Photomodeler calibration sheet used by CAMCALDEMO,
%   except that all observations of object points 13 and 60 have been
%   removed to trigger a bundle problem due to rank deficient
%   Jacobian.
%
%   References:
%       [1] Börlin and Grussenmeyer (2013). "Bundle adjustment with
%       and without damping", Photogrammetric Record,
%       vol. 28(144):396-415.
%
%See also: CAMCALDEMO.

if nargin<1, damping='gna'; end

if nargin<2, doPause='on'; end

switch damping
  case {'none','gm','gna','lm','lmp'}
    % Do nothing.
  otherwise
    error('Bad damping');
end

% Extract name of current directory.
curDir=fileparts(mfilename('fullpath'));

% Base dir with input files.
inputDir=fullfile(curDir,'data','dbat');

% PhotoModeler text export file and report file.
inputFile=fullfile(inputDir,'pmexports','camcal-pmexport-missing-obs.txt');
% Report file name.
reportFile=fullfile(inputDir,'dbatexports','camcal-dbatreport-missing-obs.txt');;
% Control point file
cptFile=fullfile(inputDir,'ref','camcal-fixed.txt');

fprintf('Loading data file %s...',inputFile);
prob=loadpm(inputFile);
if ~isstruct(prob)
    error('Loading error.');
end

probRaw=prob;
if any(isnan(cat(2,prob.images.imSz)))
    error('Image sizes unknown!');
end
disp('done.')

% Convert loaded PhotoModeler data to DBAT struct.
s0=prob2dbatstruct(prob);

% Switch to lens distortion model that supports skew/aspect.
s0.IO.model.distModel(:)=3;

% Load control points
pts=loadcpt(cptFile);

% Match control points with loaded info.
[i,j]=matchcpt(s0,pts);

% Set control points.
s0=setcpt(s0,pts,i,j);

saves0=s0;

% Set initial IO values to EXIF sensor, pp at center of censor.
s0.IO.val=zeros(size(s0.IO.val));
% c = EXIF value.s
s0.IO.val(1,:)=7.3;
% px,py = center of sensor (sign flip is due to camera model).
s0.IO.val(2:3,:)=0.5*diag([1,-1])*s0.IO.sensor.ssSize;

% Don't use any prior estimates of the IO parameters.
%s0.prior.IO.use=...
% Estimate c,px,py,c,aspect,K1-K3,P1-P2, but not skew (row 5).
s0.bundle.est.IO=repmat((1:10)'~=5,1,size(s0.bundle.est.IO,2));

% Don't use any prior EO paramters.
%s0.prior.EO.use=false(size(s0.EO.val));

% Estimate all EO parameters.
% s0.bundle.est.EO=...

% EO parameters will be computed by resection. Clear all values to NaN
% to catch any errors.
s0.EO.val(:)=NaN;

% Non-ctrl pts are computed by forward intersection. Clear values to
% catch any errors.
s0.OP.val(:,~s0.prior.OP.isCtrl)=NaN;

% Get initial camera positions by spatial intersection.
cpId=s0.OP.id(s0.prior.OP.isCtrl);
[s1,~,fail]=resect(s0,'all',cpId,1,0,cpId);
if fail
    error('Resection failed.');
end
% Get initial object points positions by forward intersection.
s2=forwintersect(s1,'all',true);

s2.proj.x0desc='Camera calibration from EXIF value';

% Plot initial camera network.
s=s2;
h=plotnetwork(s,'title','Initial network',...
              'axes',tagfigure('camcal'),'camsize',0.1);

fprintf('Running the bundle with damping %s...\n',damping);

% Run the bundle.
[result,ok,iters,sigma0,E]=bundle(s,damping,'trace');
    
if ok
    fprintf('Bundle ok after %d iterations with sigma0=%.2f (%.2f pixels)\n',...
            iters,sigma0,result.post.sigmas(1));
else
    fprintf(['Bundle failed after %d iterations (code=%d). Last sigma0 estimate=%.2f ' ...
             '(%.2f pixels)\n'],iters,E.code,sigma0,sigma0*s0.IP.sigmas(1));
end

% Pre-factorize posterior covariance matrix for speed.
E=bundle_cov(result,E,'prepare');

[COP,result]=bundle_result_file(result,E,reportFile);

fprintf('\nBundle result file %s generated.\n',reportFile);
