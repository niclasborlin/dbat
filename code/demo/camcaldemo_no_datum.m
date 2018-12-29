function [rr,s0,prob]=camcaldemo_no_datum(damping,doPause)
%CAMCALDEMO Camera calibration demo for DBAT.
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
inputFile=fullfile(inputDir,'pmexports','camcal-pmexport.txt');
% Report file name.
reportFile=fullfile(inputDir,'dbatexports','camcal-dbatreport-no-datum.txt');;
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

% Set initial IO values to EXIF sensor, pp at center of censor.
s0=setcamvals(s0,'default',7.3);

% Estimate everything but skew.
s0=setcamest(s0,'all','not','sk');

% Estimate all EO parameters.
s0=seteoest(s0,'all');

% Forget to set any control points.

% Save original structure.
saves0=s0;

% Plot initial camera network.
s=s0;
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
