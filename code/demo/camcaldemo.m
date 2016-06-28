function [rr,s0,prob]=camcaldemo(doPause)
%CAMCALDEMO Camera calibration demo for DBAT.
%
%   CAMCALDEMO runs a camera calibration bundle on a PhotoModeler
%   export file. The project is a 21-image data set of the single 2D
%   Photomodeler calibration sheet. The Camera is a Olympus Camedia
%   C4040Z with a 7.25-by-5.44 mm sensor size. As a starting
%   approximation, the EXIF focal length value of 7.3mm is used.

if nargin<1, doPause='on'; end

% Extract name of current directory.
curDir=fileparts(mfilename('fullpath'));

% Base dir with input files.
inputDir=fullfile(curDir,'data','dbat');

% PhotoModeler text export file and report file.
inputFile=fullfile(inputDir,'pmexports','camcal-pmexport.txt');
% Report file name.
reportFile=fullfile(inputDir,'dbatexports','camcal-dbatreport.txt');;

fprintf('Loading data file %s...',inputFile);
prob=loadpm(inputFile);
probRaw=prob;
if any(isnan(cat(2,prob.images.imSz)))
    error('Image sizes unknown!');
end
disp('done.')

% Set control points to nominal coordinates.
prob.ctrlPts=[1001,0,1,0,0,0,0
              1002,1,1,0,0,0,0
              1003,0,0,0,0,0,0
              1004,1,0,0,0,0,0];
              
% Convert loaded PhotoModeler data to DBAT struct.
s0=prob2dbatstruct(prob);
ss0=s0;

% Set initial IO values. Copy image size, sensor size, etc. from
% prior values.
s0.IO=s0.prior.IO;
s0.IO(1)=s0.IO(11)/2;  % px = center of sensor
s0.IO(2)=-s0.IO(12)/2; % py = center of sensor (sign is due to camera model)
s0.IO(3)=7.3;          % c = EXIF value.
s0.IO(4:8)=0;          % K1-K3, P1-P2 = 0.

% Don't use any prior estimates of the IO parameters.
s0.useIOobs=false(size(s0.IO));
% Estimate px,py,c,K1-K3,P1-P2.
s0.estIO=false(size(s0.IO));
s0.estIO(1:8,:)=true;

% Don't use any prior EO paramters.
s0.useEOobs=false(size(s0.EO));
% Estimate all EO parameters (except for last line which is axis
% sequence indicator).
s0.estEO(1:6,:)=true;
% Clear all EO parameter to NaN to catch any errors.
s0.EO(1:6,:)=NaN;

% Fix the bundle datum by fixing all control points.
s0.estOP=repmat(~s0.isCtrl(:)',3,1);
% Clear all object points values to catch any errors.
s0.OP(s0.estOP)=NaN;

% Get initial camera positions by spatial intersection.
cpId=s0.OPid(s0.isCtrl);
s1=resect(s0,'all',cpId,1,0,cpId);
% Get initial object points positions by forward intersection.
s2=forwintersect(s1,'all',true);

s2.x0desc='Camera calibration from EXIF value';

% Plot initial camera network.
s=s2;
h=plotnetwork(s,'title','Initial network',...
              'axes',tagfigure('camcal'),'camsize',0.1);

dampings={'none','gna','lm','lmp'};

dampings=dampings(2);

result=cell(size(dampings));
ok=nan(size(dampings));
iters=nan(size(dampings));
sigma0=nan(size(dampings));
E=cell(size(dampings));

for i=1:length(dampings)
    fprintf('Running the bundle with damping %s...\n',dampings{i});

    % Run the bundle.
    [result{i},ok(i),iters(i),sigma0(i),E{i}]=bundle(s,dampings{i},'trace');
    
    if ok(i)
        fprintf('Bundle ok after %d iterations with sigma0=%.2f pixels\n', ...
                iters(i),sigma0(i));
    else
        fprintf(['Bundle failed after %d iterations. Last sigma0 estimate=%.2f ' ...
                 'pixels\n'],iters(i),sigma0(i));
    end
end

COP=bundle_result_file(result{1},E{1},reportFile);

fprintf('\nBundle result file %s generated.\n',reportFile);

h=plotparams(result{1},E{1});

h=plotcoverage(result{1},true);

h=plotimagestats(result{1},E{1});

h=plotopstats(result{1},E{1},COP);

for i=1:length(E)
    h=plotnetwork(result{i},E{i},'title',...
                  ['Damping: ',dampings{i},'. Iteration %d of %d'], ...
                  'axes',tagfigure(sprintf('network%d',i)),...
                  'pause',doPause,'camsize',0.1); 
end
