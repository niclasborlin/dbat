function [rr,s0,prob]=camcaldemo(damping,doPause)
%CAMCALDEMO Camera calibration demo for DBAT.
%
%   CAMCALDEMO runs a camera calibration bundle on a PhotoModeler
%   export file. The project is a 21-image data set of the single 2D
%   Photomodeler calibration sheet. The Camera is a Olympus Camedia
%   C4040Z with a 7.25-by-5.44 mm sensor size. As a starting
%   approximation, the EXIF focal length value of 7.3mm is used.
%
%   CAMCALDEMO uses the Gauss-Newton-Armijo damping scheme of [1]
%   by default. Use CAMCALDEMO(DAMPING), where DAMPING is one of
%   - 'none' or 'gm' for classical Gauss-Markov iterations,
%   - 'gna'          Gauss-Newton with Armijo linesearch,
%   - 'lm'           Levenberg-Marquardt, or
%   - 'lmp'          Levenberg-Marquardt with Powell dogleg.
%
%   Use CAMCALDEMO(DAMPING,'off') to visualize the iteration
%   sequence without waiting for a keypress.
%
%   References:
%       [1] Börlin and Grussenmeyer (2013). "Bundle adjustment with
%       and without damping", Photogrammetric Record,
%       vol. 28(144):396-415.

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

fprintf('Running the bundle with damping %s...\n',damping);

% Run the bundle.
[result,ok,iters,sigma0,E]=bundle(s,damping,'trace');
    
if ok
    fprintf('Bundle ok after %d iterations with sigma0=%.2f (%.2f pixels)\n',...
            iters,sigma0,sigma0*s0.prior.sigmas(1));
else
    fprintf(['Bundle failed after %d iterations. Last sigma0 estimate=%.2f ' ...
             '(%.2f pixels)\n'],iters,sigma0,sigma0*s0.prior.sigmas(1));
end

COP=bundle_result_file(result,E,reportFile);

fprintf('\nBundle result file %s generated.\n',reportFile);

h=plotparams(result,E);

h=plotcoverage(result,true);

h=plotimagestats(result,E);

h=plotopstats(result,E,COP);

fig=tagfigure(sprintf('damping=%s',damping));
fprintf('Displaying bundle iteration playback for method %s in figure %d.\n',...
        E.damping.name,double(fig));
h=plotnetwork(result,E,'title',...
              ['Damping: ',damping,'. Iteration %d of %d'], ...
              'axes',fig,'pause',doPause,'camsize',0.1); 

if nargout>0
    rr=result;
end

imName='';
imNo=1;
% Check if image files exist.
if exist(s.imDir,'dir')
    % Handle both original-case and lower-case file names.
    imNames={s.imNames{imNo},lower(s.imNames{imNo}),upper(s.imNames{imNo})};    
    imNames=fullfile(s.imDir,imNames);
    imExist=cellfun(@(x)exist(x,'file')==2,imNames);
    if any(imExist)
        imName=imNames{find(imExist,1,'first')};
    end
else
    warning('Image directory %s does not exist.',s0.imDir);
end

if exist(imName,'file')
    fprintf('Plotting measurements on image %d.\n',imNo);
    imFig=tagfigure('image');
    h=[h;imshow(imName,'parent',gca(imFig))];
    pts=s.markPts(:,s.colPos(s.vis(:,imNo),imNo));
    ptsId=s.OPid(s.vis(:,imNo));
    isCtrl=s.isCtrl(s.vis(:,imNo));
    % Plot non-control points as red crosses.
    if any(~isCtrl)
        line(pts(1,~isCtrl),pts(2,~isCtrl),'marker','x','color','r',...
             'linestyle','none','parent',gca(imFig));
    end
    % Plot control points as black-yellow triangles.
    if any(isCtrl)
        line(pts(1,isCtrl),pts(2,isCtrl),'marker','^','color','k',...
             'markersize',2,'linestyle','none','parent',gca(imFig));
        line(pts(1,isCtrl),pts(2,isCtrl),'marker','^','color','y',...
             'markersize',6,'linestyle','none','parent',gca(imFig));
    end
    for i=1:length(ptsId)
        text(pts(1,i),pts(2,i),int2str(ptsId(i)),'horizontal','center',...
             'vertical','bottom','color','b','parent',gca(imFig));
    end
end
