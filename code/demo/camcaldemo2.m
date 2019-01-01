function [rr,s0,prob]=camcaldemo2(damping,doPause)
%CAMCALDEMO Camera calibration demo for DBAT.
%
%   CAMCALDEMO runs a camera calibration bundle on a PhotoModeler
%   export file. The project is a 5-image data set of the single 2D
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
inputFile=fullfile(inputDir,'pmexports','camcal-pmexport5.txt');
% Report file name.
reportFile=fullfile(inputDir,'dbatexports','camcal-dbatreport5.txt');;
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

% Initial EO parameters will be computed by resection. Clear all
% values (set to NaN) to catch any errors.
s0=cleareo(s0);

if ~any(s0.prior.OP.isCtrl)
    % No control points specified in PM file, but we know that CPs
    % have ID>1000.
    s0.prior.OP.isCtrl=s0.OP.id>1000;
end

% Load control points
pts=loadcpt(cptFile);

% Match control points with loaded info.
[i,j]=matchcpt(s0,pts);

% Set control points.
s0=setcpt(s0,pts,i,j);

% Save original structure.
saves0=s0;

% Non-ctrl pts are computed by forward intersection. Clear values to
% catch any errors.
s0=clearop(s0);


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
isAbsPath=~isempty(s0.proj.imDir) && ismember(s0.proj.imDir(1),'\\/') || ...
          length(s0.proj.imDir)>1 && s0.proj.imDir(2)==':';
if ~isAbsPath && exist(fullfile(curDir,s0.proj.imDir),'dir')
    % Expand path relative to current dir for this file.
    s0.proj.imDir=fullfile(curDir,s0.proj.imDir);
end
if exist(s0.proj.imDir,'dir')
    % Handle both original-case and lower-case file names.
    imNames={s0.EO.name{imNo},lower(s0.EO.name{imNo}),upper(s0.EO.name{imNo})};    
    imNames=fullfile(s0.proj.imDir,imNames);
    imExist=cellfun(@(x)exist(x,'file')==2,imNames);
    if any(imExist)
        imName=imNames{find(imExist,1,'first')};
    end
else
    warning('Image directory %s does not exist.',s0.proj.imDir);
end

if exist(imName,'file')
    fprintf('Plotting measurements on image %d.\n',imNo);
    imFig=tagfigure('image');
    h=[h;imshow(imName,'parent',gca(imFig))];
    pts=s0.IP.val(:,s0.IP.ix(s0.IP.vis(:,imNo),imNo));
    ptsId=s0.OP.id(s0.IP.vis(:,imNo));
    isCtrl=s0.prior.OP.isCtrl(s0.IP.vis(:,imNo));
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
