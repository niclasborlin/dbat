function [rr,s0,prob]=romabundledemo(damping,doPause)
%ROMABUNDLEDEMO Bundle demo for DBAT.
%
%   ROMABUNDLEDEMO runs the bundle on the PhotoModeler export file
%   of the ROMA data set. The PhotoModeler EO values are used as
%   initial values, except that the EO position are disturbed by
%   random noise with sigma=0.1 m. The OP initial values are
%   computed by forward intersection. The datum is defined by
%   fixing the EO parameters of the first camera and the X
%   coordinate of another camera. The other camera is chosen to
%   maximize the baseline.
%
%   ROMABUNDLEDEMO uses the Gauss-Newton-Armijo damping scheme of [1]
%   by default. Use CAMCALDEMO(DAMPING), where DAMPING is one of
%   - 'none' or 'gm' for classical Gauss-Markov iterations,
%   - 'gna'          Gauss-Newton with Armijo linesearch,
%   - 'lm'           Levenberg-Marquardt, or
%   - 'lmp'          Levenberg-Marquardt with Powell dogleg.
%
%   Use ROMABUNDLEDEMO(DAMPING,'off') to visualize the iteration
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
inputFile=fullfile(inputDir,'pmexports','roma-pmexport.txt');
% Report file name.
reportFile=fullfile(inputDir,'dbatexports','roma-dbatreport.txt');;

fprintf('Loading data file %s...',inputFile);
prob=loadpm(inputFile);
probRaw=prob;
if any(isnan(cat(2,prob.images.imSz)))
    error('Image sizes unknown!');
end
disp('done.')

% Convert loaded PhotoModeler data to DBAT struct.
s0=prob2dbatstruct(prob);
ss0=s0;

% Don't estimate IO data, i.e. treat it as exact.  This block is not
% really necessary, but may serve as a starting point if IO
% parameters are to be estimated.
s0.IO=s0.prior.IO;
s0.estIO=false(size(s0.IO));
s0.useIOobs=false(size(s0.IO));

% Noise sigma [m].
noiseLevel=0.1;

% Reset random number generator.
rng('default');

% Use supplied EO data as initial values. Again, this block is not
% really necessary but may serve as a starting point for modifications.
s0.EO=s0.prior.EO;
s0.estEO(1:6,:)=true; % 7th element is just the axis convention.
s0.useEOobs=false(size(s0.EO));
s0.EO(1:3,:)=s0.EO(1:3,:)+randn(3,size(s0.EO,2))*noiseLevel;

% Copy CP values and treat them as fixed.
s0.OP(:,s0.isCtrl)=s0.prior.OP(:,s0.isCtrl);
s0.estOP=repmat(~s0.isCtrl(:)',3,1);
s0.useOPobs=repmat(s0.isCtrl(:)',3,1);
% Compute initial OP values by forward intersection.
correctedPt=reshape(pm_multilenscorr1(diag([1,-1])*s0.markPts,s0.IO,3,2,...
                                      s0.ptCams,size(s0.IO,2)),2,[]);
s0.OP(:,~s0.isCtrl)=pm_multiforwintersect(s0.IO,s0.EO,s0.cams,s0.colPos,correctedPt,find(~s0.isCtrl));

% Warn for non-uniform mark std.
uniqueSigmas=unique(s0.markStd(:));

if length(uniqueSigmas)~=1
    uniqueSigmas
    error('Multiple mark point sigmas')
end

if all(uniqueSigmas==0)
    warning('All mark point sigmas==0. Using sigma==1 instead.');
    s0.prior.sigmas=1;
    s0.markStd(:)=1;
end

% Fix the datum by fixing camera 1...
s0.estEO(:,1)=false;
% ...and the largest other absolute camera coordinate.
camDiff=abs(s0.EO(1:3,:)-repmat(s0.EO(1:3,1),1,size(s0.EO,2)));
[i,j]=find(camDiff==max(camDiff(:)));
s0.estEO(i,j)=false;

fprintf('Running the bundle with damping %s...\n',damping);

% Run the bundle.
[result,ok,iters,sigma0,E]=bundle(s0,damping,'trace');
    
if ok
    fprintf('Bundle ok after %d iterations with sigma0=%.2f (%.2f pixels)\n',...
            iters,sigma0,sigma0*s0.prior.sigmas(1));
else
    fprintf(['Bundle failed after %d iterations. Last sigma0 estimate=%.2f ' ...
             '(%.2f pixels)\n'],iters,sigma0,sigma0*s0.prior.sigmas(1));
end

COP=bundle_result_file(result,E,reportFile);

fprintf('\nBundle result file %s generated.\n',reportFile);

% Don't plot iteration history for the 26000+ object points.
h=plotparams(result,E,'noop');

h=plotcoverage(result,true);

h=plotimagestats(result,E);

h=plotopstats(result,E,COP);

fig=tagfigure(sprintf('damping=%s',damping));
fprintf('Displaying bundle iteration playback for method %s in figure %d.\n',...
        E.damping.name,double(fig));
plotnetwork(result,E,'trans','up','align',1,'title',...
            ['Damping: ',damping,'. Iteration %d of %d'], ...
            'axes',fig,'pause',doPause);

if nargout>0
    rr=result;
end

imName='';
imNo=1;
% Check if image files exist.
if exist(s0.imDir,'dir')
    % Handle both original-case and lower-case file names.
    imNames={s0.imNames{imNo},lower(s0.imNames{imNo}),upper(s0.imNames{imNo})};    
    imNames=fullfile(s0.imDir,imNames);
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
    pts=s0.markPts(:,s0.colPos(s0.vis(:,imNo),imNo));
    ptsId=s0.OPid(s0.vis(:,imNo));
    isCtrl=s0.isCtrl(s0.vis(:,imNo));
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
