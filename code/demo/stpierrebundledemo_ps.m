function [rr,s0,prob]=stpierrebundledemo_ps(damping,doPause)
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
inputDir=fullfile(curDir,'data','hamburg2017','stpierre');

% Name of control point file.
cpName=fullfile(inputDir,'ctrl_StPierre_weighted.txt');
cpName=fullfile(inputDir,'ctrl_close.txt');

% PhotoModeler text export file and report file.
inputFile=fullfile(inputDir,'psprojects','C5.psz');
inputFile=fullfile(getenv('HOME'),'photoscan','ari-close-range','S0.psz');
% Report file name.
reportFile=fullfile(inputDir,'dbatexports','close-S0-dbatreport.txt');

fprintf('Loading PhotoScan project file %s...',inputFile);
psz=loadpsz(inputFile);
fprintf('done.\n');

[prob,pmReport,pts3d,pts2d]=ps2pmstruct(psz);

probRaw=prob; %#ok<NASGU>
if any(isnan(cat(2,prob.images.imSz)))
    error('Image sizes unknown!');
end
disp('done.')

markPtsStd=unique(prob.markPts(:,5:6));
if isscalar(markPtsStd)
    disp('Same stdev on all mark points:')
    disp(markPtsStd);
else
    disp('Different stdev on some mark points:')
    disp(markPtsStd);
end    
%fprintf('Setting all mark pts std to 1.\n');
%prob.markPts(:,5:6)=1;

fprintf('Loading control point file %s...',cpName);
ctrlPts=loadcpt(cpName);
fprintf('done.\n');

% Remap ctrl ids from PM to PS using the labels.
for i=1:length(ctrlPts.id)
    n=ctrlPts.name{i};
    if isempty(n)
        error('Cannot remap unlabelled control/check point');
    end
    % Where is it among the PM points?
    ix=find(strcmp(n,prob.OPlabels));
    if isscalar(ix)
        % Found exactly one match. Copy its id.
        ctrlPts.id(i)=prob.objPts(ix,1);
    elseif isempty(ix)
        % Found none - remove it.
        warning('Control point ''%s'' not found in PSZ file. Removing it.',n);
        % Mark point for removal.
        ctrlPts.id(i)=nan;
    elseif length(ix)>1
        % Found more than one - cannot handle this.
        error('Duplicate matches for control point ''%s'' in PSZ file.',n);
    end
end

% Remove any missing ctrl pts.
if any(isnan(ctrlPts.id))
    i=~isnan(ctrlPts.id);
    ctrlPts.id=ctrlPts.id(i);
    ctrlPts.name=ctrlPts.name(i);
    ctrlPts.pos=ctrlPts.pos(:,i);
    ctrlPts.std=ctrlPts.std(:,i);
end

% Verify all CPs used by PM are given in CP file.
if ~all(ismember(prob.ctrlPts(:,1),ctrlPts.id))
    pmCtrlPtsId=prob.ctrlPts(:,1)' %#ok<NOPRT,NASGU>
    cpFileId=ctrlPts.id %#ok<NOPRT,NASGU>
    error('Control point id mismatch.');
end

% Determine check point allocation among control points.
[ctrlIds,~,ibCP]=intersect(prob.ctrlPts(:,1),ctrlPts.id);
% Assume that PM object points that are not among the PM control
% points but are listed in the control point file are check points.
checkIds=setdiff(intersect(prob.objPts(:,1),ctrlPts.id),ctrlIds);
[~,~,ibCC]=intersect(checkIds,ctrlPts.id);

% Original coordinates.
refCtrlData=[ctrlPts.id(ibCP);ctrlPts.pos(:,ibCP);ctrlPts.std(:,ibCP)]';
refCheckData=[ctrlPts.id(ibCC);ctrlPts.pos(:,ibCC);ctrlPts.std(:,ibCC)]';

% Expand with 0 stdev for fixed points.
if ~isempty(refCtrlData) && size(refCtrlData,2)<7
    refCtrlData(1,7)=0;
end

if ~isempty(refCheckData) && size(refCheckData,2)<7
    refCheckData(1,7)=0;
end

% Coordinates from PM.
[~,ia]=intersect(prob.objPts(:,1),ctrlIds);
pmCtrlData=prob.objPts(ia,:);
[~,ia]=intersect(prob.objPts(:,1),checkIds);
pmCheckData=prob.objPts(ia,:);

% Compute rigid-body transformation for control points.
TCP=rigidalign(pmCtrlData(:,2:4)',refCtrlData(:,2:4)');
if length(checkIds)>=3
    TCC=rigidalign(pmCheckData(:,2:4)',refCheckData(:,2:4)');
else
    TCC=[];
end

% Display info about CPT
fprintf('\nFound %d control points:\n',length(ctrlIds))
disp(refCtrlData(:,1)')
disp(ctrlPts.name(ismember(ctrlPts.id,ctrlIds)));
angCP=acosd(min(max((trace(TCP(1:3,1:3))-1)/2,0),1));
fprintf('Rotation %.1f degress, translation %.1f object units.\n',...
        angCP,norm(TCP(1:3,4)));
fprintf('Max abs pos diff=%g\n',max(max(abs(refCtrlData(:,2:4)-pmCtrlData(:,2:4)))));
fprintf('Max rel std diff=%.1f%%\n',(exp(max(max(abs(log(pmCtrlData(:,5:7))-log(refCtrlData(:,5:7))))))-1)*100)

% Display info about CC
fprintf('\nFound %d check points:\n',length(checkIds))
disp(refCheckData(:,1)')
disp(ctrlPts.name(ismember(ctrlPts.id,checkIds)));
if ~isempty(TCC)
    angCC=acosd(min(max((trace(TCC(1:3,1:3))-1)/2,0),1));
    fprintf('Rotation %.1f degress, translation %.1f object units.\n',...
            angCC,norm(TCC(1:3,4)));
end
fprintf('Max abs pos diff=%g\n',max(max(abs(refCheckData(:,2:4)-pmCheckData(:,2:4)))));
fprintf('Max rel std diff=%.1f%%\n',(exp(max(max(abs(log(pmCheckData(:,5:7))-log(refCheckData(:,5:7))))))-1)*100)

% Replace PM ctrl pt with prior.
prob.ctrlPts=refCtrlData;

% Insert check points.
prob.checkPts=refCheckData;

% Set OP labels to original CP or CCP names.
[~,ia,ib]=intersect(prob.objPts(:,1),ctrlPts.id);
% Clear first.
prob.OPlabels=cell(size(prob.OPlabels));
prob.OPlabels(ia)=ctrlPts.name(ib);

% Convert loaded PhotoModeler data to DBAT struct.
s0=prob2dbatstruct(prob);
ss0=s0;

% Don't estimate IO data, i.e. treat it as exact.  This block is not
% really necessary, but may serve as a starting point if IO
% parameters are to be estimated.
s0.IO=s0.prior.IO;
s0.estIO=false(size(s0.IO));
s0.useIOobs=false(size(s0.IO));

% Add self-calibration for all non-zero parameters...
%s0.estIO(1:3+s0.nK+s0.nP)=s0.IO(1:3+s0.nK+s0.nP)~=0;

% Estimate f, pp, K1-K3, P1-P2
selfCal=true(8,1);
%selfCal(7:8)=false;
% Indicate that we want to estimate these parameters.
s0.estIO(find(selfCal))=true;
% Zero any unused lens distortion parameters.
s0.IO(find(~selfCal))=0;

% ...or for all lens distortion parameters.
%s0.estIO(1:8)=true;

% Switch to Forward/Computer Vision lens distortion model for all cameras.
s0.IOdistModel(:)=-1;

% Use supplied EO data as initial values. Again, this block is not
% really necessary but may serve as a starting point for modifications.
s0.EO=s0.prior.EO;
s0.estEO(1:6,:)=true; % 7th element is just the axis convention.
s0.useEOobs=false(size(s0.EO));

% Use estimated OP values as initial.

% Warn for non-uniform mark std.
uniqueSigmas=unique(s0.markStd(:));

if length(uniqueSigmas)~=1
    uniqueSigmas
    %error('Multiple mark point sigmas')
end

if all(uniqueSigmas==0)
    warning('All mark point sigmas==0. Using sigma==1 instead.');
    s0.prior.sigmas=1;
    s0.markStd(:)=1;
end

% Datum is given by CPs.
% Do nothing.

fprintf('\nRunning the bundle with damping %s...\n',damping);

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
