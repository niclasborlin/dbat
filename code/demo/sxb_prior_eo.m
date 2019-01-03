function [rr,s0,prob]=sxb_prior_eo(usePriorEO)
%SXB_PRIOR_EO Aerial demo with prior EO position data.
%
%   SXB_PRIOR_EO(TRUE) runs the PRAGUE2016_PM demo with weighted
%   control points and prior observations of the EO positions. The
%   prior EO observations are fake and only used to illustrate how
%   EO observations could be used by the bundle.
%
%   SXB_PRIOR_EO(FALSE) runs the demo without using any prior EO
%   positions.
%
%See also: PRAGUE2016_PM.

if nargin==0, help(mfilename), return, end

% Extract name of current directory.
curDir=fileparts(mfilename('fullpath'));

% Base dir with input files for these projects.
inputDir=fullfile(curDir,'data','prague2016','sxb');

cpWeighted=true;
stub='wsmart';
orientStr='-with-orient';

% PhotoModeler text export file and report file.
inputFile=fullfile(inputDir,'pmexports',[stub,orientStr,'-pmexport.txt']);

% Control point file.
cpName=fullfile(inputDir,'ref','ctrlpts-weighted.txt');

% EO file.
EOname=fullfile(inputDir,'ref','fake-camera-positions.txt');

if usePriorEO
    prior='-prior-eo';
else
    prior='-no-prior-eo';
end

% Report file name.
reportFile=fullfile(inputDir,'dbatexports',['sxb',prior,'-dbatreport.txt']);

fprintf('Loading data file %s...',inputFile);
prob=loadpm(inputFile);
probRaw=prob;
if any(isnan(cat(2,prob.images.imSz)))
    error('Image sizes unknown!');
end
disp('done.')

% Convert loaded PhotoModeler data to DBAT struct.
s0=prob2dbatstruct(prob);
% Store raw version of the struct.
s0raw=s0;

% Fixed camera parameters.
s0=setcamvals(s0,'loaded');
s0=setcamest(s0,'none');

fprintf('Loading control point file %s...',cpName);
ctrlPts=loadcpt(cpName);
fprintf('done.\n');

% Match control points with loaded info. Overwrite loaded CP names.
[i,j]=matchcpt(s0,ctrlPts,'id');

% Set control points.
s0=setcpt(s0,ctrlPts,i,j);

% Load, match and set prior EO observations.
if usePriorEO
    fprintf('Loading EO file %s...',EOname);
    EOtbl=loadeotable(EOname,[false,true]);
    fprintf('done.\n');

    [i,j]=matcheo(s0,EOtbl);

    s0=setprioreo(s0,EOtbl,i,j);
end

% Clear any non-control EO and OP parameters.
s0=cleareo(s0);
s0=clearop(s0);

% Compute EO parameters by spatial resection.
cpId=s0.OP.id(s0.prior.OP.isCtrl);
s1=resect(s0,'all',cpId,1,0,cpId);
% Compute OP parameters by forward intersection.
s2=forwintersect(s1,'all',true);

%s2.IOdistModel(:)=3;

s=s2;
h=plotnetwork(s,'title','Initial network (EO, OP computed from CP, IO, MP)',...
              'axes',tagfigure(mfilename),'camsize',0.1);

% Set up to run the bundle.
damping='gna';

fprintf('Running the bundle with damping %s...\n',damping);

% Run the bundle.
[result,ok,iters,sigma0,E]=bundle(s,damping,'trace','dofverb');
    
if ok
    fprintf('Bundle ok after %d iterations with sigma0=%.2f (%.2f pixels)\n',...
            iters,sigma0,result.post.sigmas(1));
else
    fprintf(['Bundle failed after %d iterations (code=%d). Last sigma0 estimate=%.2f ' ...
             '(%.2f pixels)\n'],iters,E.code,sigma0,sigma0*s0.IP.sigmas(1));
end

% Pre-factorize posterior covariance matrix for speed.
E=bundle_cov(result,E,'prepare');

% Write report file.
[COP,result]=bundle_result_file(result,E,reportFile);

fprintf('\nBundle report file %s generated.\n',reportFile);

h=plotparams(result,E);

h=plotcoverage(result,true);

h=plotimagestats(result,E);

h=plotopstats(result,E,COP);

fig=tagfigure('networkplayback');

fprintf('Displaying bundle iteration playback for method %s in figure %d.\n',...
        E.damping.name,double(fig));
h=plotnetwork(result,E,...
              'title',['Damping: ',E.damping.name,'. Iteration %d of %d'], ...
              'axes',fig,'camsize',0.1); 

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
