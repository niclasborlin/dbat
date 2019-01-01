function [rr,s0,prob,psz]=prague2016_ps(l,doPause)
%PRAGUE2016_PS Run PhotoScan demo in Prague'16 paper.
%
%   PRAGUE2016_PS(LABEL), where LABEL is 'S5', runs the PhotScan
%   experiment of [1].
%
%   PRAGUE2016_PS(LABEL,PAUSE) runs the demo with pause mode
%   PAUSE. See PLOTNETWORK for pause modes.
%
%   References:
%       [1] Börlin and Grussenmeyer (2016), "External Verification
%           of the Bundle Adjustment in Photogrammetric Software
%           using the Damped Bundle Adjustment Toolbox",
%           International Archives of the Photogrammetry, Remote
%           Sensing and Spatial Information Sciences, XLI-B5,
%           p. 7-14. Paper presented at the 2016 ISPRS Congress in
%           Prague, Czech Republic, 12-17 July 2016.
%
%See also: PLOTNETWORK.

if nargin==0, help(mfilename), return, end

if nargin<2, doPause='off'; end

% Extract name of current directory.
curDir=fileparts(mfilename('fullpath'));

switch lower(l(1))
  case 's'
    % Base dir with input files for these projects.
    inputDir=fullfile(curDir,'data','prague2016','sxb');
  otherwise
    error('Bad experiment label');
end

switch lower(l)
  case 's5'
    cpWeighted=true;
    stub='sxb';
  otherwise
    error('Bad experiment label.')
end

% PhotoModeler text export file and report file.
inputFile=fullfile(inputDir,'psprojects',[stub,'.psz']);

fprintf('Loading PhotoScan project file %s...',inputFile);
psz=loadpsz(inputFile);
fprintf('done.\n');

% Control point file.
cpName=fullfile(inputDir,'ref','ctrlpts-weighted-raw.txt');

fprintf('Loading control point file %s...',cpName);
ctrlPts=loadcpt(cpName);
fprintf('done.\n');

[prob,pmReport,pts3d,pts2d]=ps2pmstruct(psz);

s0=prob2dbatstruct(prob);

%TODO: Offset estimation.
meanOffset=zeros(3,1);

% Warn for non-uniform mark std.
uniqueSigmas=unique(s0.IP.std(:));

% Don't warn for multiple sigmas for PhotoScan projects only if
% some sigma is zero.
if length(uniqueSigmas)~=1 && any(s0.IP.std(:)==0)
    uniqueSigmas
    warning('Multiple mark point sigmas')
    s0.IP.std(s0.IP.std==0)=1;
end

% Fixed camera parameters.
s0=setcamvals(s0,'loaded');
s0=setcamest(s0,'none');

% Match control points with loaded info.
[i,j]=matchcpt(s0,ctrlPts);

% Set control points.
s0=setcpt(s0,ctrlPts,i,j);

% Clear EO and OP parameters.
s0=cleareo(s0);
s0=clearop(s0);

% Compute EO parameters by spatial resection.
cpId=s0.OP.id(s0.prior.OP.isCtrl);
s1=resect(s0,'all',cpId,1,0,cpId);
% Compute OP parameters by forward intersection.
s2=forwintersect(s1,'all',true);

s=s2;
h=plotnetwork(s,'title','Initial network (EO, OP computed from CP, IO, MP)',...
              'axes',tagfigure(mfilename),'camsize',100);


% Set up to run the bundle.
damping='gna';

fprintf('Running the bundle with damping %s...\n',damping);

% Run the bundle.
[result,ok,iters,sigma0,E]=bundle(s,damping,'trace','dofverb');
    
if ok
    fprintf('Bundle ok after %d iterations with sigma0=%.2f (%.2f pixels)\n',...
            iters,sigma0,sigma0*s.IP.sigmas(1));
else
    fprintf(['Bundle failed after %d iterations. Last sigma0 estimate=%.2f ' ...
             '(%.2f pixels)\n'],iters,sigma0,sigma0*s.IP.sigmas(1));
end


% Write report file and store computed OP covariances.
reportFile=fullfile(inputDir,'dbatexports',[stub,'-dbatreport.txt']);

[COP,result]=bundle_result_file(result,E,reportFile);

OPstd=full(reshape(sqrt(diag(COP)),3,[]));
CEO=bundle_cov(result,E,'CEO');
EOstd=reshape(full(sqrt(diag(CEO))),6,[]);
EOposStd=EOstd(1:3,:);
EOangStd=EOstd(4:6,:)*180/pi;

fprintf('\nBundle report file %s generated.\n',reportFile);

% Input statistics. Number of images, CP, OP, a priori CP sigma,
% number of observations, number of parameters, redundacy, ray
% count and angle min+max+avg.
nImages=size(s.EO.val,2);
nCP=nnz(s.prior.OP.isCtrl);
nOP=nnz(~s.prior.OP.isCtrl);
sigmaCP=unique(s.prior.OP.std(:,s.prior.OP.isCtrl)','rows')';
if all(sigmaCP==0)
    sigmaCPstr='fixed';
else
    % Determine unit of sigmaCP.
    ls=min(floor(log10(sigmaCP)));
    switch ls
      case -3
        unit='mm';
        base=1e-3;
      case -2
        unit='cm';
        base=1e-2;
      otherwise
        unit='m';
        base=1;
    end
    % Isotropic or not?
    if isscalar(unique(sigmaCP))
        sigmaCPstr=sprintf('%g %s',sigmaCP(1)/base,unit);
    else
        sigmaCPstr=sprintf('(%g,%g,%g) %s',sigmaCP/base,unit);
    end
end
    
m=E.numObs;
n=E.numParams;
r=E.redundancy;

rayCount=full(sum(s.IP.vis,2));
rayAng=angles(result,'Computing ray angles')*180/pi;

OPrayCount=rayCount(~s.prior.OP.isCtrl);
OPrayAng=rayAng(~s.prior.OP.isCtrl);
if isempty(OPrayCount), OPrayCount=0; end
if isempty(OPrayAng), OPrayAng=0; end

CPrayCount=rayCount(s.prior.OP.isCtrl);
CPrayAng=rayAng(s.prior.OP.isCtrl);
if isempty(CPrayCount), CPrayCount=0; end
if isempty(CPrayAng), CPrayAng=0; end

% Compute EO max abs differences.
EOdiff=result.EO.val(1:6,:)-pmReport.EO(1:6,:);
EOstdDiff=nan-pmReport.EOstd(1:6,:);
maxEOposDiff=max(max(abs(EOdiff(1:3,:))));
maxEOangDiff=max(max(abs(EOdiff(4:6,:))))*180/pi;
maxEOposStdDiff=max(max(abs(EOstdDiff(1:3,:))));
maxEOangStdDiff=max(max(abs(EOstdDiff(4:6,:))));

% Compute OP/CP max abs differences.
[~,i,j]=intersect(pts3d.id,result.OP.id);
if length(i)~=length(pts3d.id)
    warning('OP mismatch, disagreeing OP id:');
    disp(setxor(pts3d.id,result.OP.id));
end

OPisCP=result.prior.OP.isCtrl(j);

pmOP=pts3d.pos(:,i);
pmOPstd=pts3d.std(:,i);
dbatOP=result.OP.val(:,j)-repmat(meanOffset,1,size(pmOP,2));
dbatOPstd=OPstd(:,j);
OPdiff=abs(dbatOP-pmOP);
OPstdDiff=abs(dbatOPstd-pmOPstd);
maxOPdiff=max(max(OPdiff(:,~OPisCP)));
maxCPdiff=max(max(OPdiff(:,OPisCP)));
maxOPstdDiff=max(max(OPstdDiff(:,~OPisCP)));
maxCPstdDiff=max(max(OPstdDiff(:,OPisCP)));

% Compare 2d residuals.

% Find mapping from pts2d.id to OPid.
[~,j]=ismember(pts2d.id,result.OP.id);
% Find columns in markPts that correspond to pts2d id, imNo.
cols=result.IP.ix(sub2ind(size(result.IP.ix),j,pts2d.imNo));
res2d=nan(size(result.post.res.IP));
res2d(:,cols)=pts2d.res;
maxRes2d=max(max(abs(res2d-result.post.res.IP)));

fprintf(['\nExperiment %s:\n%d images, %d CP, %d OP, sigmaCP=%s, m=%d, ' ...
         'n=%d, r=%d.\n'],l,nImages,nCP,nOP,sigmaCPstr,m,n,r);
fprintf(['  OP  ray count=%.0f-%.0f (%.1f avg), ray angle=%.0f-%.0f ' ...
         '(%.1f avg) deg\n'],min(OPrayCount),max(OPrayCount),...
        mean(OPrayCount),min(OPrayAng),max(OPrayAng),mean(OPrayAng));       
fprintf(['  CP  ray count=%.0f-%.0f (%.1f avg), ray angle=%.0f-%.0f ' ...
         '(%.1f avg) deg\n'],min(CPrayCount),max(CPrayCount),...
        mean(CPrayCount),min(CPrayAng),max(CPrayAng),mean(CPrayAng));       
fprintf(['  All ray count=%.0f-%.0f (%.1f avg), ray angle=%.0f-%.0f ' ...
         '(%.1f avg) deg\n'],min(rayCount),max(rayCount),...
        mean(rayCount),min(rayAng),max(rayAng),mean(rayAng));       

fprintf('\nResults (project units/degrees/pixels):\n');
noYes={'no','yes'};

if isfield(pmReport,'procOpts')
    fprintf('  PM orient         : %s\n',noYes{pmReport.procOpts.orient+1});
    fprintf('  PM stages/iters   : %d/%d\n',pmReport.totError.numStages,...
            pmReport.totError.numIters);,
end

if isfield(pmReport,'totError') && isfield(pmReport.totError,'lastErr')
    sigma0=pmReport.totError.lastErr;
else
    sigma0=nan;
end

fprintf(['  sigma0 (PM/DBAT)  : %g/%g=%g\n'], sigma0,E.s0,sigma0/E.s0);

fprintf('  EO max pos diff   : %g.\n',maxEOposDiff);
fprintf('  EO max angle diff : %g.\n',maxEOangDiff);
fprintf('  EO max pos std df : %g.\n',maxEOposStdDiff);
fprintf('  EO max ang std df : %g.\n',maxEOangStdDiff);
fprintf('  OP max diff       : %g.\n',maxOPdiff);
fprintf('  CP max diff       : %g.\n',maxCPdiff);
fprintf('  OP max std diff   : %g.\n',maxOPstdDiff);
fprintf('  CP max std diff   : %g.\n',maxCPstdDiff);
fprintf('  2D res max diff   : %g.\n',maxRes2d);

h=plotparams(result,E);

h=plotcoverage(result,true);

h=plotimagestats(result,E);

h=plotopstats(result,E,COP);

fig=tagfigure('networkplayback');

fprintf('Displaying bundle iteration playback for method %s in figure %d.\n',...
        E.damping.name,double(fig));
h=plotnetwork(result,E,...
              'title',['Damping: ',E.damping.name,'. Iteration %d of %d'], ...
              'axes',fig,'pause',doPause,'camsize',100); 

if nargout>0
    rr=result;
end

imName='';
imNo=1;
% Check if image files exist.
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
