function [rr,s0,prob]=prague2016_pm(l,orient,doPause)
%PRAGUE2016_PM Run PhotoModeler demos in Prague'16 paper.
%
%   PRAGUE2016_PM(LABEL,ORIENT), where LABEL is 'C1', 'C2', 'S1',
%   'S2', 'S3',' or 'S4' and ORIENT is logical, runs the respective
%   experiments of [1].
%
%   PRAGUE2016_PM(LABEL,ORIENT,PAUSE) runs the demos with pause
%   mode PAUSE. See PLOTNETWORK for pause modes.
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

if nargin<3, doPause='off'; end

% Extract name of current directory.
curDir=fileparts(mfilename('fullpath'));

switch lower(l(1))
  case 'c'
    % Base dir with input files for these projects.
    inputDir=fullfile(curDir,'data','prague2016','cam');
  case 's'
    % Base dir with input files for these projects.
    inputDir=fullfile(curDir,'data','prague2016','sxb');
  otherwise
    error('Bad experiment label');
end

switch lower(l)
  case 'c1'
    cpWeighted=false;
    stub='fixed';
  case 'c2'
    cpWeighted=true;
    stub='weighted';
  case 's1'
    cpWeighted=false;
    stub='f-op0';
  case 's2'
    cpWeighted=true;
    stub='w-op0';
  case 's3'
    cpWeighted=true;
    stub='w-op1';
  case 's4'
    cpWeighted=true;
    stub='wsmart';
  otherwise
    error('Bad experiment label.')
end

if orient
    orientStr='-with-orient';
else
    orientStr='-no-orient';
end

% PhotoModeler text export file and report file.
inputFile=fullfile(inputDir,'pmexports',[stub,orientStr,'-pmexport.txt']);
reportFile=fullfile(inputDir,'pmexports',[stub,orientStr,'-pmreport.txt']);
% PhotoModeler dump files for 3D and 2D points.
input3dFile=fullfile(inputDir,'pmexports',[stub,orientStr,'-3dpts.txt']);
input2dFile=fullfile(inputDir,'pmexports',[stub,orientStr,'-2dpts.txt']);
% PhotoModeler dump files for 3D and 2D smartpoints.
input3dSmartFile=fullfile(inputDir,'pmexports',[stub,orientStr,'-3dsmartpts.txt']);
input2dSmartFile=fullfile(inputDir,'pmexports',[stub,orientStr,'-2dsmartpts.txt']);

% Control point file.
if cpWeighted
    cpName=fullfile(inputDir,'ref','ctrlpts-weighted.txt');
else
    cpName=fullfile(inputDir,'ref','ctrlpts-fixed.txt');
end

fprintf('Loading data file %s...',inputFile);
prob=loadpm(inputFile);
probRaw=prob;
if any(isnan(cat(2,prob.images.imSz)))
    error('Image sizes unknown!');
end
disp('done.')
fprintf('Loading 3D point table %s...',input3dFile);
pts3dNormal=loadpm3dtbl(input3dFile);
fprintf('done.\n');

fprintf('Loading 2D point table %s...',input2dFile);
pts2dNormal=loadpm2dtbl(input2dFile);
fprintf('done.\n');

pts3d=pts3dNormal;
pts2d=pts2dNormal;

if exist(input3dSmartFile,'file')
    fprintf('Loading 3D smart point table %s...',input3dSmartFile);
    pts3dSmart=loadpm3dtbl(input3dSmartFile,true);
    fprintf('done.\n');
    % Merge normal and smart tables.
    pts3d.id=cat(2,pts3dNormal.id,pts3dSmart.id);
    pts3d.name=cat(2,pts3dNormal.name,pts3dSmart.name);
    pts3d.pos=cat(2,pts3dNormal.pos,pts3dSmart.pos);
    pts3d.std=cat(2,pts3dNormal.std,pts3dSmart.std);
    n1=size(pts3dNormal.vis,2);
    n2=size(pts3dSmart.vis,2);
    if n1<n2
        pts3dNormal.vis(1,n2)=0;
    elseif n1>n2
        pts3dSmart.vis(1,n1)=0;
    end
    pts3d.vis=cat(1,pts3dNormal.vis,pts3dSmart.vis);
end

if exist(input2dSmartFile,'file')
    fprintf('Loading 2D smart point table %s...',input2dSmartFile);
    pts2dSmart=loadpm2dtbl(input2dSmartFile);
    fprintf('done.\n');
    % Merge normal and smart tables.
    pts2d.id=cat(2,pts2dNormal.id,pts2dSmart.id);
    pts2d.imNo=cat(2,pts2dNormal.imNo,pts2dSmart.imNo);
    pts2d.pos=cat(2,pts2dNormal.pos,pts2dSmart.pos);
    pts2d.res=cat(2,pts2dNormal.res,pts2dSmart.res);
end

fprintf('Loading PM report file %s...',reportFile);
pmReport=loadpmreport(reportFile);
fprintf('done.\n');

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

% Verify all CPs used by PM are given in CP file.
if ~all(ismember(prob.ctrlPts(:,1),ctrlPts.id))
    pmCtrlPtsId=prob.ctrlPts(:,1)'
    cpFileId=ctrlPts.id
    error('Control point id mismatch.');
end

% Estimate the offset between the world coordinate system and the PM
% bundle coordinate system. The offset range for fixed control points
% should be as small as the difference between the number of digits
% used, typically 1e-3 object units. For weighted control points we
% can also expect a deviation between the a posteriori CP positions
% from the PM export file and the a priori CP positions in the CP
% file.

% Compute the actual offset between CP coordinates from PM and the
% CP file.
[~,ia,ib]=intersect(prob.ctrlPts(:,1),ctrlPts.id);
offset=prob.ctrlPts(ia,2:4)'-ctrlPts.pos(:,ib);
meanOffset=mean(offset,2);
offsetRange=max(offset,[],2)-min(offset,[],2);

% Compute average a priori and a posteriori CP stdev.
avgPreCPStd=mean(ctrlPts.std(:,ib),2);
avgPostCPStd=mean(max(0,prob.ctrlPts(ia,5:7))',2);
avgCPStd=(avgPostCPStd+avgPreCPStd)/2;

% Warn if offset range is above 1e-3 + 2*average CP std.
if max(offsetRange)>1e-3 + 2*(avgPostCPStd+avgPreCPStd)/2
    warning('Large offset range:')
    offsetRange
end

% Adjust a priori control point positions by the offset.
ctrlPts.pos=ctrlPts.pos+repmat(meanOffset,1,size(ctrlPts.pos,2));

% Match control points with loaded info. Overwrite loaded CP names.
[i,j]=matchcpt(s0,ctrlPts,'id');

% Set control points.
s0=setcpt(s0,ctrlPts,i,j);

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

% Write report file and store computed OP covariances.
reportFile=fullfile(inputDir,'dbatexports',[stub,orientStr,'-dbatreport.txt']);

[COP,result]=bundle_result_file(result,E,reportFile);

OPstd=result.post.std.OP;
EOposStd=result.post.std.EO(1:3,:);
EOangStd=result.post.std.EO(4:6,:);

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
EOstdDiff=result.post.std.EO(1:6,:)-pmReport.EOstd(1:6,:);
maxEOposDiff=max(max(abs(EOdiff(1:3,:))));
maxEOangDiff=max(max(abs(EOdiff(4:6,:))))*180/pi;
maxEOposStdDiff=max(max(abs(EOstdDiff(1:3,:))));
maxEOangStdDiff=max(max(abs(EOstdDiff(4:6,:))))*180/pi;

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

fprintf('  PM orient         : %s\n',noYes{pmReport.procOpts.orient+1});
fprintf('  PM stages/iters   : %d/%d\n',pmReport.totError.numStages,...
        pmReport.totError.numIters);,

fprintf(['  sigma0 (PM/DBAT)  ' ...
         ': %g/%g=%g\n'], pmReport.totError.lastErr,E.s0,pmReport.totError.lastErr/E.s0);

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
