function [rr,s0,prob]=riva2014_camcal(l,doPause)
%RIVA2014_CAMCAL Run PhotoModeler demos in Riva'14 paper.
%
%   RIVA2014_CAMCAL(LABEL), where LABEL is 'C1' to 'C5' runs
%   the respective camera calibration networks of [1].
%
%   RIVA2014_CAMCAL(LABEL,PAUSE) runs the demos with pause mode
%   PAUSE. See PLOTNETWORK for pause modes.
%
%   References:
%       [1] Börlin and Grussenmeyer (2014), "Camera Calibration
%           using the Damped Bundle Adjustment Toolbox",
%           ISPRS Annals of the Photogrammetry, Remote
%           Sensing and Spatial Information Sciences, II-5,
%           p. 89-96. Paper presented at the 2014 ISPRS Technical
%           Commision V Symposium, 23-25 June 2014, Riva del Garda,
%           Italy.
%
%See also: PLOTNETWORK.

if nargin<2, doPause='off'; end

% Extract name of current directory.
curDir=fileparts(mfilename('fullpath'));

switch lower(l(1))
  case 'c'
    % Base dir with input files for these projects.
    inputDir=fullfile(curDir,'data','riva2014');
  otherwise
    error('Bad experiment label');
end


switch lower(l)
  case 'c1'
    stub='net1-C404Z-2d';
    cpName='ctrlpts-2d.txt';
  case 'c2'
    stub='net2-canon40d-2d';
    cpName='ctrlpts-2d.txt';
  case 'c3'
    stub='net3-canon7d-2d';
    cpName='ctrlpts-2d.txt';
  case 'c4'
    stub='net4-canon7d-3d';
    cpName='ctrlpts-3d-1.txt';
  case 'c5'
    stub='net5-canon5d-3d';
    cpName='ctrlpts-3d-2.txt';
  otherwise
    error('Bad experiment label.')
end

% PhotoModeler text export file and report file.
inputFile=fullfile(inputDir,'pmexports',[stub,orientStr,'-pmexport.txt']);
reportFile=fullfile(inputDir,'pmexports',[stub,orientStr,'-pmreport.txt']);

fprintf('Loading data file %s...',inputFile);
prob=loadpm(inputFile);
probRaw=prob;
if any(isnan(cat(2,prob.images.imSz)))
    error('Image sizes unknown!');
end
disp('done.')
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

% Replace a posteriori ctrl positions and std by a priori values.
prob.ctrlPts(ia,2:4)=ctrlPts.pos(:,ib)';
prob.ctrlPts(ia,5:7)=ctrlPts.std(:,ib)';

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

% Warn for non-uniform mark std.
uniqueSigmas=unique(s0.markStd(:));

if length(uniqueSigmas)~=1
    uniqueSigmas
    warning('Multiple mark point sigmas')
    s0.markStd(s0.markStd==0)=1;
end

% Clear EO and OP parameters.
s0.EO(s0.estEO)=nan;
s0.OP(s0.estOP)=nan;

% Insert any prior obs to use.
s0.EO(s0.useEOobs)=s0.prior.EO(s0.useEOobs);
s0.OP(s0.useOPobs)=s0.prior.OP(s0.useOPobs);

% Use specified sigma as first approximation.
s0.markStd(:)=s0.prior.sigmas(1);

% Compute EO parameters by spatial resection.
cpId=s0.OPid(s0.isCtrl);
s1=resect(s0,'all',cpId,1,0,cpId);
% Compute OP parameters by forward intersection.
s2=forwintersect(s1,'all',true);

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
            iters,sigma0,sigma0*s.prior.sigmas(1));
else
    fprintf(['Bundle failed after %d iterations. Last sigma0 estimate=%.2f ' ...
             '(%.2f pixels)\n'],iters,sigma0,sigma0*s.prior.sigmas(1));
end

% Write report file and store computed OP covariances.
reportFile=fullfile(inputDir,'dbatexports',[stub,orientStr,'-dbatreport.txt']);

COP=bundle_result_file(result,E,reportFile);

OPstd=full(reshape(sqrt(diag(COP)),3,[]));
CEO=bundle_cov(result,E,'CEO');
EOstd=reshape(full(sqrt(diag(CEO))),6,[]);
EOposStd=EOstd(1:3,:);
EOangStd=EOstd(4:6,:)*180/pi;

fprintf('\nBundle report file %s generated.\n',reportFile);

% Input statistics. Number of images, CP, OP, a priori CP sigma,
% number of observations, number of parameters, redundacy, ray
% count and angle min+max+avg.
nImages=size(s.EO,2);
nCP=nnz(s.isCtrl);
nOP=nnz(~s.isCtrl);
sigmaCP=unique(s.prior.OPstd(:,s.isCtrl)','rows')';
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

rayCount=full(sum(s.vis,2));
rayAng=angles(result,'Computing ray angles')*180/pi;

OPrayCount=rayCount(~s.isCtrl);
OPrayAng=rayAng(~s.isCtrl);
if isempty(OPrayCount), OPrayCount=0; end
if isempty(OPrayAng), OPrayAng=0; end

CPrayCount=rayCount(s.isCtrl);
CPrayAng=rayAng(s.isCtrl);
if isempty(CPrayCount), CPrayCount=0; end
if isempty(CPrayAng), CPrayAng=0; end

% Compute EO max abs differences.
EOdiff=result.EO(1:6,:)-pmReport.EO(1:6,:);
EOstdDiff=result.EOstd(1:6,:)-pmReport.EOstd(1:6,:);
maxEOposDiff=max(max(abs(EOdiff(1:3,:))));
maxEOangDiff=max(max(abs(EOdiff(4:6,:))))*180/pi;
maxEOposStdDiff=max(max(abs(EOstdDiff(1:3,:))));
maxEOangStdDiff=max(max(abs(EOstdDiff(4:6,:))));

% Compute OP/CP max abs differences.
[~,i,j]=intersect(pts3d.id,result.OPid);
if length(i)~=length(pts3d.id)
    warning('OP mismatch, disagreeing OP id:');
    disp(setxor(pts3d.id,result.OPid));
end

OPisCP=result.isCtrl(j);

pmOP=pts3d.pos(:,i);
pmOPstd=pts3d.std(:,i);
dbatOP=result.OP(:,j)-repmat(meanOffset,1,size(pmOP,2));
dbatOPstd=OPstd(:,j);
OPdiff=abs(dbatOP-pmOP);
OPstdDiff=abs(dbatOPstd-pmOPstd);
maxOPdiff=max(max(OPdiff(:,~OPisCP)));
maxCPdiff=max(max(OPdiff(:,OPisCP)));
maxOPstdDiff=max(max(OPstdDiff(:,~OPisCP)));
maxCPstdDiff=max(max(OPstdDiff(:,OPisCP)));

% Compare 2d residuals.

% Find mapping from pts2d.id to OPid.
[~,j]=ismember(pts2d.id,result.OPid);
% Find columns in markPts that correspond to pts2d id, imNo.
cols=result.colPos(sub2ind(size(result.colPos),j,pts2d.imNo));
res2d=nan(size(result.residuals.markPt));
res2d(:,cols)=pts2d.res;
maxRes2d=max(max(abs(res2d-result.residuals.markPt)));

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
    result=rr;
end
