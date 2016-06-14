function prague2016_cam(l,doPause)
%PRAGUE2016_CAM
%
%   PRAGUE2016_CAM(LABEL), where LABEL is 'C1' or 'C2' runs the
%   respective experiments of [1].
%
%   References:
%       [1] Börlin and Grussenmeyer, "External Verification of the
%           Bundle Adjustment in Photogrammetric Software using the
%           Damped Bundle Adjustment Toolbox", presented at the
%           2016 ISPRS Congress in Prague, Czech Republic, 12-17
%           July 2016.

if nargin<2, doPause='off'; end

switch lower(l)
  case 'c1'
    weighted=false;
  case 'c2'
    weighted=true;
  otherwise
    error('Bad experiment label.')
end

% Extract name of current directory.
curDir=fileparts(mfilename('fullpath'));

if weighted
    stub='weighted-1mm';
else
    stub='fixed';
end

% Base dir with input files for these projects.
inputDir=fullfile(curDir,'data','prague2016','cam');

% PhotoModeler text export file.
inputFile=fullfile(inputDir,'pmexports',[stub,'-pmexport.txt']);
% PhotoModeler dump files for 3D and 2D points.
input3dFile=fullfile(inputDir,'pmexports',[stub,'-3dpts.txt']);
input2dFile=fullfile(inputDir,'pmexports',[stub,'-2dpts.txt']);

% Control point file.
cpName=fullfile(inputDir,['ctrlpts-',stub,'.txt']);

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
pts3d=loadpm3dtbl(input3dFile);
fprintf('done.\n');

fprintf('Loading 2D point table %s...',input2dFile);
pts2d=loadpm2dtbl(input2dFile);
fprintf('done.\n');

% Convert loaded PhotoModeler data to DBAT struct.
s0=prob2dbatstruct(prob);
% Store raw version of the struct.
s0raw=s0;

% Warn for non-uniform mark std.
uniqueSigmas=unique(s0.markStd(:));

if length(uniqueSigmas)~=1
    uniqueSigmas
    error('Multiple mark point sigmas')
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
reportFile=fullfile(inputDir,'dbatexports',[stub,'-dbatreport.txt']);

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
nImages=length(s.EO);
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
rayMin=min(full(sum(s.vis,2)));
rayMax=max(full(sum(s.vis,2)));
rayAvg=mean(full(sum(s.vis,2)));
rayAng=angles(result,'Computing ray angles')*180/pi;
rayAngMin=min(rayAng);
rayAngMax=max(rayAng);
rayAngAvg=mean(rayAng);

fprintf(['Experiment %s: %d images, %d CP, %d OP, sigmaCP=%s, m=%d, ' ...
         'n=%d, r=%d, ray count=%.0f-%.0f (%.1f avg), ray angle=%.0f-%.0f ' ...
         '(%.1f avg) deg\n'],l,nImages,nCP,nOP,sigmaCPstr,m,n,r,...
        rayMin,rayMax,rayAvg,rayAngMin,rayAngMax,rayAngAvg);

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
