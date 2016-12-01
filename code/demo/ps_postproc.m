function [rr,E,s0,prob,psz]=ps_postproc(fileName,nRays,minAngle,pauseMode)
%PS_POSTPROC Post-process a PhotoScan project.
%
%   PS_POSTPROC(FILENAME), loads the PhotoScan .psz file in
%   FILENAME and runs the bundle adjustment using the PhotoScan
%   results as initial values.
%
%   PS_POSTPROC(FILENAME,NRAYS), removes all measurements of object
%   points with NRAYS rays or less before processing.
%
%   PS_POSTPROC(FILENAME,NRAYS,ANGLE), removes all measurements of
%   object points with an intersection angle below ANGLE degrees
%   before processing. The intersection angle is computed from
%   Photoscan EO/OP values.
%
%   PS_POSTPROC(FILENAME,NRAYS,ANGLE,PMODE) runs the demo in pause
%   mode PMODE. See PLOTNETWORK for pause modes.
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

if nargin<4, pauseMode='off'; end
if nargin<3, minAngle=0; end
if nargin<2, nRays=0; end

% Extract dir of input file.
[inputDir,inputName,inputExt]=fileparts(fileName);

fprintf('Loading PhotoScan project file %s...',fileName);
psz=loadpsz(fileName);
fprintf('done.\n');

% Conver to Photomodeler structure.
prob=ps2pmstruct(psz);

% Convert to DBAT structure.
s0=prob2dbatstruct(prob);

if nRays>0 || minAngle>0
    if nRays>0
        tooFewRayPts=sum(s0.vis,2)<=nRays & ~s0.isCtrl;
    else
        tooFewRayPts=false;
    end
    
    if minAngle>0
        rayAng=angles(s0,'Computing ray angles')*180/pi;

        tooNarrowAnglePts=rayAng<minAngle & ~s0.isCtrl;
    else
        tooNarrowAnglePts=false;
    end

    % Remove bad points.
    badPts=tooFewRayPts | tooNarrowAnglePts;
    ids2remove=s0.OPid(badPts);
    prob.objPts(ismember(prob.objPts(:,1),ids2remove),:)=[];
    prob.markPts(ismember(prob.markPts(:,2),ids2remove),:)=[];

    % Re-convert to DBAT structure.
    s0=prob2dbatstruct(prob);
end

if psz.camera.isAdjusted
    % Auto-calibration
    s0.estIO(3)=psz.camera.adjustedParams.f; % f
    s0.estIO(1:2)=psz.camera.adjustedParams.cxcy; % cx,cy
    s0.estIO(4:6)=psz.camera.adjustedParams.k(1:3); % K1,K2,K3
    s0.estIO(7:8)=psz.camera.adjustedParams.p(1:2); % P1,P2
    if any(s0.estIO(4:8))
        warning(['Ki/Pi values estimated by Photoscan used as initial ' ...
                 'values for Photomodeler lens distortion model.']);
    end
end

%TODO: Offset estimation.
meanOffset=zeros(3,1);

% Warn for non-uniform mark std.
uniqueSigmas=unique(s0.markStd(:));

% Warn for multiple sigmas for PhotoScan projects only if
% some sigma is zero.
if length(uniqueSigmas)~=1 && any(s0.markStd(:)==0)
    uniqueSigmas
    warning('Multiple mark point sigmas')
    s0.markStd(s0.markStd==0)=1;
end

s=s0;
h=plotnetwork(s,'title','Initial network from PhotoScan',...
              'axes',tagfigure(mfilename),'camsize',1);

% Set up to run the bundle.
damping='gna';

fprintf('Running the bundle with damping %s...\n',damping);

% Run the bundle.
[result,ok,iters,sigma0,E]=bundle(s,damping,20,'trace');
    
if ok
    fprintf('Bundle ok after %d iterations with sigma0=%.2f (%.2f pixels)\n',...
            iters,sigma0,sigma0*s.prior.sigmas(1));
else
    fprintf(['Bundle failed after %d iterations. Last sigma0 estimate=%.2f ' ...
             '(%.2f pixels)\n'],iters,sigma0,sigma0*s.prior.sigmas(1));
end


% Write report file and store computed OP covariances.
reportFile=fullfile(inputDir,[inputName,'-dbatreport.txt']);

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

fprintf('\nPost-processing of PhotoScan project %s.\n',fileName);
fprintf(['%d images, %d CP, %d OP, sigmaCP=%s, m=%d, ' ...
         'n=%d, r=%d.\n'],nImages,nCP,nOP,sigmaCPstr,m,n,r);
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

fprintf(['  sigma0 (DBAT)    : %g\n'], E.s0);

h=plotparams(result,E,'noop','noeo');

h=plotcoverage(result,true);

h=plotimagestats(result,E);

%h=plotopstats(result,E,COP);

fig=tagfigure('networkplayback');

fprintf('Displaying bundle iteration playback for method %s in figure %d.\n',...
        E.damping.name,double(fig));
h=plotnetwork(result,E,...
              'title',['Damping: ',E.damping.name,'. Iteration %d of %d'], ...
              'axes',fig,'pause',pauseMode,'camsize',1); 

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
