function [rr,E,s0,prob,psz]=ps_postproc(fileName,sLocal,minRays,minAngle,pauseMode)
%PS_POSTPROC Post-process a PhotoScan project.
%
%   PS_POSTPROC(FILENAME), loads the PhotoScan .psz file in FILENAME
%   and runs the bundle adjustment using the PhotoScan results as
%   initial values. The processing is performed in global
%   coordinates. For processing in semi-local coordinates (translation
%   and scaling user by Photoscan, but no rotation), use
%   PS_POSTPROC(FILENAME,SLOCAL) with SLOCAL==TRUE.
%
%   If FILENAME is blank, the SXB project from [1] is used.
%
%   PS_POSTPROC(FILENAME,SLOCAL,MINRAYS), removes all measurements of
%   object points with less than MINRAYS rays before processing.
%
%   PS_POSTPROC(FILENAME,SLOCAL,MINRAYS,ANGLE), removes all measurements
%   of object points with an intersection angle below ANGLE degrees
%   before processing. The intersection angle is computed from
%   Photoscan EO/OP values.
%
%   PS_POSTPROC(FILENAME,SLOCAL,MINRAYS,ANGLE,PMODE) runs the demo in
%   pause mode PMODE. See PLOTNETWORK for pause modes.
%
%   Use PS_POSTPROC(PSZ,...) if the Photoscan project has already
%   been loaded by loadpsz.
%
%   To run a self-calibration post-processing, please modify the
%   code block near line 86.
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

if nargin<2, sLocal=false; end
if nargin<3, minRays=0; end
if nargin<4, minAngle=0; end
if nargin<5, pauseMode='off'; end

if isstruct(fileName)
    psz=fileName;
    fileName=psz.fileName;
elseif isempty(fileName)
    curDir=fileparts(mfilename('fullpath'));
    inputDir=fullfile(curDir,'data','prague2016','sxb');
    fileName=fullfile(inputDir,'psprojects','sxb.psz');
    fprintf('Loading PhotoScan project file %s...',fileName);
    psz=loadpsz(fileName);
    fprintf('done.\n');
else
    fprintf('Loading PhotoScan project file %s...',fileName);
    psz=loadpsz(fileName);
    fprintf('done.\n');
end

% Extract dir of input file.
[inputDir,inputName]=fileparts(fileName);

% Load project file.
[psz,prob,s0,s0PreFilt]=loadplotpsz(psz,sLocal,[minRays,minAngle]);
% Set the Forward Brown lens distortion model for all cameras.
s0.IO.model.distModel(:)=-1;

resultFile1=fullfile(inputDir,[inputName,'-psstats-prefilt.txt']);
fprintf('Writing report file %s...',resultFile1);
desc1='Initial, unfiltered statitistics';
writestats(s0PreFilt,resultFile1,desc1);
fprintf('done.\n');

resultFile2=fullfile(inputDir,[inputName,'-psstats-postfilt.txt']);
fprintf('Writing report file %s...',resultFile2);
desc2=sprintf('Filtered statitistics with minRays=%d, minAngle=%g',...
              minRays,minAngle);
s0=writestats(s0,resultFile2,desc2);
fprintf('done.\n');

% Should the camera be self-calibrated?
if psz.camera.isAdjusted
    % Self-calibration.
    s0=setcamest(s0,'none');
    if psz.camera.givenParams.f || psz.camera.optimizedParams.f
        s0=setcamest(s0,'cc');
    end
    if psz.camera.givenParams.cxcy(1) || psz.camera.optimizedParams.cxcy(1)
        s0=setcamest(s0,'px');
    end
    if psz.camera.givenParams.cxcy(2) || psz.camera.optimizedParams.cxcy(2)
        s0=setcamest(s0,'py');
    end
    for i=1:3
        if psz.camera.givenParams.k(i) || psz.camera.optimizedParams.k(i) 
            s0=setcamest(s0,sprintf('K%d',i));
        end
    end
    for i=1:2
        if psz.camera.givenParams.p(i) || psz.camera.optimizedParams.p(i)
            s0=setcamest(s0,sprintf('P%d',i));
        end
    end
end

%TODO: Offset estimation.
%meanOffset=zeros(3,1);

% Warn for non-uniform mark std.
uniqueSigmas=unique(s0.IP.std(:));

% Warn for multiple sigmas for PhotoScan projects only if
% some sigma is zero.
if length(uniqueSigmas)~=1 && any(s0.IP.std(:)==0)
    uniqueSigmas %#ok<NOPRT>
    warning('Multiple mark point sigmas')
    s0.IP.std(s0.IP.std==0)=1;
end

s=s0;
h=plotnetwork(s,'title','Pre-bundle network from PhotoScan',...
              'axes',tagfigure(mfilename),'camsize',1); %#ok<NASGU>

pause(0.1)
% Set up to run the bundle.
damping='gna';

fprintf('Running the bundle with damping %s...\n',damping);

% Run the bundle.
[result,ok,iters,sigma0,E]=bundle(s,damping,20,'trace');
    
if ok
    fprintf('Bundle ok after %d iterations with sigma0=%.2f (%.2f pixels)\n',...
            iters,sigma0,sigma0*s.IP.sigmas(1));
else
    fprintf(['Bundle failed after %d iterations (code=%d). Last sigma0 estimate=%.2f ' ...
             '(%.2f pixels)\n'],iters,E.code,sigma0,sigma0*s0.IP.sigmas(1));
end

% Pre-factorize posterior covariance matrix for speed.
E=bundle_cov(result,E,'prepare');

% Write report file and store computed OP covariances.
reportFile=fullfile(inputDir,[inputName,'-dbatreport.txt']);

[COP,result]=bundle_result_file(result,E,reportFile);

OPstd=full(reshape(sqrt(diag(COP)),3,[])); %#ok<NASGU>
CEO=bundle_cov(result,E,'CEO');
EOstd=reshape(full(sqrt(diag(CEO))),6,[]);
EOposStd=EOstd(1:3,:); %#ok<NASGU>
EOangStd=EOstd(4:6,:)*180/pi; %#ok<NASGU>

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

fprintf('  sigma0 (DBAT)    : %g\n', E.s0);

h=plotparams(result,E,'noop','noeo'); %#ok<NASGU>

h=plotcoverage(result,true); %#ok<NASGU>

h=plotimagestats(result,E); %#ok<NASGU>

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
    h=[h;imshow(imName,'parent',gca(imFig))]; %#ok<NASGU>
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
