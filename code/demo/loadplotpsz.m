function [ppsz,prob,s0,s0PreFilt]=loadplotpsz(fName,sLocal,minVals,plotVals)
%LOADPLOTPSZ Load and plot Photoscan PSZ file.
%
%   PSZ=LOADPLOTPSZ(F), where F is the name of the Photoscan .psz
%   file, loads and plots the camera network and OP that were computed
%   by PhotoScan together with OP and image statistics. The loaded
%   data is returned in the structure PSZ. The processing is performed
%   in global coordinates. For processing in semi-local coordinates
%   (translation and scaling user by Photoscan, but no rotation), use
%   LOADPLOTPSZ(F,SLOCAL) with SLOCAL==TRUE.
%
%   To filter the OP, use LOADPLOTPSZ(F,SLOCAL,MINVALS), where
%   MINVALS=[minRays,minAngle] indicate the lower ray count and angle (in
%   degrees) limits for the OP.
%
%   LOADPLOTPSZ(F,SLOCAL,MINVALS,SZ) sets the width of the
%   camera icon to SZ in object units. SZ defaults to unity.
%
%   If PSZ has been pre-loaded, LOADPSZ(PSZ,...) will process the
%   PSZ structure instead of loading from a file.
%
%   [PSZ,PROB,S0POSTFILT,S0PREFILT]=... also returns the Photomodeler
%   structure PROB and DBAT structures S0* before and after filtering.
%
%See also: PLOTNETWORK, PLOTIMAGESTATS, LOADPLOTDEMO.

if nargin<2, sLocal=false; end
if nargin<3, minVals=[]; end
if nargin<4, plotVals=[]; end

if length(minVals)<1
    minRays=0;
else
    minRays=minVals(1);
end
if length(minVals)<2
    minAngle=0;
else
    minAngle=minVals(2);
end

if length(plotVals)<1
    alignCam=0;
else
    alignCam=plotVals(1);
end

if length(plotVals)<2
    camSz=1;
else
    camSz=plotVals(2);
end

if isstruct(fName)
    % Assume psz structure.
    psz=fName;
    fName='Pre-loaded';
else
    fprintf('Loading PhotoScan project file %s...',fName);
    psz=loadpsz(fName);
    fprintf('done.\n');
end

% Convert to Photomodeler structure.
fprintf('Converting to Photomodeler structure...');
prob=ps2pmstruct(psz,sLocal);
fprintf('done.\n');

% Convert to DBAT structure.
fprintf('Converting to DBAT structure...');
s0=prob2dbatstruct(prob);
fprintf('done.\n');

h=plotnetwork(s0,'title','Initial network from PhotoScan', ...
              'axes',tagfigure([mfilename,'-initnetwork']),'camsize',camSz);
h=plotimagestats(tagfigure([mfilename,'-initimstat']),s0);
pause(0.01);

s0PreFilt=s0;

if minRays>0 || minAngle>0
    fprintf('Filtering OP...');
    if minRays>0
        tooFewRayPts=sum(s0.vis,2)<minRays & ~s0.isCtrl;
    else
        tooFewRayPts=false;
    end
    
    if minAngle>0
        rayAng=angles(s0,'Computing OP ray angles')*180/pi;

        tooNarrowAnglePts=rayAng<minAngle & ~s0.isCtrl;
    else
        tooNarrowAnglePts=false;
    end

    % Remove bad points.
    badPts=tooFewRayPts | tooNarrowAnglePts;
    ids2remove=s0.OPid(badPts);
    prob.objPts(ismember(prob.objPts(:,1),ids2remove),:)=[];
    prob.markPts(ismember(prob.markPts(:,2),ids2remove),:)=[];
    fprintf('done.\n');

    % Re-convert to DBAT structure.
    fprintf('Converting to DBAT structure...');
    s0=prob2dbatstruct(prob);
    fprintf('done.\n');

    titleStr=sprintf(['Filtered network from PhotoScan (minRays=%d, ' ...
                      'minAngle=%.1f)'],minRays,minAngle);
    h=plotnetwork(s0,'title',titleStr,...
                  'axes',tagfigure([mfilename,'-filnetwork']),'camsize',camSz);
    h=plotimagestats(tagfigure([mfilename,'-filtimstat']),s0);
    pause(0.01);
end

if nargout>0, ppsz=psz; end
