function [ppsz,prob,s0,rayAng,camRayAng]=loadplotpsz(fName,sLocal,minVals,plotVals,resultFile)
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
%   MINVALS=[minRays,minAngle] indicate lower ray count and angle (in
%   degrees) limits for the OP.
%
%   LOADPLOTPSZ(F,SLOCAL,MINVALS,SZ) sets the width of the
%   camera icon to SZ in object units. SZ defaults to unity.
%
%   LOADPLOTPSZ(F,SLOCAL,MINVALS,SZ,RF) writes project statistics
%   into a report file RF instead of to the terminal.
%
%   If PSZ has been pre-loaded, LOADPSZ(PSZ,...) will process the
%   PSZ structure instead of loading from a file.
%
%   [PSZ,PROB,S0]=... also returns the Photomodeler and DBAT
%   structures, optionally after filtering.
%
%   [PSZ,PROB,S0,RAYANG,CAMRAYANG]=... also returns the computed ray
%   angles of the final data set.
%
%See also: PLOTNETWORK, PLOTIMAGESTATS, LOADPLOTDEMO.

if nargin<2, sLocal=false; end
if nargin<3, minVals=[]; end
if nargin<4, plotVals=[]; end
if nargin<5, resultFile=''; end

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

if minRays>0 || minAngle>0
    fprintf('Filtering OP...');
    if minRays>0
        tooFewRayPts=sum(s0.vis,2)<=minRays & ~s0.isCtrl;
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

    h=plotnetwork(s0,'title','Filtered network from PhotoScan', ...
                  'axes',tagfigure([mfilename,'-filnetwork']),'camsize',camSz);
    h=plotimagestats(tagfigure([mfilename,'-filtimstat']),s0);
    pause(0.01);
end

if ~isempty(resultFile)
    fprintf('Opening report file %s...',resultFile);
    [fid,message]=fopen(resultFile,'wt+');
    if fid>0
        fprintf('done.\n');
    else
        error('Failed to open report file: %s.\n',message);
    end
else
    fid=1;
end

fprintf(fid,'Project file: %s\n',fName);

nCp=nnz(s0.isCtrl);
fprintf(fid,'\nTotal # OP         : %d\n',size(s0.vis,1)-nCp);
fprintf(fid,'Total # CP         : %d\n',nCp);
fprintf(fid,'Total # cams       : %d\n',size(s0.vis,2));
fprintf(fid,'Total # image marks: %d\n',nnz(s0.vis));

fprintf(fid,'\nOP ray count histogram (n, count):\n');
rayCount=ihist(sum(s0.vis,2)+1);
[n,~,count]=find(rayCount);
for i=1:length(n)
    fprintf(fid,'%2d: %7d\n',n(i)-1,count(i));
end

rayAng=angles(s0,'Computing OP ray angles')*180/pi;
fprintf(fid,'\nSmallest OP ray angles (angle, ID, cameras)\n');
[ang,i]=sort(rayAng);
for j=1:nnz(ang<ang(min(3,end))*1.1+0.1)
    camVis=find(s0.vis(i(j),:));
    s=sprintf('%4d ',camVis);
    fprintf(fid,'%5.2f: %6d (%s)\n',ang(j),s0.OPid(i(j)),s(1:end-1));
end

fprintf(fid,'\nCamera lowest ray count (count, camNo):\n');
[camCount,i]=sort(full(sum(s0.vis,1)));
for j=1:nnz(camCount<max(10,camCount(min(3,end))*1.25+1))
    fprintf(fid,'%3d: %4d\n',camCount(j),i(j));
end

fprintf(fid,'\nCamera ray count histogram (nPts, nCams):\n');
[n,edges]=histcounts(sum(s0.vis,1),10);
edges(end)=edges(end)+1;
for i=1:length(n)
    fprintf(fid,'%4d-%4d: %3d\n',edges(i),edges(i+1)-1,n(i));
end

camRayAng=camangles(s0,'Computing camera ray angles')*180/pi;
fprintf(fid,'\nSmallest camera ray angles (angle, cam ID, nPts)\n');
[ang,i]=sort(camRayAng);
for j=1:nnz(ang<ang(min(5,end))*1.1+0.1)
    nPts=nnz(s0.vis(:,i(j)));
    fprintf(fid,'%5.2f: %6d (%d)\n',ang(j),i(j),nPts);
end

if fid>2, fclose(fid); end

if nargout>0, ppsz=psz; end
