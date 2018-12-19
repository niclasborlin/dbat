function [h,s,prob]=loadplotdemo(fName,alignCam,camSize)
%LOADPLOTDEMO Load and plot PhotoModeler export file demo.
%
%   LOADPLOTDEMO(N), where N is the name of the PhotoModeler text
%   export file, load and plots the camera network and OP that were
%   computed by PhotoModeler.
%
%   LOADPLOTDEMO(N,ALIGNCAM), aligns the data set with camera ALIGNCAM
%   such that the camera +Z axis is 'up' instead of forward. Use the
%   default ALIGNCAM=0 to plot the raw coordinates instead.
%
%   LOADPLOTDEMO(N,ALIGNCAM,SZ) sets the width of the camera icon to
%   SZ in object units. SZ defaults to unity.
%
%   Use LOADPLOTDEMO('ROMA') or LOADPLOTDEMO('CAM') to load and plot
%   two test data sets shipped with DBAT. LOADPLOTDEMO without
%   arguments defaults to the ROMA data set.
%
%See also: PLOTNETWORK.

if nargin<1, fName='roma'; end

% Extract name of current directory.
curDir=fileparts(mfilename('fullpath'));

% Default cameras size.
defCamSize=1;

% Alignment camera number. Default to no alignment.
defAlignCam=0;

switch lower(fName)
  case 'roma'
    fName=fullfile(curDir,'data','dbat','pmexports','roma-pmexport.txt');
    titleStr='Roma data (computed by PhotoModeler)';
    msg='Plotting the Roma camera network.';
    defAlignCam=1;
  case 'cam'
    fName=fullfile(curDir,'data','dbat','pmexports','camcal-pmexport.txt');
    titleStr='Camera calibration data set (computed by Photomodeler)';
    msg='Plotting the camera calibration network.';
    defCamSize=0.15;
  otherwise
    if exist(fName,'file')~=2
        error('File %s does not exist.\n',fName);
    end
    msg=sprintf('Plotting %s.',fName);
    % Use file name as title.
    [~,n,e]=fileparts(fName);
    titleStr=[n,e];
end

% Use defaults unless the corresponding parameters were supplied.
if nargin<2, alignCam=defAlignCam; end
if nargin<3, camSize=defCamSize; end

if alignCam==0
    alignCmd={};
else
    alignCmd={'align',alignCam,'trans','up'};
end

% Load data...
fprintf('Loading data file %s...',fName);
prob=loadpm(fName);
disp('done.')

% ...convert to useful struct...
s=prob2dbatstruct(prob);

% ...plot it.
disp(msg);
if alignCam==0
    disp('Using raw coordinates.');
else
    fprintf('Aligning with camera %d.\n',alignCam);
end
hh=plotnetwork(s,'title',titleStr,alignCmd{:},'camsize',camSize,...
              'axes',tagfigure(fName));

imName='';
imNo=1;
% Check if image files exist.
isAbsPath=~isempty(s.proj.imDir) && ismember(s.proj.imDir(1),'\\/') || ...
          length(s.proj.imDir)>1 && s.proj.imDir(2)==':';
if ~isAbsPath && exist(fullfile(curDir,s.proj.imDir),'dir')
    % Expand path relative to current dir for this file.
    s.proj.imDir=fullfile(curDir,s.proj.imDir);
end
if exist(s.proj.imDir,'dir')
    % Handle both original-case and lower-case file names.
    imNames={s.EO.name{imNo},lower(s.EO.name{imNo}),upper(s.EO.name{imNo})};    
    imNames=fullfile(s.proj.imDir,imNames);
    imExist=cellfun(@(x)exist(x,'file')==2,imNames);
    if any(imExist)
        imName=imNames{find(imExist,1,'first')};
    end
else
    warning('Image directory %s does not exist.',s.proj.imDir);
end

if exist(imName,'file')
    fprintf('Plotting measurements on image %d.\n',imNo);
    imFig=tagfigure('image');
    hh=[hh;imshow(imName,'parent',gca(imFig))];
    pts=s.IP.val(:,s.IP.ix(s.IP.vis(:,imNo),imNo));
    ptsId=s.OP.id(s.IP.vis(:,imNo));
    isCtrl=s.prior.OP.isCtrl(s.IP.vis(:,imNo));
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

if nargout>0, h=hh; end
