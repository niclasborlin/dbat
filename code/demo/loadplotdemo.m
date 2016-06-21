function h=loadplotdemo(fName,alignCam,camSize)
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
    fName=fullfile(curDir,'data','phor2013','pmexports','roma-pmexport.txt');
    titleStr='Roma data (computed by PhotoModeler)';
    msg='Plotting the Roma camera network.';
    defAlignCam=1;
  case 'cam'
    fName=fullfile(curDir,'data','phor2013','pmexports','C4040Z-2272x1704-pmexport.txt');
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
plotnetwork(s,'title',titleStr,alignCmd{:},'camsize',camSize,...
            'axes',tagfigure(fName));
