function s=loadpmreport(fileName)
%LOADPMREPORT Parse PM processing report file.
%
%   S=LOADPMREPORT(FILE) parses the PhotoModeler processing report
%   file FILE. The returned structure S contains the following
%   fields:
%   
%   projName      - project name, string.
%   runDate       - Processing timestamp, struct with fields
%       string        - date as string,
%       num           - parsed date as datenum.
%   pmVersion     - PhotoModeler version, string.
%   status        - success/failure status, logical.
%   procOpts      - processing options, struct with logical fields
%       orient        - orientation on/off,
%       globalOpt     - global optimization on/off,
%       calibration   - calibration on/off,
%       constraints   - constraints on/off.
%   totError      - iteration information, struct with fields
%       numStages     - number of processing stages, integer,
%       numIters      - number of iterations in last stage, integer,
%       firstErr      - first error, float,
%       lastErr       - last error, float.
%   imNames       - image file names, N-cell string array
%   EO            - estimated EO values [Xc;Yc;Zc;Omega;Phi;Kappa],
%                   6-by-N array. Returned positions are in project units,
%                   angles in radians.
%   EOstd         - estimated EO standard deviations, 6-by-N array
%   EOcorr        - estimated high EO correlations, 6N-by-6N sparse array
%   imageCount    - struct with integer fields
%       total         - total number,
%       bad           - number of images classified as bad,
%       weak          - number of images classified as weak,
%       ok            - number of images classified as ok,
%       oriented      - number of oriented images,
%       invCam        - number of images with inverse camera flag set.
%   cameras       - struct array with fields
%       name          - camera name, string,
%       calibrated    - calibrated status, logical,
%       usedInImages  - number of image using this camera, integer,
%       coverage      - average image point coverage, float percent.
%   ptsUncalibrated - struct with points outside calibrated camera
%                     region with fields
%       ptId            - point id,
%       imNum           - image number.
%   markPtResiduals - struct with fields
%       overallRMS      - overall 2D residual RMS
%       markMax2dRMS,
%       markMin2dRMS    - max/min 2D per-mark-pt RMS, struct with fields
%           rms             - 2D rms value, float
%           id              - id of mark point with max/min RMS
%           imNo            - image number where the max/min point was measured
%       objMax2dRMS,
%       objMin2dRMS     - max/min 2D per-obj-pt RMS, struct with fields
%           rms             - 2D rms value, float
%           id              - id of object point with max/min RMS
%   tightness       - max/min tightness, struct with fields
%       max             - maximum tightness value, float   
%       maxId           - id of point with maximum tightness, integer
%       min             - minimum tightness value, float
%       minId           - id of point with minimum tightness, integer
%   ptPrecision     - point precision statistics, struct with fields
%       overall3DRMS    - overall RMS vector length, float
%       max3Dvector     - maximum 3D vector length, float
%       max3DvectorId   - id of point with maximum 3D vector length, integer
%       min3Dvector     - minimum 3D vector length, float
%       min3DvectorId   - id of point with minimum 3D vector length, integer
%       max             - 3-by-N float with X/Y/Z maximum values,
%       min             - 3-by-N float with X/Y/Z minimum values.
%   ptAngles        - ray intersection angle statistics, struct with fields
%       avg             - average angle (degrees), float,
%       max             - maximum angle (degrees), float,
%       maxId           - id of point with maximum angle, integer,
%       min             - minimum angle (degrees), float,
%       minId           - id of point with minimum angle, integer,


[fid,msg]=fopen(fileName,'rt');
if fid<0
    error('Failed to open file %s: %s.\n',fileName,msg);
end

projName=ReadProjName(fid);
runDate=ReadRunDate(fid);
pmVersion=ReadPMVersion(fid);
status=ReadProcessingStatus(fid);
procOpts=ReadProcessingOptions(fid);
totError=ReadTotalError(fid);
[imNames,EO,EOstd,EOcorr,l]=ReadEO(fid);
imageCount=ReadPhotoCount(fid,l);
[cameras,l]=ReadCameras(fid);
[ptsUncalibrated,l]=ReadUncalibratedPts(fid,l);
markPtResiduals=ReadMarkPtResiduals(fid,l);
tightness=ReadTightness(fid);
ptPrecision=ReadPrecision(fid);
ptAngles=ReadAngles(fid);

fclose(fid);

s=struct('projName',projName,'runDate',runDate,'pmVersion',pmVersion,...
         'status',status,'procOpts',procOpts,'totError',totError,...
         'imNames',{imNames},'EO',EO,'EOstd',EOstd,'EOcorr',EOcorr,...
         'imageCount',imageCount,'cameras',cameras,...
         'ptsUncalibrated',ptsUncalibrated,'markPtResiduals',markPtResiduals,...
         'tightness',tightness,'ptPrecision',ptPrecision,'ptAngles',ptAngles);


% Read lines from the opened text file fid until eof(fid) or the line
% matches the regular expression pat. The pattern should contain a
% grouping operator (paranthesis). The matching string inside the
% parantheses are returned as match. If supplied, the match is tried
% first on the string s. No match before eof(fid) returns ''.
%
% The matching is performed case INsensitive.
%
% Example:
%     pat='^\s*Project Name:\s*(\S*)\s*$'
%     will return the sequence of non-whitespace characters after the
%     'Project Name:' prefix.
function match=ReadUntilMatch(fid,pat,s)

if nargin<3, s=fgetl(fid); end
    
while ~feof(fid)
    % Try to find a match.
    t=regexpi(s,pat,'tokens');
    if ~isempty(t)
        % Found a match, return it.
        match=t{1};
        if isempty(match)
            match='';
        else
            match=match{1};
        end
        return;
    else
        % Otherwise, get next line and try again.
        s=fgetl(fid);
    end
end
match='';


% Parse the PM report file until the project name is found.
function projName=ReadProjName(fid)

projName=ReadUntilMatch(fid,'^\s*Project Name:\s*(\S+)\s*$');
    

% Parse the PM report file until the 'Last Processing Attempt' is found.
function runDate=ReadRunDate(fid)

% Get the full Last Processing Attempt string as a string.
runDateString=ReadUntilMatch(fid,'\s*Last Processing Attempt:\s*(.*)\s*$');

% Try to parse the string, sans leading weekday name. May fail for
% date strings are written in non-English locales.
try
    % Remove first word
    if any(isspace(runDateString))
        num=datenum(runDateString(find(isspace(runDateString),1)+1:end));
    else
        num=nan;
    end
catch
    % Failed to parse a datenum.
    num=nan;
end

runDate=struct('string',runDateString,'num',num);


% Parse the PM report file until the PM version string is found.
function ver=ReadPMVersion(fid)

ver=ReadUntilMatch(fid,'^\s*Version:\s*(.*)\s*$');
    

% Parse the PM report file until the processing status string is found.
function status=ReadProcessingStatus(fid)

s=ReadUntilMatch(fid,'^\s*Status:\s*(.*)\s*$');

status=strcmpi(s,'successful');


% Parse the PM report file for the Processing Options.
function opts=ReadProcessingOptions(fid)

% Skip until 'Processing Options'
ReadUntilMatch(fid,'^\s*Processing Options\s*$');

orientStr=ReadUntilMatch(fid,'^\s*Orientation:\s*(\S+)\s*$');
orient=strcmpi(orientStr,'on');
globalOptStr=ReadUntilMatch(fid,'^\s*Global Optimization:\s*(\S+)\s*$');
globalOpt=strcmpi(globalOptStr,'on');
calibrationStr=ReadUntilMatch(fid,'^\s*Calibration:\s*(\S+)\s*$');
calibration=strcmpi(calibrationStr,'on');
constraintsStr=ReadUntilMatch(fid,'^\s*Constraints:\s*(\S+)\s*$');
constraints=strcmpi(constraintsStr,'on');

opts=struct('orient',orient,'globalOpt',globalOpt,'calibration',calibration,...
            'constraints',constraints);


% Parse the PM report file for the Total Error.
function opts=ReadTotalError(fid)

% Skip until 'Total Error'
ReadUntilMatch(fid,'^\s*Total Error\s*$');

numIters=str2double(ReadUntilMatch(fid,'^\s*Number of Processing Iterations:\s*(\S+)\s*$'));
numStages=str2double(ReadUntilMatch(fid,'^\s*Number of Processing Stages:\s*(\S+)\s*$'));
firstErr=str2double(ReadUntilMatch(fid,'^\s*First Error:\s*(\S+)\s*$'));
lastErr=str2double(ReadUntilMatch(fid,'^\s*Last Error:\s*(\S+)\s*$'));

opts=struct('numIters',numIters,'numStages',numStages,'firstErr',firstErr,...
            'lastErr',lastErr);



% Return photo number imNum and image file name imName if line l is
% a photo header matching 'Photo X: IMFILENAME'. Otherwise returns
% empty arrays.
function [imNum,imName]=ImageNumName(l)

imNum=[];
imName='';

pat='^\s*Photo\s*(\d+)\s*:\s*(\S*)\s*$';
t=regexpi(l,pat,'tokens');

if length(t)==1 && length(t{1})==2
    imNum=str2double(t{1}{1});
    imName=t{1}{2};
end


% Parse correlation information from line l. Return percentage
% value in double val and parameter name in string param if a
% match, otherwise return empty arrays.
function [val,param]=CorrelationValues(l)
    
val=[];
param='';

pat='^\s*Correlations over[^:]+:\s*(.*)\s*$';
t=regexpi(l,pat,'tokens');
if length(t)==1 && iscell(t) && length(t{1})==1
    % Correlation string to parse.
    s=t{1}{1};
    pat2='^(\w+):([^%]+)%\s*(.*)$';
    t2=regexpi(s,pat2,'tokens');
    if length(t2)==1 && iscell(t2)
        val=str2double(t2{1}{2});
        param=t2{1}{1};
    end
    if ~isempty(t2{1}{3})
        % TODO: Get PM report file with multiple high correlations
        % are parse all correlations.
        warning('Only first correlation parsed: line=%s.\n',l)
    end
end


% Parse coverage information from line l. Return percentage
% value in double val if a match, otherwise return empty array.
function val=AveragePhotoPointCoverage(l)
    
val=[];

pat='^\s*Average Photo Point Coverage:\s*(\d+)%\s*$';
t=regexpi(l,pat,'tokens');
if length(t)==1 && iscell(t) && length(t{1})==1
    % Correlation string to parse.
    val=str2double(t{1}{1});
end


% Parse one EO element. l is first unprocessed line before/after the
% call.
function [val,dev,corrVal,corrName,l]=ReadOneEOVal(fid,l,name,unit)

val=[];
dev=[];
corrVal=[];
corrName='';

% Skip until name is found. Should be in l on first line from file.
pat=['^\s*',name,'\s*$'];
ReadUntilMatch(fid,pat,l);

% Next line should contain the value.
l=fgetl(fid);
pat=['^\s*Value:\s*(\S+)\s*',unit,'\s*$'];
t=regexpi(l,pat,'tokens');
if length(t)==1 && iscell(t) && length(t{1})==1
    val=str2double(t{1}{1});
else
    return;
end

% Next line should contain the standard deviation.
l=fgetl(fid);
pat=['^\s*Deviation:[^:]*:\s*(\S+)\s*',unit,'\s*$'];
t=regexpi(l,pat,'tokens');
if length(t)==1 && iscell(t) && length(t{1})==1
    dev=str2double(t{1}{1});
else
    return;
end

% Next line MAY contain correlations.
l=fgetl(fid);
[corrVal,corrName]=CorrelationValues(l);

if ~isempty(corrVal)
    % If we found correlations, consider the line processed.
    l='';
end

% Read image number, image name + values, std, and correlation for
% one image. Convert angles to radians.
function [EO,EOstd,EOcorrIJV,l]=ReadOneEO(fid)

fileOrderNames={'Omega','Phi','Kappa','Xc','Yc','Zc'};
fileOrderUnits={'deg','deg','deg','','',''};
dataOrderName={'X','Y','Z','Omega','Phi','Kappa'};
dataOrderNum=[4:6,1:3];

EO=nan(6,1);
EOstd=nan(6,1);
EOcorrIJV=zeros(0,3);

l='';

for i=1:length(fileOrderNames)
    [val,dev,corrVal,corrName,l]=ReadOneEOVal(fid,l,fileOrderNames{i},...
                                              fileOrderUnits{i});
    EO(dataOrderNum(i))=val;
    EOstd(dataOrderNum(i))=dev;
    if ~isempty(corrVal)
        EOcorrIJV=[EOcorrIJV;
                   dataOrderNum(i),find(strcmp(corrName,dataOrderName)),corrVal];
    end
end

% Convert angles to radians.
EO(4:6)=EO(4:6)*pi/180;
EOstd(4:6)=EOstd(4:6)*pi/180;


% Parse the PM report file for EO pos, std, corr, image names. Also
% return first unprocessed line.
function [imNames,EO,EOstd,EOcorr,l]=ReadEO(fid);

% Skip until 'Photograph Standard Deviations'
ReadUntilMatch(fid,'^\s*Photograph Standard Deviations\s*$');

imNames=cell(1,0);
EO=zeros(6,0);
EOstd=zeros(6,0);
EOcorrIJV=zeros(0,3);

% Get first line.
l=fgetl(fid);
% Check if line matches Photo N: NAME
[imNum,imName]=ImageNumName(l);
while isscalar(imNum)
    imNames{imNum}=imName;
    % Read one block of EO info.
    [EO(:,imNum),EOstd(:,imNum),newEOcorrIJV,l]=ReadOneEO(fid);
    EOcorrIJV=[EOcorrIJV;newEOcorrIJV+...
               repmat([[1,1]*(imNum-1)*6,0],size(newEOcorrIJV,1),1)];
    if isempty(l)
        l=fgetl(fid);
    end
    [imNum,imName]=ImageNumName(l);
end

EOcorr=sparse(EOcorrIJV(:,1),EOcorrIJV(:,2),EOcorrIJV(:,3),...
              size(EO,2)*6,size(EO,2)*6);


% Parse the PM report file for the photo count.
function s=ReadPhotoCount(fid,l)

% Skip until 'Photographs'
ReadUntilMatch(fid,'^\s*Photographs\s*$');

total=str2double(ReadUntilMatch(fid,'^\s*Total Number:\s*(\S+)\s*$'));
bad=str2double(ReadUntilMatch(fid,'^\s*Bad Photos:\s*(\S+)\s*$'));
weak=str2double(ReadUntilMatch(fid,'^\s*Weak Photos:\s*(\S+)\s*$'));
ok=str2double(ReadUntilMatch(fid,'^\s*Ok Photos:\s*(\S+)\s*$'));
oriented=str2double(ReadUntilMatch(fid,'^\s*Number Oriented:\s*(\S+)\s*$'));
invCam=str2double(ReadUntilMatch(fid,['^\s*Number with inverse camera ' ...
                    'flags set:\s*(\S+)\s*$']));

s=struct('total',total,'bad',bad,'weak',weak,'ok',ok,'oriented',oriented,...
         'invCam',invCam);


% Return camera number camNum and camera name camName if line l is a
% camera header matching 'CameraX: Camera name'. Otherwise returns
% empty arrays.
function [camNum,camName]=CameraNumName(l)

camNum=[];
camName='';

pat='^\s*Camera(\d+)\s*:\s*(.*)\s*$';
t=regexpi(l,pat,'tokens');

if length(t)==1 && length(t{1})==2
    camNum=str2double(t{1}{1});
    camName=t{1}{2};
end


% Parse the PM report file for the camera info. l is the first
% unprocessed line.
function [s,l]=ReadCameras(fid)

% Skip until 'Cameras'
ReadUntilMatch(fid,'^\s*Cameras\s*$');

camNames=cell(1,0);
isCalibrated=false(1,0);
usedInImages=zeros(1,0);
coverage=zeros(1,0);

% Get first line.
l=fgetl(fid);
% Check if line matches CameraN: NAME
[camNum,camName]=CameraNumName(l);
while isscalar(camNum)
    camNames{camNum}=camName;
    isCalibrated(camNum)=strcmpi(...
        ReadUntilMatch(fid,'^\s*Calibration:\s*(\S*)\s*$'),'yes');
    usedInImages(camNum)=str2double(...
        ReadUntilMatch(fid,'^\s*Number of photos using camera:\s*(\d+)\s*$'));
    % Get next lnie.
    l=fgetl(fid);
    % Parse coverage if it is there.
    val=AveragePhotoPointCoverage(l);
    if ~isempty(val)
        % Detected. Store value and read next line.
        coverage(camNum)=val;
        l=fgetl(fid);
    end
    [camNum,camName]=CameraNumName(l);
end

s=struct('name',{camNames},'calibrated',isCalibrated,...
         'usedInImages',usedInImages,'coverage',coverage);


% Parse list of points outside calibrated region.
function [s,l]=ReadUncalibratedPts(fid,l);

ptId=zeros(1,0);
imNum=zeros(1,0);

% Skip until 'Photo Coverage'
ReadUntilMatch(fid,'^\s*Photo Coverage\s*$',l);

% Skip until 'Photo Coverage'
ReadUntilMatch(fid,'^\s*Referenced points outside of the camera''s calibrated coverage region:\s*$');

pat='^\s*Point Marking Residuals\s*$';
l=fgetl(fid);
t=regexpi(l,pat,'tokens');
if isempty(t) && ~feof(fid)
    % TODO: Get file with uncalibrated points and fix this function
    % to handle it.
    warning('Uncalibrated points: line=%s.\n',l);
    l=fgetl(fid);
    t=regexpi(l,pat,'tokens');
end

s=struct('ptId',ptId,'imNum',imNum);

% Parse point marking residual statistics.
function markPtResiduals=ReadMarkPtResiduals(fid,l);

% Skip until 'Point Marking Residuals'
ReadUntilMatch(fid,'^\s*Point Marking Residuals\s*$',l);

overallRMS=str2double(...
    ReadUntilMatch(fid,'^\s*Overall RMS\s*:\s*(\S+)\s*pixels\s*$'));

max=str2double(ReadUntilMatch(fid,'^\s*Maximum\s*:\s*(\S+)\s*pixels\s*$'));
l=fgetl(fid);
pat='^\s*Point\s*(\d+)\s*on Photo\s*(\d+)\s*$';
t=regexpi(l,pat,'tokens');
if length(t)==1 && iscell(t) && length(t{1})==2
    maxId=str2double(t{1}{1});
    maxImNum=str2double(t{1}{2});
else
    maxId=nan;
    maxImNum=nan;
end
markMax2dRMS=struct('rms',max,'id',maxId,'imNo',maxImNum);

min=str2double(ReadUntilMatch(fid,'^\s*Minimum\s*:\s*(\S+)\s*pixels\s*$'));
l=fgetl(fid);
pat='^\s*Point\s*(\d+)\s*on Photo\s*(\d+)\s*$';
t=regexpi(l,pat,'tokens');
if length(t)==1 && iscell(t) && length(t{1})==2
    minId=str2double(t{1}{1});
    minImNum=str2double(t{1}{2});
else
    minId=nan;
    minImNum=nan;
end
markMin2dRMS=struct('rms',min,'id',minId,'imNo',minImNum);

max=str2double(ReadUntilMatch(fid,'^\s*Maximum RMS\s*:\s*(\S+)\s*pixels\s*$'));
maxId=str2double(ReadUntilMatch(fid,'^\s*Point\s*(\d+)\s*$'));
objMax2dRMS=struct('rms',max,'id',maxId);

min=str2double(ReadUntilMatch(fid,'^\s*Minimum RMS\s*:\s*(\S+)\s*pixels\s*$'));
minId=str2double(ReadUntilMatch(fid,'^\s*Point\s*(\d+)\s*$'));
objMin2dRMS=struct('rms',min,'id',minId);

markPtResiduals=struct('overallRMS',overallRMS,'markMax2dRMS',markMax2dRMS,...
                       'markMin2dRMS',markMin2dRMS,...
                       'objMax2dRMS',objMax2dRMS,...
                       'objMin2dRMS',objMin2dRMS);


% Parse tightness.
function tightness=ReadTightness(fid);

% Skip until 'Point Tightness'
ReadUntilMatch(fid,'^\s*Point Tightness\s*$');

max=str2double(ReadUntilMatch(fid,'^\s*Maximum\s*:\s*(\S+)\s*\S+\s*$'));
maxId=str2double(ReadUntilMatch(fid,'^\s*Point\s*(\d+)\s*$'));

min=str2double(ReadUntilMatch(fid,'^\s*Minimum\s*:\s*(\S+)\s*\S+\s*$'));
minId=str2double(ReadUntilMatch(fid,'^\s*Point\s*(\d+)\s*$'));

tightness=struct('max',max,'maxId',maxId,'min',min,'minId',minId);


% Parse point precisions.
function ptPrecision=ReadPrecision(fid);

% Skip until 'Point Precisions'
ReadUntilMatch(fid,'^\s*Point Precisions\s*$');

overall3DRMS=str2double(ReadUntilMatch(fid,['^\s*Overall RMS Vector ' ...
                    'Length\s*:\s*(\S+)\s*\S+\s*$']));

max3Dvector=str2double(ReadUntilMatch(fid,'^\s*Maximum Vector Length\s*:\s*(\S+)\s*\S+\s*$'));
max3DvectorId=str2double(ReadUntilMatch(fid,'^\s*Point\s*(\d+)\s*$'));

min3Dvector=str2double(ReadUntilMatch(fid,'^\s*Minimum Vector Length\s*:\s*(\S+)\s*\S+\s*$'));
min3DvectorId=str2double(ReadUntilMatch(fid,'^\s*Point\s*(\d+)\s*$'));

maxX=str2double(ReadUntilMatch(fid,'^\s*Maximum X\s*:\s*(\S+)\s*\S+\s*$'));
maxY=str2double(ReadUntilMatch(fid,'^\s*Maximum Y\s*:\s*(\S+)\s*\S+\s*$'));
maxZ=str2double(ReadUntilMatch(fid,'^\s*Maximum Z\s*:\s*(\S+)\s*\S+\s*$'));
minX=str2double(ReadUntilMatch(fid,'^\s*Minimum X\s*:\s*(\S+)\s*\S+\s*$'));
minY=str2double(ReadUntilMatch(fid,'^\s*Minimum Y\s*:\s*(\S+)\s*\S+\s*$'));
minZ=str2double(ReadUntilMatch(fid,'^\s*Minimum Z\s*:\s*(\S+)\s*\S+\s*$'));

ptPrecision=struct('overall3DRMS',overall3DRMS,'max3Dvector',max3Dvector,...
                   'max3DvectorId',max3DvectorId,'min3Dvector',min3Dvector,...
                   'min3DvectorId',min3DvectorId,'max',[maxX;maxY;maxZ],...
                   'min',[minX;minY;minZ]);

% Parse point angles.
function ptAngles=ReadAngles(fid);

% Skip until 'Point Angles'
ReadUntilMatch(fid,'^\s*Point Angles\s*$');

max=str2double(ReadUntilMatch(fid,'^\s*Maximum\s*:\s*(\S+)\s*\S+\s*$'));
maxId=str2double(ReadUntilMatch(fid,'^\s*Point\s*(\d+)\s*$'));

min=str2double(ReadUntilMatch(fid,'^\s*Minimum\s*:\s*(\S+)\s*\S+\s*$'));
minId=str2double(ReadUntilMatch(fid,'^\s*Point\s*(\d+)\s*$'));

avg=str2double(ReadUntilMatch(fid,'^\s*Average\s*:\s*(\S+)\s*\S+\s*$'));

ptAngles=struct('avg',avg,'max',max,'maxId',maxId,'min',min,'minId',minId);
