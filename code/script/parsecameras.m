function cams=parsecameras(cameras,baseDir)
%PARSECAMERAS Parse the input/cameras section of a DBAT XML file.
%
%   CAMS=PARSECAMERAS(CAMERAS,BASEDIR) parses the CAMERAS block of the
%   INPUT section of a DBAT XML script file. The result is returned
%   in the struct CAMS.
%
%   The camera data can either be present in the CAMERA field of
%   CAMERAS or in a DBAT camera XML file specified by the FILE
%   field.
%
%See also: LOADCAMERAS.

narginchk(2,2)
 
% Test for extra fields.
[ok,msg]=checkxmlfields(cameras,{'file','camera','c'},false);
if ~ok
    allFields=join(fieldnames(cameras),', ');
    error('DBAT XML input/cameras error: %s. Read fields are: %s.', ...
          msg,allFields{1});
end

if isfield(cameras,'file')
    % Cameras are stored in a file.
    file=cameras.file;
    
    % Check fields.
    [ok,msg]=checkxmlfields(file,'Text');
    if ~ok
        allFields=join(fieldnames(file),', ');
        error('DBAT XML input/cameras/file error: %s. Read fields are: %s.', ...
              msg,allFields{1});
    end
    
    fileName=parsepath(file.Text,baseDir);
    if ~exist(fileName,'file')
        warning('Cameras file %s does not exist',fileName);
    end

    cams=loadcameras(fileName);
elseif isfield(cameras,'camera')
    % Camera data is hardcoded in the XML file.
    cams=parsedbatxmlcamstruct(cameras.camera);
else
    % No camera data present.
    error('Missing camera data');
end

% Set id of single camera without id to 1.
if length(cams)==1 && isnan(cams{1}.Id)
    cams{1}.Id=1;
end

% Currently no support for multiple cameras with different camera
% units.
if length(unique(cellfun(@(x)x.Unit,cams,'UniformOutput',false)))>1
    error('Multiple cameras must have the same camera unit');
end

