function pts=parseimagepts(imagePts,baseDir)
%PARSEIMAGEPTS Parse the input/image_pts section of a DBAT XML file.
%
%   PTS=PARSEIMAGEPTS(IMAGEPTS,BASEDIR) parses the IMAGE_PTS block of
%   the INPUT section of a DBAT XML script file. The result is
%   returned in the struct PTS.
%
%   The image point data is loaded from the file specified in the
%   IMAGEPTS.file.Text field and with the format specified by
%   IMAGEPTS.file.Attributes.format.
%
%See also: PARSEIMAGES, LOADIMAGEPTS.

narginchk(2,2)

% Image list must be in a file.
[ok,msg]=checkxmlfields(imagePts,'file');
if ~ok
    allFields=join(fieldnames(imagePts),', ');
    error('DBAT XML input/image_pts error: %s. Read fields are: %s.', ...
          msg,allFields{1});
end

file=imagePts.file;

[ok,msg]=checkxmlfields(file,{'Attributes','Text'});
if ~ok
    allFields=join(fieldnames(file),', ');
    error('DBAT XML input/image_pts/file error: %s. Read fields are: %s.', ...
          msg,allFields{1});
end

fileName=parsepath(file.Text,baseDir);
if ~exist(fileName,'file')
    warning('image_pts file %s does not exist',fileName);
end

[ok,msg]=checkxmlfields(file.Attributes,{'format','sxy','sx','sy'}, ...
                        [true,false(1,3)]);
if ~ok
    allFields=join(fieldnames(file),', ');
    error('DBAT XML input/image_pts/file attribute error: %s. Read fields are: %s.', ...
          msg,allFields{1});
end

% Extract format and load file using the format.
format=file.Attributes.format;

pts=loadimagepts(fileName,format);

% Apply standard deviations supplied via Attributes.
if isfield(file.Attributes,'sxy')
    sxy=sscanf(file.Attributes.sxy,'%f');
    pts.std(1:2,:)=sxy;
end

if isfield(file.Attributes,'sx')
    sx=sscanf(file.Attributes.sx,'%f');
    pts.std(1,:)=sx;
end

if isfield(file.Attributes,'sy')
    sy=sscanf(file.Attributes.sy,'%f');
    pts.std(2,:)=sy;
end
