function ims=parseimages(images,baseDir)
%PARSEIMAGES Parse the input/images section of a DBAT XML file.
%
%   IMS=PARSEIMAGES(IMAGES,BASEDIR) parses the IMAGES block of the
%   INPUT section of a DBAT XML script file. The result is returned
%   in the struct IMS.
%
%   The image data is loaded from the file specified in the
%   IMAGES.file.Text field and with the format specified by
%   IMAGES.file.Attributes.format.
%
%See also: LOADIMAGETABLE.

narginchk(2,2)

% Image list must be in a file.
[ok,msg]=checkxmlfields(images,'file');
if ~ok
    allFields=join(fieldnames(images),', ');
    error('DBAT XML input/images error: %s. Read fields are: %s.', ...
          msg,allFields{1});
end

file=images.file;

[ok,msg]=checkxmlfields(file,{'Attributes','Text'});
if ~ok
    allFields=join(fieldnames(file),', ');
    error('DBAT XML input/images/file error: %s. Read fields are: %s.', ...
          msg,allFields{1});
end

fileName=parsepath(file.Text,baseDir);
if ~exist(fileName,'file')
    warning('Images file %s does not exist',fileName);
end

[ok,msg]=checkxmlfields(file.Attributes,'format');
if ~ok
    allFields=join(fieldnames(file),', ');
    error('DBAT XML input/images/file attribute error: %s. Read fields are: %s.', ...
          msg,allFields{1});
end

% Extract format and load file using the format.
format=file.Attributes.format;

ims=loadimagetable(fileName,format);

% If all camera numbers are unspecified, use camera id=1 for all
% images.
if all(isnan(ims.cam))
    ims.cam(:)=1;
end
