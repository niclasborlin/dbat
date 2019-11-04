function ims=parseimages(images,baseDir,docFile)
%PARSEIMAGES Parse the input/images section of a DBAT XML file.
%
%   IMS=PARSEIMAGES(IMAGES,DOCFILE,BASEDIR) parses the IMAGES block of
%   the INPUT section of a DBAT XML script file. The result is
%   returned in the struct IMS. The BASEDIR is used as the image base
%   dir unless it is overridden by an image_base_dir attribute in
%   images. The string DOCFILE should contain the path name of the
%   source XML file.
%
%   The image data is loaded from the file specified in the
%   IMAGES.file.Text field and with the format specified by
%   IMAGES.file.Attributes.format.
%
%See also: LOADIMAGETABLE.

narginchk(3,3)

% Known fields are 'file' and optionally 'Attributes'.
[ok,msg]=checkxmlfields(images,{'file','Attributes'},[true,false]);
if ~ok
    allFields=join(fieldnames(images),', ');
    error('DBAT XML input/images error: %s. Read fields are: %s.', ...
          msg,allFields{1});
end

% Check if an image base directory was specified.
[imageBaseDir,rawImageBaseDir]=getattrbasedir(images,docFile,'image_base_dir');

if isempty(imageBaseDir)
    imageBaseDir=baseDir;
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

% Load the file.
ims=loadimagetable(fileName,format);

% Apply substitution to the image file names.
if ~isempty(imageBaseDir)
    for i=1:length(ims.path)
        ims.path{i}=parsepath(ims.path{i},imageBaseDir);
    end
end

% If all camera numbers are unspecified, use camera id=1 for all
% images.
if all(isnan(ims.cam))
    ims.cam(:)=1;
end
