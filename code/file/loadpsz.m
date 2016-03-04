function [p,obj,pts]=loadpsz(fName)
%LOADPSZ Load Photoscan 1.1.6 PSZ camera calibration file.
%
%

% Unpack zip archive into temporary folder.
tmpDir=tempname;
files=unzip(fName,tmpDir);

% Find doc.xml file.
docFile=files{strcmp(fullfile(tmpDir,'doc.xml'),files)};

% Read the doc.xml file.
obj=xml2struct2(docFile);

pFile=obj.document.chunks.chunk.frames.frame.point_cloud.points.Attributes.path;
fullpFile=files{strcmp(fullfile(tmpDir,pFile),files)};
[~,pts]=ply_read(fullpFile,'tri');

p=[];

for i=1:length(files)
    delete(files{i});
end
rmdir(tmpDir);
