function s=parseoutputfiles(s,outputFiles,docFile)
%PARSEOUTPUTFILES Parse output files section of a DBAT XML file.
%
%    S=PARSEOUTPUTFILES(FILES,DOCFILE,S) parses the output/files XML
%    block FILES from a DBAT XML script file and generates the
%    requested files. The string DOCFILE should contain the path name
%    of the XML file and is used to determine base directories.
%
%See also: PARSEINPUT.

narginchk(3,3);

knownFields={'Attributes','report','io','image_residuals','c'};
[ok,msg]=checkxmlfields(outputFiles,knownFields,false(size(knownFields)));
if ~ok, error('DBAT XML script output/files error: %s',msg); end

% Check if a base directory was specified.
[baseDir,rawBaseDir]=getattrbasedir(outputFiles,docFile);

% Warn if base dir does not exist.
if ~isempty(baseDir) && ~exist(baseDir,'dir')
    warning('DBAT output base directory %s (%s) does not exist!',baseDir, ...
            rawBaseDir);
end

files=fieldnames(outputFiles);

for i=1:length(files)
    file=files{i};
    switch file
      case 'Attributes'
        % Do nothing
      case 'report'
        WriteReportFile(s,outputFiles.report,baseDir);
      case 'io'
        warning('%s output file not implemented yet',file)
      case 'image_residuals'
        warning('%s output file not implemented yet',file)
      otherwise
        error('DBAT XML script output/file error: Unknown file %s',file);
    end
end


function WriteReportFile(s,report,baseDir)
% Write a report file to 

[ok,msg]=checkxmlfields(report,'file');
if ~ok, error('DBAT XML script output/files/report error: %s',msg); end

file=report.file;

% File name should be in the 'Text' field.
[ok,msg]=checkxmlfields(file,'Text');
if ~ok, error('DBAT XML script output/files/report/file error: %s',msg); end

% Get filename.
fileName=parsepath(file.Text,baseDir);

% Verify that directory exists.
if ~exist(fileparts(fileName),'dir')
    error('DBAT XML error: Output dir %s does not exist',fileparts(fileName));
end

s.bundle.info=bundle_cov(s,s.bundle.info,'prepare');

fprintf('\nBundle result file %s generated.\n',fileName);
s=bundle_result_file(s,s.bundle.info,fileName);
