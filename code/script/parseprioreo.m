function camStations=parseprioreo(eo,baseDir)
%PARSEPRIOREO Parse the input/prior_eo section of a DBAT XML file.
%
%   CAMSTATIONS=PARSEPRIOREO(EO,BASEDIR) parses the prior_eo block of
%   the INPUT section of a DBAT XML script file. The result is
%   returned in the struct CAMSTATION.
%
%   The EO data is loaded from the file specified in the EO.file.Text
%   field and with the format specified by EO.file.Attributes.format.
%
%See also: LOADEOTABLE.

narginchk(2,2)

[ok,msg]=checkxmlfields(eo,{'file','c'},[true,false]);
if ~ok
    allFields=join(fieldnames(eo),', ');
    error('DBAT XML input/prior_eo error: %s. Read fields are: %s.', ...
          msg,allFields{1});
end

file=eo.file;

[ok,msg]=checkxmlfields(file,{'Attributes','Text'});
if ~ok
    allFields=join(fieldnames(file),', ');
    error('DBAT XML input/eo/file error: %s. Read fields are: %s.', ...
          msg,allFields{1});
end

fileName=parsepath(file.Text,baseDir);
if ~exist(fileName,'file')
    warning('prior_eo file %s does not exist',fileName);
end

[ok,msg]=checkxmlfields(file.Attributes,{'format','units'},[true,false]);
if ~ok
    allFields=join(fieldnames(file),', ');
    error('DBAT XML input/eo/file attribute error: %s. Read fields are: %s.', ...
          msg,allFields{1});
end

% Extract format and load file using the format.
format=file.Attributes.format;

camStations=loadeotable(fileName,format);

scaling=nan;
if isfield(file.Attributes,'units')
    switch (file.Attributes.units)
      case 'radian'
        scaling=1;
      case 'degrees'
        scaling=pi/180;
      case 'gon'
        scaling=pi/200;
      otherwise
        error('DBAT XML input/eo/file attribute error: Unknown unit %s.',...
              file.Attributes.units);
    end
end

if any(any(~isnan(camStations.ang)))
    % Angle information was specified. Verify that angle units were
    % specified.
    if ~isnan(scaling)
        camStations.ang=camStations.ang*scaling;
        camStations.angStd=camStations.angStd*scaling;
    else
        error('DBAT XML input/eo error: Angles read but no unit specified');
    end
end
