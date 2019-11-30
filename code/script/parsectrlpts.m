function pts=parsectrlpts(ctrlPts,baseDir)
%PARSECTRLPTS Parse the input/ctrl_pts section of a DBAT XML file.
%
%   PTS=PARSECTRLPTS(CTRLPTS,BASEDIR) parses the ctrl_pts block of the
%   INPUT section of a DBAT XML script file. The result is returned
%   in the struct PTS.
%
%   The image data is loaded from the file specified in the
%   CTRLPTS.file.Text field and with the format specified by
%   CTRLPTS.file.Attributes.format.
%
%See also: LOADIMAGETABLE.

narginchk(2,2)

[ok,msg]=checkxmlfields(ctrlPts,{'file','filter','c'},[true,false,false]);
if ~ok
    allFields=join(fieldnames(ctrlPts),', ');
    error('DBAT XML input/ctrl_pts error: %s. Read fields are: %s.', ...
          msg,allFields{1});
end

file=ctrlPts.file;

[ok,msg]=checkxmlfields(file,{'Attributes','Text'});
if ~ok
    allFields=join(fieldnames(file),', ');
    error('DBAT XML input/ctrlPts/file error: %s. Read fields are: %s.', ...
          msg,allFields{1});
end

fileName=parsepath(file.Text,baseDir);
if ~exist(fileName,'file')
    warning('ctrl_pts file %s does not exist',fileName);
end

[ok,msg]=checkxmlfields(file.Attributes,'format');
if ~ok
    allFields=join(fieldnames(file),', ');
    error('DBAT XML input/ctrlPts/file attribute error: %s. Read fields are: %s.', ...
          msg,allFields{1});
end

% Extract format and load file using the format.
format=file.Attributes.format;

pts=loadctrlpts(fileName,format);

% Parse and execute any filtering commands
if isfield(ctrlPts,'filter')
    filters=ctrlPts.filter;
    if ~iscell(filters)
        filters={filters};
    end
    for i=1:length(filters)
        filter=filters{i};

        [ok,msg]=checkxmlfields(filter,{'Text','Attributes'},[true,true]);
        if ~ok
            allFields=join(fieldnames(filter),', ');
            error(['DBAT XML input/ctrl_pts/filter error: %s. Read ' ...
                   'fields are: %s.'], msg,allFields{1});
        end
       
        [ok,msg]=checkxmlfields(filter.Attributes,'id');
        if ~ok
            allFields=join(fieldnames(filter.Atrributes),', ');
            error(['DBAT XML input/ctrl_pts/filter/Attributes error: %s. Read ' ...
                   'fields are: %s.'], msg,allFields{1});
        end
        
        % TODO: Improve parsing to allow ranges
        id=sscanf(filter.Attributes.id,'%d,');
        
        switch strip(filter.Text)
          case 'remove'
            keep=~ismember(pts.id,id);
          case 'keep'
            keep=ismember(pts.id,id);
          otherwise
            error(['DBAT XML input/ctrl_pts/filter error: Unknown ' ...
                   'filter %s'],strip(filter.Text));
        end
        
        % Apply the actual filtering
        pts.id=pts.id(keep);
        pts.name=pts.name(keep);
        pts.pos=pts.pos(:,keep);
        pts.std=pts.std(:,keep);
    end
end
