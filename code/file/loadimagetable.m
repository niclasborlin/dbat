function ims=loadimagetable(fName,fmt,sep,cmt)
%LOADIMAGETABLE Load image table from text file.
%
%   IMS=LOADIMAGETABLE(FNAME,FORMAT) loads image information from the
%   text file FNAME. The format of the text lines are given by FORMAT,
%   see below. The data is returned in a struct IMS with fields
%       id       - 1-by-M array with id numbers.
%       name     - 1-by-M cell array with names.
%       path     - 1-by-M cell array with file paths.
%       fileName - string with FNAME.
%
%   Unspecified values are NaNs or the empty string, except name
%   that defaults to the last path component.
%
%   The string FORMAT contain the specification for the information on
%   each line in FNAME. FORMAT can contain any number of the following
%   strings, separated by commas:
%       id         - integer id for the control point,
%       label      - string with the name of the control point,
%       path       - string with the name of the control point,
%       dummy      - string/numerical field to be ignored,
%   Whitespaces in the format string is ignored.
%
%   Blank lines and lines starting with a first non-whitespace
%   character '#' are treated as comments and are ignored. All other
%   lines are expected to match the format string, i.e., to contain
%   the data in the expected format.
%
%   Note that any string fields cannot contain the separator
%   character. Use IMS=LOADIMAGETABLE(FNAME,FORMAT,SEP) to specify
%   that another separator character. The separator character should
%   be used both in the FORMAT string and the text file.
%
%   Use IMS=LOADIMAGETABLE(FNAME,FORMAT,SEP,CMT) to specify another
%   comment character.

if nargin<3, sep=','; end
if nargin<4, cmt='#'; end

[fid,msg]=fopen(fName,'rt');
if fid<0
    error('%s: Could not open %s for reading: %s.',mfilename,fName,msg);
end

% Parse the format string.
fmtParts=strip(strsplit(fmt,sep));

% Verify all parts are known.
knownParts={'id','label','path','dummy'};

if ~isempty(setdiff(fmtParts,knownParts))
    bad=join(setdiff(fmtParts,knownParts),', ');
    error('%s: Invalid format parts: %s',mfilename,bad{1});
end

id=zeros(1,0);
name=cell(1,0);
path=cell(1,0);

lineNo=0;

while ~feof(fid)
    % Read one line and clean from whitespace.
    s=strip(fgets(fid));
    lineNo=lineNo+1;
    % Skip blank or line starting with #.
    if isempty(s) || s(1)==cmt
        continue;
    end
    % Split line into parts
    parts=strip(strsplit(s,sep));
    
    if length(parts)~=length(fmtParts)
        error(['%s: %s, line %d: Wrong number of elements ',...
               '(got %d, expected %d)'], mfilename, fName, lineNo, ...
              length(parts), length(fmtParts));
    end

    % Determine if default labels should be used.
    defaultLabel=~ismember('label',fmtParts);
        
    ii=nan;
    n='';
    p='';
    for i=1:length(fmtParts)
        switch fmtParts{i}
          case 'id'
            ii=sscanf(parts{i},'%d');
          case 'label'
            n=parts{i};
          case 'path'
            p=parts{i};
            if defaultLabel
                [~,n1,n2]=fileparts(p);
                n=[n1,n2];
            end
          case 'dummy'
            % Do nothing
          otherwise
            % Should never happen
            error('%s: Invalid format part: ''%s''',mfilename,fmtParts{i});
        end
    end
    id(end+1)=ii;
    name{end+1}=n;
    path{end+1}=p;
end

fclose(fid);

ims=struct('id',id,'name',{name},'path',{path},'fileName',fName);
