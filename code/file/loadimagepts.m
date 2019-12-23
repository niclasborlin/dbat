function pts=loadimagepts(fName,fmt,sep,cmt)
%LOADIMAGEPTS Load image points from text file.
%
%   PTS=LOADIMAGEPTS(FNAME,FORMAT) loads 2D image point information
%   from the text file FNAME. The format of the text lines are given
%   by FORMAT, see below. The data is returned in a struct PTS with
%   fields
%       id       - 1-by-M array with point id numbers.
%       im       - 1-by-M array with image id numbers.
%       pos      - 2-by-M array of positions.
%       std      - 2-by-M array with prior standard deviations.
%       fileName - string with FNAME.
%
%   Unspecified values are NaNs or the empty string, except std that
%   defaults to zero.
%
%   The string FORMAT contain the specification for the information on
%   each line in FNAME. FORMAT can contain any number of the following
%   strings, separated by commas:
%       id      - integer id for the point,
%       im      - integer id for the image number,
%       ignored - string/numerical field to be ignored,
%       x, y    - numeric fields with the x or y coordinates of the
%                 image point,
%       sx, sy  - numeric fields with the standard deviations of the
%                 individual coordinates, 
%       sxy     - numeric fields with the XY deviation.
%   Whitespaces in the format string is ignored.
%
%   Blank lines and lines starting with a first non-whitespace
%   character '#' are treated as comments and are ignored. All other
%   lines are expected to match the format string, i.e., to contain
%   the data in the expected format.
%
%   Use PTS=LOADIMAGEPTS(FNAME,FORMAT,SEP) to specify that another
%   separator character. The separator character should be used both
%   in the FORMAT string and the text file.
%
%   Use PTS=LOADIMAGEPTS(FNAME,FORMAT,SEP,CMT) to specify another
%   comment character.

if nargin<3, sep=','; end
if nargin<4, cmt='#'; end

[fid,msg]=fopen(fName,'rt');
if fid<0
    error('%s: Could not open %s for reading: %s.',mfilename,fName,msg);
end

% Get file size for progress bar.
if fseek(fid,0,'eof')~=0
	err='Failed to get file size';
	if nargout<2, error(err); else return; end
end
sz=ftell(fid);
% Rewind.
fseek(fid,0,'bof');

% Parse the format string.
fmtParts=strip(strsplit(fmt,sep));

% Verify all parts are known.
knownParts={'id','im','ignored','x','y','sx','sy','sxy'};

if ~isempty(setdiff(fmtParts,knownParts))
    bad=join(setdiff(fmtParts,knownParts),', ');
    error('%s: Invalid format parts: %s',mfilename,bad{1});
end

id=zeros(1,0);
im=zeros(1,0);
pos=zeros(2,0);
std=zeros(2,0);

lineNo=0;

% Initialize progress bar.
h=waitbar(0,'Loading image points');

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

    if rem(lineNo+1,1000)==0
        if ~ishandle(h)
            % Waitbar closed, abort.
            fclose(fid);
            err='Aborted by user';
            if nargout<2, error(err); else return; end
        else
            % Update progressbar every 1000 input lines.
            waitbar(ftell(fid)/sz,h);
        end
    end
    
    ii=nan;
    mm=nan;
    p=nan(2,1);
    s=zeros(2,1);
    for i=1:length(fmtParts)
        switch fmtParts{i}
          case 'id'
            ii=sscanf(parts{i},'%d');
          case 'im'
            mm=sscanf(parts{i},'%d');
          case {'x', 'y'}
            p(abs(fmtParts{i})-abs('x')+1)=sscanf(parts{i},'%f');
          case {'sx', 'sy'}
            s(abs(fmtParts{i}(2))-abs('x')+1)=sscanf(parts{i},'%f');
          case 'sxy'
            s(:)=sscanf(parts{i},'%f');
          case 'ignored'
            % Do nothing
          otherwise
            % Should never happen
            error('%s: Invalid format part: ''%s''',mfilename,fmtParts{i});
        end
    end
    id(end+1)=ii;
    im(end+1)=mm;
    pos(:,end+1)=p;
    std(:,end+1)=s;
end

if ishandle(h), waitbar(ftell(fid)/sz,h); end

fclose(fid);

pts=struct('id',id,'im',im,'pos',pos,'std',std,'fileName',fName);

if ishandle(h), close(h); end
