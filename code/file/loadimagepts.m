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

% Convert parts to scanf strings
idIx=0;
imIx=0;
posIx=[0,0];
stdIx=[0,0];
fieldsToRead=0;

scanfStrings=cell(size(fmtParts));
for i=1:length(fmtParts)
    switch fmtParts{i}
      case 'id'
        scanfStrings{i}='%u';
        idIx=i;
        fieldsToRead=fieldsToRead+1;
      case 'im'
        scanfStrings{i}='%u';
        imIx=i;
        fieldsToRead=fieldsToRead+1;
      case 'ignored'
        scanfStrings{i}='%*f';
      case 'x'
        scanfStrings{i}='%g';
        posIx(1)=i;
        fieldsToRead=fieldsToRead+1;
      case 'y'
        scanfStrings{i}='%g';
        posIx(2)=i;
        fieldsToRead=fieldsToRead+1;
      case 'sx'
        scanfStrings{i}='%g';
        stdIx(1)=i;
        fieldsToRead=fieldsToRead+1;
      case 'sy'
        scanfStrings{i}='%g';
        stdIx(2)=i;
        fieldsToRead=fieldsToRead+1;
      case 'sxy'
        scanfStrings{i}='%g';
        stdIx(1:2)=i;
        fieldsToRead=fieldsToRead+1;
    end        
end
scanfString=join(scanfStrings,sep);
scanfString=scanfString{1};

% Initialize progress bar.
h=waitbar(0,'Loading image points');

% Initialize file pointer.
pos=ftell(fid);

% Read initial part that may contain comments.
while ~feof(fid)
    % Read one line and clean from whitespace.
    pos=ftell(fid);
    s=strip(fgets(fid));
    % Skip blank or line starting with #.
    if ~isempty(s) && s(1)~=cmt
        break;
    end
end

if ~ishandle(h)
    % Waitbar closed, abort.
    fclose(fid);
    err='Aborted by user';
    if nargout<2, error(err); else return; end
else
    % Update progressbar.
    waitbar(ftell(fid)/sz,h);
end

% Seek back to start of first data line.
fseek(fid,pos,'bof');

% Read 10 lines to estimate bytes/line.
for i=1:10
    fgetl(fid);
end

bytesPerLine=(ftell(fid)-pos)/10;

% Seek back to start of first data line.
fseek(fid,pos,'bof');

% Approximate number of data lines.
numDataLines=ceil((sz-pos)/bytesPerLine);

% Pre-allocate output data.
id=nan(1,numDataLines);
im=nan(1,numDataLines);
pos=nan(2,numDataLines);
std=zeros(2,numDataLines);

% Last occupied data element.
i=0;

% Read data blocks from file
blockSize=10000;
while ~feof(fid)
    [A,cnt]=fscanf(fid,scanfString,[fieldsToRead,blockSize]);
    rows=floor(cnt/fieldsToRead);
    if rem(cnt,fieldsToRead)~=0
        error('Bad data after reading %d rows + %d elements in file %s',...
              rows,rem(cnt,fieldToRead),fName);
    end
    if ~ishandle(h)
        % Waitbar closed, abort.
        fclose(fid);
        err='Aborted by user';
        if nargout<2, error(err); else return; end
    else
        % Update progressbar.
        waitbar(ftell(fid)/sz,h);
    end

    if idIx~=0
        id(i+(1:rows))=A(idIx,:);
    end
    if imIx~=0
        im(i+(1:rows))=A(imIx,:);
    end
    if all(posIx)~=0
        pos(:,i+(1:rows))=A(posIx,:);
    elseif posIx(1)~=0
        pos(1,i+(1:rows))=A(posIx(1),:);
    elseif posIx(2)~=0
        pos(2,i+(1:rows))=A(posIx(2),:);
    end
    if all(stdIx)~=0
        std(:,i+(1:rows))=A(stdIx,:);
    elseif stdIx(1)~=0
        std(1,i+(1:rows))=A(stdIx(1),:);
    elseif stdIx(2)~=0
        std(2,i+(1:rows))=A(stdIx(2),:);
    end
    
    i=i+rows;
end

if ishandle(h)
    % Update progressbar.
    waitbar(ftell(fid)/sz,h);
end

id=id(1:i);
im=im(1:i);
pos=pos(:,1:i);
std=std(:,1:i);

fclose(fid);

pts=struct('id',id,'im',im,'pos',pos,'std',std,'fileName',fName);

if ishandle(h), close(h); end
