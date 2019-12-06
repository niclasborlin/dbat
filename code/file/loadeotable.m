function pts=loadeotable(fName,fmt,sep,cmt)
%LOADEOTABLE Load external orientation data from text file.
%
%   EO=LOADEOTABLE(FNAME,FORMAT) loads external orientation information
%   from the text file FNAME. The format of the text lines are given
%   by FORMAT, see below. The data is returned in a struct EO with
%   fields
%       id       - 1-by-M array with id numbers.
%       pos      - 3-by-M array of positions.
%       ang      - 3-by-M array with omega-phi-kappa angles.
%       std      - 3-by-M array with prior position standard deviations.
%       angStd   - 3-by-M array with prior angular standard deviations.
%       fileName - string with FNAME.
%
%   Unspecified values are NaNs or the empty string, except std
%   that defaults to zero.
%
%   The string FORMAT contain the specification for the information on
%   each line in FNAME. FORMAT can contain any number of the following
%   strings, separated by commas:
%       id         - integer id for the control point,
%       ignored    - string/numerical field to be ignored,
%       x, y, z    - numeric fields with the x, y, or, z
%                    coordinates of the control point,
%       sx, sy, sz - numeric fields with the standard deviations of
%                    the individual coordinates,
%       sxy, sxyz  - numeric fields with the XY and XYZ standard
%                    deviations, respectively,
%       omega,phi,
%       kappa      - numeric fields with the omega, phi, or kappa, angles,
%       so, sp, sk - numeric fields with the standard deviations of
%                    the individual angles.
%   Whitespaces in the format string is ignored.
%
%   Blank lines and lines starting with a first non-whitespace
%   character '#' are treated as comments and are ignored. All other
%   lines are expected to match the format string, i.e., to contain
%   the data in the expected format.
%
%   Note that the label cannot contain the separator character. Use
%   PTS=LOADCTRLPTS(FNAME,FORMAT,SEP) to specify that another
%   separator character. The separator character should be used both
%   in the FORMAT string and the text file.
%
%   Use PTS=LOADCTRLPTS(FNAME,FORMAT,SEP,CMT) to specify another
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
knownParts={'id','ignored','x','y','z','sx','sy','sz','sxy', ...
            'sxyz','omega','phi','kappa','so','sp','sk','sang'};

if ~isempty(setdiff(fmtParts,knownParts))
    bad=join(setdiff(fmtParts,knownParts),', ');
    error('%s: Invalid format parts: %s',mfilename,bad{1});
end

id=zeros(1,0);
pos=zeros(3,0);
ang=zeros(3,0);
std=zeros(3,0);
angStd=zeros(3,0);

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

    ii=nan;
    n='';
    p=nan(3,1);
    s=zeros(3,1);
    a=nan(3,1);
    as=nan(3,1);
    for i=1:length(fmtParts)
        switch fmtParts{i}
          case 'id'
            ii=sscanf(parts{i},'%d');
          case 'label'
            n=parts{i};
          case {'x', 'y', 'z'}
            p(abs(fmtParts{i})-abs('x')+1)=sscanf(parts{i},'%f');
          case {'sx', 'sy', 'sz'}
            s(abs(fmtParts{i}(2))-abs('x')+1)=sscanf(parts{i},'%f');
          case 'sxy'
            s(1:2)=sscanf(parts{i},'%f');
          case 'sxyz'
            s(:)=sscanf(parts{i},'%f');
          case 'omega'
            a(1)=sscanf(parts{i},'%f');
          case 'phi'
            a(2)=sscanf(parts{i},'%f');
          case 'kappa'
            a(3)=sscanf(parts{i},'%f');
          case 'so'
            as(1)=sscanf(parts{i},'%f');
          case 'sp'
            as(2)=sscanf(parts{i},'%f');
          case 'sk'
            as(3)=sscanf(parts{i},'%f');
          case 'sang'
            as(:)=sscanf(parts{i},'%f');
          case 'ignored'
            % Do nothing
          otherwise
            % Should never happen
            error('%s: Invalid format part: ''%s''',mfilename,fmtParts{i});
        end
    end
    id(end+1)=ii;
    pos(:,end+1)=p;
    std(:,end+1)=s;
    ang(:,end+1)=a;
    angStd(:,end+1)=as;
end

fclose(fid);

pts=struct('id',id,'pos',pos,'ang',ang,'std',std,'angStd',angStd, ...
           'fileName',fName);
