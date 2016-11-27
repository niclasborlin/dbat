function pts=loadcpt(fName)
%LOADCPT Load control points from text file. 
%
%   PTS=loadcpt(FNAME) loads control point information
%   from the text file FNAME. The data is returned in a struct PTS
%   with fields
%       id   - M-array with id numbers.
%       name - M-cell array with names.
%       vis  - M-by-N visibility array
%       pos  - N-by-3 array of estimated positions
%       std  - N-by-3 array with posteriori standard deviations
%
%   Each line is expected to contain a comma-separated list of an
%   integer id, a name, X, Y, Z positions, and optionally X, Y, Z standard
%   devations. Blank lines and lines starting with # are ignored.


[fid,msg]=fopen(fName,'rt');
if fid<0
    error('Could not open %s for reading: %s.',fName,msg);
end

id=zeros(1,0);
pos=zeros(3,0);
std=zeros(3,0);
name=cell(1,0);

while ~feof(fid)
    % Read one line.
    s=fgets(fid);
    % Remove leading blanks.
    s=fliplr(deblank(fliplr(s)));
    % Skip blank or line starting with #.
    if isempty(s) || s(1)=='#'
        continue;
    end
    % Parse id.
    [i,n,msg,ni]=sscanf(s,'%d,');
    s=s(ni:end);
    % Parse name.
    [nn,n,msg,ni]=sscanf(s,'%s,');
    if ~isempty(nn) && nn(end)==','
        nn=nn(1:end-1);
    end
    s=s(ni:end);
    [a,n,msg,ni]=sscanf(s,'%g,');
    id(end+1)=i;
    name{end+1}=nn;
    switch n
      case 3
        pos(:,end+1)=a;
        std(:,end+1)=0;
      case 6
        pos(:,end+1)=a(1:3);
        std(:,end+1)=a(4:6);
      otherwise
        error('Bad number of items on CP line.');
    end
end

fclose(fid);

pts=struct('id',id,'name',{name},'pos',pos,'std',std);
