function pts=loadpm2dtbl(fName)
%LOADPM2DTBL Load 2D point table exported from Photomodeler.
%
%   PTS=LOADPM2DTBL(FNAME) loads 2D point data from the exported 2D
%   point table in the file FNAME. The data is returned in a struct
%   PTS with fields
%       id   - N-array with id numbers.
%       imNo - N-array with image numbers.
%       pos  - 2-by-N array of measured positions.
%       res  - 2-by-N array of residuals.
% 
%       N is the number of loaded points.
%
%   The function assumes the file is a save of a 2D Point table or a
%   2D SmartPoint table. Only lines starting with a numeric ID are
%   processed.


[fid,msg]=fopen(fName,'rt');
if fid<0
    error('Failed to open %s for reading: %s.\n',fName,msg); ...
end

digits='0':'9';

id=zeros(1,0);
imNo=zeros(1,0);
pos=zeros(2,0);
res=zeros(2,0);
while ~feof(fid)
    % Read one line.
    sRaw=fgetl(fid);
    if ~ischar(sRaw)
        % End of file.
        break;
    end
    % Remove whitespace.
    s=strtrim(sRaw);
    
    % Skip unless first non-white-space is a digit.
    if isempty(s) || ~any(s(1)==digits)
        continue
    end

    [ii,mm,pp,rr]=ParsePointLine(s);
    
    id(end+1)=ii;
    imNo(end+1)=mm;
    pos(:,end+1)=pp;
    res(:,end+1)=rr;
end

fclose(fid);

pts=struct('id',id,'imNo',imNo,'pos',pos,'res',res);


function [id,imNo,pos,res]=ParsePointLine(s)
% Parse line from 2D point table.

% Id, imNo, x, y, resX, resY
v=sscanf(s,'%d,%d,%g,%g,%g,%g');

if length(v)~=6
    error('Unable to parse string: %s',s);
end

id=v(1);
imNo=v(2);
pos=v(3:4);
res=v(5:6);
