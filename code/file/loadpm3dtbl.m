function pts=loadpm3dtbl(fName,isSmartTable)
%LOADPM3DTBL Load 3D point table exported from Photomodeler.
%
%   PTS=LOADPM3DTBL(FNAME) loads 3D point data from the exported 3D
%   point table in the file FNAME. The data is returned in a struct
%   PTS with fields
%       id   - N-array with id numbers.
%       name - N-cell array with names.
%       vis  - N-by-P visibility array.
%       pos  - 3-by-N array of estimated positions.
%       std  - 3-by-N array with posteriori standard deviations.
% 
%       N is the number of loaded points. P is the highest seen
%       image number.
%
%   PTS=LOADPM3DTBL(FNAME,TRUE) load 3D smart point data instead.
%
%   The function assumes the file is a save of a Photomodeler point
%   table or a SMART point table. Only lines starting with a
%   numeric ID are processed.


if nargin<2, isSmartTable=false; end

[fid,msg]=fopen(fName,'rt');
if fid<0
    error('Failed to open %s for reading: %s.\n',fName,msg); ...
end

digits='0':'9';

id=zeros(1,0);
name=cell(1,0);
vis=sparse(0,0);
pos=zeros(3,0);
std=zeros(3,0);
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

    if isSmartTable
        [ii,nn,vv,pp,ss]=ParseSmartPointLine(s);
    else
        [ii,nn,vv,pp,ss]=ParsePointLine(s);
    end
    
    id(end+1)=ii;
    name{end+1}=nn;
    vis(end+1,vv)=true;
    pos(:,end+1)=pp;
    std(:,end+1)=ss;
end

fclose(fid);

pts=struct('id',id,'name',{name},'vis',vis,'pos',pos,'std',std);


function [id,name,vis,pos,std]=ParsePointLine(s)
% Parse line from ordinary point table.

% Use regexp for first 3 components
% id, "name", "1,3,4",
v=regexp(s,'^\s*(\d+)\s*,\s*"([^"]*)"\s*,\s*"([^"]*)"\s*,(.*)$',...
         'tokens','once');
if length(v)~=4
    error('Unable to parse string: %s',s);
end
% Integer id
id=sscanf(v{1},'%d');
% Trim name and store as is.
name=strtrim(v{2});
% Images that the point was measured in (and used in the computations)
vis=sscanf(v{3},'%d,');
    
% Pos, std are floating point numbers.
[n,~,~,ni]=sscanf(v{4},'%g,');
if length(n)==3
    % Special handling of 'fixed' points.
    w=regexp(s,'(\s*fixed\s*,)+(.*)','tokens','once');
    n=[n;zeros(nnz(w{1}=='f'),1);sscanf(w{2},'%g,')];
end
if length(n)<6
    error('Could not parse position part of string: %s',v{4});
end
pos=n(1:3);
std=n(4:6);


function [id,name,vis,pos,std]=ParseSmartPointLine(s)
% Parse line from SMART point table.

% Use regexp for first 3 components

% id, "name"
v=regexp(s,'^\s*(\d+)\s*,\s*"([^"]*)"\s*,(.*)$','tokens','once');
if length(v)~=3
        error('Unable to parse SMART string: %s',s);
end
% Integer id
id=sscanf(v{1},'%d');
% Trim name and store as is.
name=strtrim(v{2});
    
% Pos, std are floating point numbers.
[n,~,~,ni]=sscanf(v{3},'%g,');
if length(n)<6
    error('Could not parse position part of SMART string: %s',v{3});
end

pos=n(1:3);
std=n(4:6);

% "1,3,4",
w=regexp(v{3}(ni:end),'\s*"([^"]*)"\s*,(.*)$','tokens','once');
% Images that the point was measured in (and used in the computations)
vis=sscanf(w{1},'%d,');
