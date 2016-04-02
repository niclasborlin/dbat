function pts=loadpm3dtbl(fName)
%LOADPM3DTBL Load 3D point table exported from Photomodeler.
%
%   PTS=LOADPM3DTBL(FNAME) loads 3D point data from the exported 3D
%   point table in the file FNAME. The data is returned in a struct
%   PTS with fields
%       id   - N-array with id numbers.
%       name - N-cell array with names.
%       vis  - N-by-P visibility array.
%       pos  - N-by-3 array of estimated positions.
%       std  - N-by-3 array with posteriori standard deviations.
% 
%       N is the number of loaded points. P is the highest seen
%       image number.
%
%   The function assumes the file is a save of a Photomodeler point
%   table. Only lines starting with a numeric ID are processed.

% $Id$

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
    s=strtrim(fgetl(fid));
    
    % Skip unless first non-white-space is a digit.
    if isempty(s) || ~any(s(1)==digits)
        continue
    end

    % Use regexp for first 3 components
    % id, "name", "1,3,4",
    v=regexp(s,'^\s*(\d+)\s*,\s*"([^"]*)"\s*,\s*"([^"]*)"\s*,(.*)$',...
             'tokens','once');
    if length(v)~=4
        error('Unable to parse string: %s',s);
    end
    % Integer id
    id(end+1)=sscanf(v{1},'%d');
    % Trim name and store as is.
    name{end+1}=strtrim(v{2});
    % Images that the point was measured in (and used in the
    % computations)
    j=sscanf(v{3},'%d,');
    vis(end+1,j)=true;
    
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
    pos(:,end+1)=n(1:3);
    std(:,end+1)=n(4:6);
end

fclose(fid);

pts=struct('id',id,'name',{name},'vis',vis,'pos',pos,'std',std);
