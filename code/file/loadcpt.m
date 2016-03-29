function [id,XYZ,sigma,name]=loadcpt(fName)
%LOADCPT Load control points from text file. 
%
%   [ID,XYZ,SIGMA,NAME]=loadcpt(FNAME) loads control point information
%   from the text file FNAME. The returned information is the
%   N-vector ID with point ids, the 3-by-N array XYZ with
%   coordinates, the 3-by-N array SIGMA with coordinate standard
%   deviations, and the N-cell array NAME with point names.  If no
%   standard deviations are given, SIGMA=[].

% $Id$

[fid,msg]=fopen(fName,'rt');
if fid<0
    error('Could not open %s for reading: %s.',fName,msg);
end

id=zeros(0,1);
XYZ=zeros(3,0);
sigma=zeros(3,0);
name=cell(0);

while ~feof(fid)
    % Read one line.
    s=fgets(fid);
    % Remove leading blanks.
    s=fliplr(deblank(fliplr(s)));
    % Skip if starts with #
    if ~isempty(s) && s(1)=='#'
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
        XYZ(:,end+1)=a;
      case 6
        XYZ(:,end+1)=a(1:3);
        sigma(:,end+1)=a(4:6);
      otherwise
        error('Bad number of items on CP line.');
    end
end

fclose(fid);
