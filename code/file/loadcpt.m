function pts=loadcpt(fName,has)
%LOADCPT Load control points from text file. 
%
%   PTS=LOADCPT(FNAME) loads control point information
%   from the text file FNAME. The data is returned in a struct PTS
%   with fields
%       id       - 1-by-M array with id numbers, or NaN if no id.
%       name     - 1-by-M cell array with names, or blank if no name.
%       pos      - 3-by-M array of positions.
%       std      - 3-by-M array with prior standard deviations.
%       cov      - 3-by-3-by-M array with prior covariance matrices.
%       fileName - string with filename.
%
%   Each line is expected to contain a comma-separated list of an
%   integer id, a name, X, Y, Z positions, and optional std/covariance
%   information. Blank lines and lines starting with # are ignored.
%
%   PTS=LOADCPT(FNAME,[HAS_ID,HAS_NAME]), where HAS_ID and
%   HAS_NAME are logical, indicates whether the ID and/or NAME are
%   present in the file. At least one of HAS_ID and HAS_NAME must
%   be TRUE.
%
%   The std/covariance information may consist of:
%     - no value - pt is fixed,
%     - 1 value  - sigma_xyz,
%     - 2 values - sigma_xy, sigma_z,
%     - 3 values - sigma_x, sigma_y, sigma_z,
%     - 9 values - full covariance matrix.

if nargin<2, has=[true,true]; end

hasId=has(1);
hasName=has(2);

[fid,msg]=fopen(fName,'rt');
if fid<0
    error('Could not open %s for reading: %s.',fName,msg);
end

id=zeros(1,0);
pos=zeros(3,0);
std=zeros(3,0);
cov=zeros(3,3,0);
name=cell(1,0);

% True if any full covariance matrix is specified.
anyCov=false;

while ~feof(fid)
    % Read one line.
    s=fgets(fid);
    % Remove leading blanks.
    s=fliplr(deblank(fliplr(s)));
    % Skip blank or line starting with #.
    if isempty(s) || s(1)=='#'
        continue;
    end
    if hasId
        % Parse id if present.
        [i,~,msg,ni]=sscanf(s,'%d,',1); %#ok<ASGLU>
    s=s(ni:end);
    else
        i=nan;
    end
    if hasName
        % Parse name if present.
        [nn,s]=strtok(s,','); %#ok<STTOK>
        % Remove any blanks.
        nn=fliplr(deblank(fliplr(nn)));
        if ~isempty(s) && s(1)==','
            s=s(2:end);
        end
    else
        nn='';
    end
    % Parse numbers.
    [a,n,msg,ni]=sscanf(s,'%g,'); %#ok<ASGLU>
    id(end+1)=i; %#ok<AGROW>
    name{end+1}=nn; %#ok<AGROW>
    pos(:,end+1)=a(1:3); %#ok<AGROW>
    st=nan(3,1);
    cc=nan(3,3);
    switch n
      case 3
        st(:)=0;
      case 4
        st(:)=a(4); % sigma_xyz
      case 5
        st(:)=a([4,4,5]); % sigma_xy, sigma_z
      case 6
        st=a(4:6);
      case 12
        cc(:)=a(4:end);
        st=sqrt(diag(cc));
        anyCov=true;
      otherwise
        error('Bad number of items on CP line.');
    end
    if all(isnan(cc(:)))
        cc=diag(st.^2);
    end
    std(:,end+1)=st; %#ok<AGROW>
    cov(:,:,end+1)=cc; %#ok<AGROW>
end

fclose(fid);

if ~anyCov
    cov=[];
end

pts=struct('id',id,'name',{name},'pos',pos,'std',std,'cov',cov,'fileName',fName);
