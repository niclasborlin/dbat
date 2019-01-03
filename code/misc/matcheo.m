function [i1,j1,i2,j2]=matcheo(s,tbl,match)
%MATCHEO Loaded control point table with DBAT struct.
%
%   [I,J]=MATCHEO(S,TBL,MATCH) matches the image table TBL with the images in
%   the DBAT structure S. Images can be matched using id and/or label,
%   depending on the MATCH parameter. If MATCH=='id', S.EO.id(I) corresponds to
%   TBL.id(J). If MATCH=='label', S.EO.label(I) corresponds to TBL.name(J). If
%   MATCH=='both', both will be true. If MATCH=='auto' (default), images are
%   matched by id if TBL.id are non-NaN. Images are matched by label if
%   TBL.name are non-empty. If both conditions are true, images are matched
%   using both criteria. An error is issued if the matchings are inconsistent.
%
%   [I1,J1,I2,J2]=... returns the matched indices in I1,J1. The index vectors
%   I2 and J2 indicate unmatched images in S and TBL, respectively.
%
%See also: LOADEOTABLE.

if nargin<3
    match='auto';
end

% How should we match?
matchById=false;
matchByLabel=false;

switch match
  case 'id'
    matchById=true;
  case 'label'
    matchByLabel=true;
  case 'both'
    matchById=true;
    matchByLabel=true;
  case 'auto'
    matchById=any(~isnan(tbl.id));
    matchByLabel=any(~cellfun(@isempty,tbl.name));
  otherwise
    error('MATCHEO: Bad match string');
end

% Keep track of matched images.
matchedInTbl=false(size(tbl.id));
matchedInS=false(size(s.EO.id));

% Match among control points.
[i1,j1]=domatch(s,tbl,matchByLabel,matchById);
matchedInS(i1)=true;
matchedInTbl(j1)=true;

if nargout>2
    i2=find(~matchedInS);
    j2=find(~matchedInTbl);
end

% Match among points with sel==true.
function [i,j]=domatch(s,tbl,matchByLabel,matchById)

il=[]; jl=[]; ii=[]; ji=[];

if matchByLabel
    [~,il,jl]=intersect(s.EO.label,tbl.name);
end

if matchById
    [~,ii,ji]=intersect(s.EO.id,tbl.id);
end

switch double(matchByLabel)*2+double(matchById)
  case 3 % Both
    % Give labels priority if we found matched for both.
    if ~isempty(il)
        i=il;
        j=jl;
        if ~isequal(il,ii) || ~isequal(jl,ji)
            error('Inconsistent match');
        end
    else
        % Use id match. Verify that all destination labels were empty.
        i=ii;
        j=ji;
    end
  case 2 % match-by-label only
    i=il;
    j=jl;
  case 1 % match-by-id only
    i=ii;
    j=ji;
  case 0
    error('Neither match-by-id nor match-by-label possible.');
end


