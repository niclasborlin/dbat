function [i,j]=matchcpt(s,pts,match)
%MATCHCPT Loaded control point table with DBAT struct.
%
%   [I,J]=MATCHCPT(S,PTS,MATCH) matches the control point table PTS with the
%   information in the DBAT structure S. Points can be matched using id and/or
%   label, depending on the MATCH parameter. If MATCH=='id', S.OP.rawId(I)
%   corresponds to PTS.id(J). If MATCH=='label', S.OP.label(I) corresponds to
%   PTS.name(J). If MATCH=='both', both will be true. If MATCH=='auto'
%   (default), points are matched by id if pts.id are non-NaN. Points are
%   matched by label if pts.name are non-empty. If both conditions are true,
%   points are matched using both criteria. An error is issued if the matchings
%   are inconsistent.
%
%See also: LOADCPT, SETCPT.

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
    matchById=any(~isnan(pts.id));
    matchByLabel=any(~cellfun(@isempty,pts.name));
  otherwise
    error('MATCHCPT: Bad match string');
end

il=[]; jl=[]; ii=[]; ji=[];

if matchByLabel
    [~,il,jl]=intersect(s.OP.label,pts.name);
end

if matchById
    [~,ii,ji]=intersect(s.OP.rawId,pts.id);
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
        if any(~cellfun(@isempty,s.OP.label(i)))
            error('Inconsistent match');
        end
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

