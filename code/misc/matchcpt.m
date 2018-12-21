function [i,j]=matchcpt(s,pts)
%MATCHCPT Loaded control point table with DBAT struct.
%
%   [I,J]=MATCHCPT(S,PTS) matches the control point table PTS with the
%   information in the DBAT structure S. Points can be matched
%   using id and/or label, depending on what information is
%   available in PTS. If points are matched by id, S.OP.rawId(I)
%   corresponds to PTS.id(J). If points are matched by label,
%   S.OP.label(I) corresponds to PTS.name(J).
%
%   Points are matched by id if pts.id are non-NaN. Points are matched
%   by label if pts.name are non-empty. If both conditions are true,
%   points are matched using both criteria. An error is issued if the
%   matchings are inconsistent.
%
%See also: LOADCPT, SETCPT.

% How should we match?
matchById=any(~isnan(pts.id));
matchByLabel=any(~cellfun(@isempty,pts.name));

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

