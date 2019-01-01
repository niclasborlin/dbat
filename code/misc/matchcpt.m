function [i1,j1,i2,j2,j3]=matchcpt(s,pts,match)
%MATCHCPT Loaded control point table with DBAT struct.
%
%   [I,J]=MATCHCPT(S,PTS,MATCH) matches the control point table PTS with the
%   control points stored in the DBAT structure S, i.e., among points with
%   s.prior.OP.isCtrl set. Points can be matched using id and/or label,
%   depending on the MATCH parameter. If MATCH=='id', S.OP.rawId(I) corresponds
%   to PTS.id(J). If MATCH=='label', S.OP.label(I) corresponds to PTS.name(J).
%   If MATCH=='both', both will be true. If MATCH=='auto' (default), points are
%   matched by id if pts.id are non-NaN. Points are matched by label if
%   pts.name are non-empty. If both conditions are true, points are matched
%   using both criteria. An error is issued if the matchings are inconsistent.
%
%   [I1,J1,I2,J2,J3]=... returns control point matches in I1,J1 and check point
%   matches in I2,J2, i.e. among points with s.prior.OP.isCheck set. The
%   index vector J3 returns the indices for unmatched points in PTS.
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

% Keep track of matched points.
matchedPts=false(size(pts.id));

% Match among control points.
[i1,j1]=domatch(s,pts,s.prior.OP.isCtrl,matchByLabel,matchById);
matchedPts(j1)=true;

if nargout>2
    [i2,j2]=domatch(s,pts,s.prior.OP.isCheck,matchByLabel,matchById);
    matchedPts(j2)=true;
    j3=find(~matchedPts);
end

% Match among points with sel==true.
function [i,j]=domatch(s,pts,sel,matchByLabel,matchById)

il=[]; jl=[]; ii=[]; ji=[];

% Only match among selected points.
ci=find(sel)';

if matchByLabel
    [~,iii,jl]=intersect(s.OP.label(ci),pts.name);
    il=ci(iii);
end

if matchById
    [~,iii,ji]=intersect(s.OP.rawId(ci),pts.id);
    ii=ci(iii);
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


