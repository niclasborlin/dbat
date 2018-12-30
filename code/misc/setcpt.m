function s=setcpt(s,pts,i,j,isCtrl)
%SETCPT Set control points in DBAT struct.
%
%   S=SETCPT(S,PTS,I,J) copies the control point information from
%   CPT to the DBAT struct S. [I,J] should have been the result of
%   a previous MATCHCPT(S,PTS).
%
%   S=SETCPT(S,PTS,I,J,FALSE) sets the points as check points instead.
%
%   SETCPT copies the prior pos, std, and/or cov to S. Furthermore,
%   the isCtrl and isCheck properties are set. The prior.use and
%   bundle.est fields are set depending on whether the control
%   points are fixed or not. The pts.fileName is recorded in s.proj.
%
%   Labels in S are overwritten by non-blank labels in PTS.
%
%See also: LOADCPT, MATCHCPT.

if nargin<5, isCtrl=true; end

s.proj.cptFile=pts.fileName;

% Copy positions.
s.prior.OP.val(:,i)=pts.pos(:,j);
s.OP.val(:,i)=pts.pos(:,j);

if isempty(pts.cov)
    s.prior.OP.std(:,i)=pts.std(:,j);
else
    error('Individual ctrl pt covariances not yet supported.');
end

% Copy non-blank labels.
for k=1:length(i)
    if ~isempty(pts.name{j(k)})
        s.OP.label{i(k)}=pts.name{j(k)};
    end
end

s.prior.OP.isCtrl(i)=isCtrl;
s.prior.OP.isCheck(i)=~isCtrl;

if isCtrl
    % Fixed control points should neither be used as observations
    % nor be estimated.
    isFixed=all(pts.std(:,j)==0,1);
    s.prior.OP.use(:,i)=repmat(~isFixed,3,1);
    s.bundle.est.OP(:,i)=repmat(~isFixed,3,1);
else
    % Check points should not be used by the bundle but should be estimated.
    s.prior.OP.use(:,i)=false;
    s.bundle.est.OP(:,i)=true;
end
