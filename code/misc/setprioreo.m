function s=setprioreo(s,tbl,i,j)
%SETPRIOREO Set prior EO positions in DBAT struct.
%
%   S=SETPRIOREO(S,TBL,I,J) copies the EO information from TBL to the
%   DBAT struct S. [I,J] should have been the result of a previous
%   MATCHEO(S,TBL).
%
%   SETPRIOREO copies the prior pos, std, and/or cov to S. The
%   prior.use and bundle.est fields are set depending on whether the
%   EO positions are fixed or not. The tbl.fileName is recorded in
%   s.proj.
%
%   Labels in S are overwritten by non-blank labels in TBL.
%
%See also: LOADEOTBL, MATCHEO.

s.proj.EOfile=tbl.fileName;

% Copy positions.
s.prior.EO.val(1:3,i)=tbl.pos(:,j);
s.EO.val(1:3,i)=tbl.pos(:,j);

if isempty(tbl.cov)
    s.prior.EO.std(1:3,i)=tbl.std(:,j);
else
    error('Individual EO pos covariances not yet supported.');
end

% Copy non-blank labels.
for k=1:length(i)
    if ~isempty(tbl.name{j(k)})
        s.EO.label{i(k)}=tbl.name{j(k)};
    end
end

% Fixed EO positions should neither be used as observations nor be
% estimated.
isFixed=tbl.std(:,j)==0;
s.prior.EO.use(1:3,i)=~isFixed;
s.bundle.est.EO(1:3,i)=~isFixed;
