function [s,id,res]=forwintersect(s0,ids,skipPrior)
%FORWINTERSECT Perform forward intersection on points in a project.
%
%   S=FORWINTERSECT(S0,ID) computes the OP coordinates of all points in
%   the project S0 with IDs listed in ID using forward intersection.
%   Points visible in too few images are given NaN coordinates.
%
%   S=FORWINTERSECT(S0,'all') computes all OP coordinates.
%
%   S=FORWINTERSECT(S0,...,TRUE) does not compute OP coordinates for points
%   with prior observations, i.e. points with either fixed coordinate(s)
%   or with prior observations.
%
%   [S,ID,RES]=... returns the IDs of each computed points in ID and the
%   corresponding rms object space residual in RES.
%
%See also: PM_MULTIFORWINTERSECT.

if any(~isfinite(s0.EO.val(:))), error('Bad or uninitialized EO data'); end
if any(~isfinite(s0.IO.val(:))), error('Bad or uninitialized IO data'); end
    
if strcmp(ids,'all'), ids=s0.OP.id; end
   
if nargin<3, skipPrior=false; end

% Remove lens distortion from measured coordinates and convert to mm.
xy=reshape(pm_multilenscorr1(diag([1,-1])*s0.IP.val,s0.IO.val,s0.IO.model,...
                             s0.IO.sensor,s0.IP.cam,size(s0.IO.val,2)),2,[]);

% Extract wanted points.
if skipPrior
    doEst=all(s0.bundle.est.OP,1)' & ~any(s0.prior.OP.use,1)';
else
    doEst=true;
end
    
i=find(ismember(s0.OP.id,ids)' & doEst);

% Compute the forward intersection.
[OP,res]=pm_multiforwintersect(s0.IO.val,s0.EO.val,1:size(s0.IO.val,2),...
                               s0.IP.ix,xy,i);

% Store result.
id=s0.OP.id(i);
s=s0;
s.OP.val(:,i)=OP;
