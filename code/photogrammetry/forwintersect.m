function [s,id,res]=forwintersect(s0,ids,skipFixed)
%FORWINTERSECT Perform forward intersection on points in a project.
%
%   S=FORWINTERSECT(S0,ID) computes the OP coordinates of all points in
%   the project S0 with IDs listed in ID using forward intersection.
%   Points visible in too few images are given NaN coordinates.
%
%   S=FORWINTERSECT(S0,'all') computes all OP coordinates.
%
%   S=FORWINTERSECT(S0,...,TRUE) does not compute OP coordinates for points
%   with any fixed coordinate(s), i.e. with a FALSE in the corresponding
%   S0.estOP entry.

%   [S,ID,RES]=... returns the IDs of each computed points in ID and the
%   corresponding rms object space residual in RES.
%
%See also: PM_MULTIFORWINTERSECT.

% $Id$

if strcmp(ids,'all'), ids=s0.OPid; end
   
if nargin<3, skipFixed=false; end

% Remove lens distortion from measured coordinates and convert to mm.
xy=reshape(pm_multilenscorr1(diag([1,-1])*s0.markPts,s0.IO,s0.nK,s0.nP, ...
                             s0.ptCams,size(s0.IO,2)),2,[]);

% Extract wanted points.
if skipFixed
    isFree=all(s0.estOP,1)';
else
    isFree=true;
end
    
i=find(ismember(s0.OPid,ids) & isFree);

% Compute the forward intersection.
[OP,res]=pm_multiforwintersect(s0.IO,s0.EO,s0.cams,s0.colPos,xy,i);

% Store result.
id=s0.OPid(i);
s=s0;
s.OP(:,i)=OP;