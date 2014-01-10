function [s,rms]=resect(s0,cams,cpId,chkId)
%RESECT Perform spatial resection on cameras in a project.
%
%   S=RESECT(S0,CAMS,CP_ID) uses the 3-point algorithm to performs spatial
%   resection on the cameras stations of the project dbatstruct S0 listed in
%   the N-vector CAMS using the object points listed in CP_ID as control
%   points. CP_ID must contain at least 3 object points IDs visible in each
%   of the camera stations. If CP_ID is longer than 3, the first 3 IDs
%   visible in each image is used. All object points in S0 not used to
%   compute the resection is used to as check points to distinguish between
%   possible solutions.
%
%   S=RESECT(S0,'all',CP_ID) does spatial resection on all cameras.
%
%   S=RESECT(S0,CAMS,CP_ID,CHK_ID) uses the object points with IDs in CHK_ID
%   as check points. CP_ID and CHK_ID may contain the same IDs.
%
%   [S,RES]=... also returns the rms RES of the residuals of the check
%   points. A failed resection is indicated by a NaN rms.
%
%References:
%Haralick et al., "Review and Analysis of Solutions of the Three Point
%   Perspective Pose Estimation Problem", Int J Comp Vis, 13(3):331-356,
%   1994.
%McGlone et al., "Manual of Photogrammetry", Chapter 11.1.3.4,
%   pp. 786-788, American Society of Photogrammetry and Remote Sensing,
%   2004.
%
%See also: PM_RESECT_3PT.

% $Id$

% Handle defaults.
if nargin<4, chkId=s0.OPid; end

if strcmp(cams,'all'), cams=1:size(s0.EO,2); end

s=s0;

% Remove lens distortion from measured coordinates.
xy=reshape(pm_multilenscorr1(diag([1,-1])*s0.markPts,s0.IO,s0.nK,s0.nP, ...
                             s0.ptCams,size(s0.IO,2)),2,[]);

rms=nan(size(cams));
% For each camera.
for i=1:length(cams)
    camIx=cams(i);

    % Create camera calibration matrix.
    IO=s0.IO(:,s0.cams(camIx));
    K=diag([-IO(3),-IO(3),1]);
    K(1:2,3)=IO(1:2);
    
    % What control points are visible in this camera?
    vis=find(ismember(cpId,s0.OPid(s0.vis(:,camIx))));

    % Pick the first 3 pts if we have more.
    if length(vis)>3, vis=vis(1:3); end
       
    if length(vis)==3
        % We have the 3 pts we need.
        
        useId=cpId(vis);

        % Normalize all measured coordinates visible in this image.
        pt2=xy(:,s0.colPos(s0.vis(:,camIx),camIx));
        pt2N=K\homogeneous(pt2);
    
        % Corresponding object pts.
        pt3=s0.OP(:,s0.vis(:,camIx));

        visId=s0.OPid(s0.vis(:,camIx));

        [P,PP,res]=pm_resect_3pt(pt3,pt2N,ismember(visId,useId),true);
        rms(i)=min(res);
        
        s.EO(1:3,camIx)=euclidean(null(P));
        s.EO(4:6,camIx)=derotmat3d(P(:,1:3));
    end
end
