function [s,rms,fail]=resect(s0,cams,cpId,n,v,chkId)
%RESECT Perform spatial resection on cameras in a project.
%
%   S=RESECT(S0,CAMS,CP_ID) uses the 3-point algorithm to performs spatial
%   resection on the cameras stations of the project dbatstruct S0 listed in
%   the N-vector CAMS using the object points listed in CP_ID as control
%   points. CP_ID must contain at least 3 object points IDs visible in each
%   of the camera stations. If CP_ID is longer than three, the three points
%   covering the largest triangle in each image is used. All object points
%   in S0 not used to compute the resection are used as check points to
%   distinguish between possible solutions.
%
%   S=RESECT(S0,CAMS,CP_ID,N), where N is an integer or Inf, returns the
%   resection based on the N largest triangles that result in the smallest
%   reprojection error.
%
%   S=RESECT(S0,CAMS,CP_ID,N,V), where V is a scalar 0<V<1, tries the N
%   largest triangles that also have an area of at least V times the largest.
%
%   S=RESECT(S0,'all',...) does spatial resection on all cameras.
%
%   S=RESECT(S0,CAMS,CP_ID,N,V,CHK_ID) uses the object points with IDs in
%   CHK_ID as check points. CP_ID and CHK_ID may contain the same IDs.
%
%   [S,RES]=... also returns the rms RES of the residuals of the check
%   points.
%
%   [S,RES,FAIL]=... also returns a flag FAIL which is TRUE if any
%   of the requested resections failed.
%
%   References:
%     Haralick, Lee, Ottenberg, Nölle (1994), "Review and Analysis of Solutions
%       of the Three Point Perspective Pose Estimation Problem".
%       Internation Journal of Computer Vision, 13(3):331-356.
%     McGlone, Mikhail, Bethel, eds. (2004), "Manual of Photogrammetry",
%       5th ed., Chapter 11.1.3.4, pp. 786-788. American Society of
%       Photogrammetry and Remote Sensing.
%
%See also: PM_RESECT_3PT.

% Handle defaults.
if nargin<4, n=1; end
if nargin<5, v=0; end
if nargin<6, chkId=s0.OP.id; end

if strcmp(cams,'all'), cams=1:size(s0.EO.val,2); end

fail=false;

s=s0;

% Remove lens distortion from measured coordinates.
xy=reshape(pm_multilenscorr1(diag([1,-1])*s0.IP.val,s0.IO.val,s0.IO.model, ...
                             s0.IO.sensor,s0.IP.cam,size(s0.IO.val,2)),2,[]);

rms=nan(size(cams));
% For each camera.
for i=1:length(cams)
    camIx=cams(i);

    % Create camera calibration matrix.
    IO=s0.IO.val(:,camIx);
    K=diag([-IO(1),-IO(1),1]);
    K(1:2,3)=IO(2:3);
    
    % What control points are visible in this camera?
    vis=find(ismember(cpId,s0.OP.id(s0.IP.vis(:,camIx))));

    % If we have more than 3 points, pick the ones covering the largest
    % measured area.
    if length(vis)>3
        % Visible measured points.
        meaIx=find(s0.IP.vis(:,camIx) & ismember(s0.OP.id,cpId)');
        mea=xy(:,s0.IP.ix(meaIx,camIx));
        [~,~,T,A]=largesttriangle(mea);
        % Try triangles among N best and above v times highest area.
        tryPt=T((1:end)'<=n & A>=v*A(1),:);
        % Ids of points to try.
        tryId=s0.OP.id(meaIx(tryPt));
    elseif length(vis)==3
        % We have max 3, use all.
        tryId=cpId(vis);
    else
        % We have too few. Use this to fall through loop below.
        tryId=[];
    end

    % Ensure that each index/ID triplet is stored row-wise.
    if isvector(tryId)
        tryId=reshape(tryId,1,[]);
    end
    
    bestRes=inf;
    bestP=[];
    
    for j=1:size(tryId,1)
        % Ids of points to try.
        useId=tryId(j,:);
        
        % Normalize all measured coordinates visible in this image.
        pt2=xy(:,s0.IP.ix(s0.IP.vis(:,camIx),camIx));
        pt2N=K\homogeneous(pt2);
    
        % Corresponding object pts.
        pt3=s0.OP.val(:,s0.IP.vis(:,camIx));
        visId=s0.OP.id(s0.IP.vis(:,camIx));
        
        % Only keep check points.
        keep=ismember(visId,union(cpId,chkId));
        pt2N=pt2N(:,keep);
        pt3=pt3(:,keep);
        visId=visId(keep);
        
        [P,PP,res]=pm_resect_3pt(pt3,pt2N,ismember(visId,useId),true); %#ok<ASGLU>
        if ~isempty(res) && min(res)<bestRes
            bestRes=min(res);
            bestP=P;
        end
    end
       
    rms(i)=bestRes;

    if ~isempty(bestP)
        s.EO.val(1:3,camIx)=euclidean(null(bestP));
        s.EO.val(4:6,camIx)=derotmat3d(bestP(:,1:3));
    else
        % Signal failure.
        fail=true;
        s.EO.val(1:6,camIx)=nan;
    end
end
