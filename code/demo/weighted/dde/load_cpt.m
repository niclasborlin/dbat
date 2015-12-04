function [idsToSimulate,idsToRecord,ctrlId,pt,ptPMerr,PMid,allId,allPMerr]=load_cpt(year,noFWpts,verb)
%[idsToSimulate,idsToRecord,ctrlId,pt,ptPMerr,PMid,allId,allPMerr]=load_cpt(year,noFWpts[,verb])
%idsToSimulate - vector of ids of control points to use as input to simulation.
%idsToRecord   - vector of ids of control points to record as output of
%                simulation. 
%pt            - 3-by-n PTGaussian with mean and uncertainty values
%                loaded from edfdata/orig/cpdata<year>.txt.
%ctrlId        - n-vector of control ids of pt.
%ptPMerr       - 3-by-n PTGaussian with mean zero and variance derived
%                from Photomodeler, corresponding to pt.
%PMid          - m-by-3 sparse array with mapping from control id to
%                Photomodeler id. PMid(I,:)=[ID,CP,FW] is Photomodeler ID
%                for control id I, CP==1 if point is input control point and
%                FW==1 if point is calculated by forward intersection.
%allId         - k-by-3 sparse array with mapping from Photomodeler id to
%                control id. allId(J,:)=[ID,CP,FW] is control ID for
%                Photomodeler id J, CP==1 if point is input control point
%                and FW==1 if point is calculated by forward
%                intersection.
%allPMerr      - 3-by-l PTGaussian with mean zero and variance derived
%                from Photomodeler. allPMerr(:,I) corresponds to control
%                id I.

% $Id$

if nargin<3, verb=false; end

if verb, disp(sprintf('Loading %d cp data',year)); end
% Exported from excel sheet.
Z=load(sprintf('edfdata/orig/cpdata_%d.txt',year));
% Mapping between control id and Photomodeler id.
tbl=load(sprintf('edfdata/orig/%d_cpids.txt',year));
    
% Control points used as control points.
isCtrl=~isnan(tbl(:,3));
% Control points used as object points.
isObj=isnan(tbl(:,3)) & tbl(:,2)>0;
if noFWpts
    isFW=false(size(tbl,1),1);
else
    % Forward intersection points.
    isFW=isnan(tbl(:,3)) & tbl(:,2)<0;
end

% Largest used id.
mxId=max(max(tbl(:,2:3)));
    
% Create mapping btn ctrl pt id and PM id. Second column=1 if point
% is used as control point. Third column=1 if point is a forward
% intersection point only.
ctrlId=sparse([tbl(isCtrl,1);tbl(isObj,1);tbl(isFW,1)],1,[tbl(isCtrl,3);tbl(isObj,2);tbl(isFW,1)+10000]);
ctrlId(tbl(isCtrl,1),2)=1;
ctrlId(1,3)=0;
ctrlId(tbl(isFW,1),3)=1;
allId=ctrlId;

PMid=sparse([tbl(isCtrl,3);tbl(isObj,2);tbl(isFW,1)+10000],1,[tbl(isCtrl,1);tbl(isObj,1);tbl(isFW,1)]);
PMid(tbl(isCtrl,3),2)=1;
PMid(1,3)=0;
PMid(tbl(isFW,1)+10000,3)=1;

cId=1:502;
pId=full(PMid(cId,1))';
cId2=zeros(size(cId));
cId2(pId~=0)=full(ctrlId(pId(pId~=0),1));
inXL=ismember(cId,Z(2:end,1));
isC=full(PMid(cId,2))';

if verb
    disp(sprintf('Excel sheet contains %d ctrl pts with id>0',...
                 nnz(Z(:,1)>0)));
end
if verb
    disp(sprintf('Of those, %d are also in project',...
                 nnz(PMid(Z(2:end,1)))));
end

if verb
    disp(sprintf('Project contains %d ctrl pts used as ctrl pts',...
                      nnz(pId~=0 & isC)));
end
if verb
    disp(sprintf('Of those, %d are also in excel sheet',...
                 nnz(pId~=0 & isC & inXL)));
end
if verb
    disp(sprintf('Project contains %d ctrl pts used as tie pts',...
                 nnz(pId~=0 & ~isC)));
end
if verb
    disp(sprintf('Of those, %d are also in excel sheet',...
                 nnz(pId~=0 & ~isC & inXL)));
end
if verb
    disp(sprintf('Project contains %d pts not in excel sheet',...
                 nnz(pId~=0 & ~inXL)));
end
if verb
    disp(cId(pId~=0 & ~inXL))
end

idsToSimulate=cId(pId~=0 & isC & inXL);
idsToRecord=[cId(pId~=0 & ~isC & inXL),tbl(isFW,1)'+10000];

% Remove points 339, 342, 470 due to excessive Z error in 2007 survey.
idsToSimulate=setdiff(idsToSimulate,[339,342,470]);

if year==1990
    % Clear some parameters due to incomplete data in excel sheet.
    Z(ismember(Z(:,1),369:400),[5:8,11])=0;
    Z(ismember(Z(:,1),470),8)=0;
end

n=size(Z,1);

[dummy,i]=sort(Z(:,1));
Z=Z(i,:);

% Create PTGaussian points and ids for all control points.
ctrlId=Z(:,1)';
pt=PTGaussian(zeros(3,n));
ptPMerr=PTGaussian(zeros(3,n));
for i=1:n
    % xyz coordinates in m for this point.
    xyz=Z(i,[3,4,9])';
    % Planimetric error, [major axis, minor axis], convert to semi-axis
    % and from mm to m
    planError=PTGaussian(zeros(2,1),diag((Z(i,[7,8])/1000/2).^2));
    % Angle between north and major semi-axis (gon), positive direction
    % clockwise.
    theta=Z(i,6);
    % Corresponding angle between positive x-axis and major semi-axis
    % (radians), positive direction anti-clockwise.
    alpha=pi/2-theta/400*2*pi;
    % Create a matrix to rotate planimetric error ellipse.
    R=rotmat2d(PTGaussian(alpha));
    R=R.mean;
    % Rotate planimetric error ellipse.
    rotPlanError=R*planError;
    % Extend to include vertical error.
    xyzError=[rotPlanError;PTGaussian(0,(Z(i,11)/1000/2).^2)];
    % Construct point from mean and errors.
    if (0)
        if (i==1 && verb), disp('Converting to xyz std error (diagonal cov matrix)'); end
        pt(:,i)=PTGaussian(xyz,diag(diag(xyzError.cov)));
    else
        if (i==1 && verb), disp('Using full planimetric cov matrix'); end
        pt(:,i)=PTGaussian(xyz,xyzError.cov);
    end        
    % Determine 3D variance from PM
    if ctrlId(i)>0
        tblRow=find(tbl(:,1)==PMid(ctrlId(i),1));
        if ~isempty(tblRow)
            PMstd=tbl(tblRow,7:9);
            ptPMerr(:,i)=PTGaussian(zeros(3,1),diag(PMstd.^2));
        end
    end
end

% Generate allPMerr
if nargout>7
    % Get variance for non-forward intersection points.

    % Get highest control id.
    m=max(allId(allId(:,3)==0,1));
    % Preallocate.
    PMstd=zeros(3,m);
    PMval=nan(3,m);
    for i=1:m
        tblRow=find(tbl(:,1)==PMid(i,1));
        if ~isempty(tblRow)
            PMstd(:,i)=tbl(tblRow,7:9)';
            PMval(:,i)=0;
        end
    end
    allPMerr=PTGaussian(PMval,diag(PMstd(:).^2));
end
