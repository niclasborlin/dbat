function s=setdbatpts(s,ctrlPts,checkPts,imagePts)
%SETDBATPTS Set control, check, and image points in a DBAT structure
%
%   S=SETDBATPTS(S,CTRLPTS,CHECKPTS,IMAGEPTS) sets the control
%   points, check points, and image points in the DBAT structure S.
%
%   An error is thrown if the number of points do not match to
%   insert does not match the number of allocated points.

% Collect all ids
OPid=unique([ctrlPts.id,checkPts.id,imagePts.id]);

% Prior values
priorVal=nan(3,length(OPid));
priorStd=nan(3,length(OPid));
label=repmat({''},1,length(OPid));

% Insert prior control points values and standard deviations.
[~,ia,ib]=intersect(OPid,ctrlPts.id);
if ~isempty(ia)
    priorVal(:,ia)=ctrlPts.pos(:,ib);
    priorStd(:,ia)=ctrlPts.std(:,ib);
    label(ia)=ctrlPts.name(ib);
end

% Ditto for check points
[~,ia,ib]=intersect(OPid,checkPts.id);
if ~isempty(ia)
    priorVal(:,ia)=checkPts.pos(:,ib);
    priorStd(:,ia)=checkPts.std(:,ib);
    label(ia)=checkPts.name(ib);
end

% Store prior values and std
s.prior.OP.val=priorVal;
s.prior.OP.std=priorStd;
s.prior.OP.isCtrl=ismember(OPid,ctrlPts.id);
s.prior.OP.isCheck=ismember(OPid,checkPts.id);
% Use prior value as observation if ctrl pt has non-zero std.
% Never use check points as priors.
s.prior.OP.use=~isnan(s.prior.OP.std) & s.prior.OP.std~=0 & ...
    ~repmat(s.prior.OP.isCheck,3,1);

% Store ids, etc.
s.OP.id=OPid;
s.OP.rawId=OPid;
s.OP.label=label;

% Image points should be sorted by image number, then point id.
[~,i]=sortrows([imagePts.im;imagePts.id]');
s.IP.id=imagePts.id(i);
s.IP.cam=imagePts.im(i);
s.IP.val=imagePts.pos(:,i);
s.IP.std=imagePts.std(:,i);
s.IP.sigmas=unique(s.IP.std(:));

nOP=length(s.OP.id);
nImages=length(s.EO.id);
nMarkPts=length(s.IP.id);

% Create id-to-row mapping
row=sparse(s.OP.id,ones(1,nOP),1:nOP);

% Create visibility matrix.
s.IP.vis=sparse(full(row(s.IP.id)),s.IP.cam,true,nOP,nImages);

% Create indices into the measurement array.
s.IP.ix=sparse(full(row(s.IP.id)),s.IP.cam,1:nMarkPts,nOP,nImages);
