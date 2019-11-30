function [s,imDir,cptFile,EOfile]=parseinput(input,docFile)
%PARSEINPUT Parse input section of a DBAT XML file.
%
%    S=PARSEINPUT(INPUT,DOCFILE) parses the input XML block INPUT from
%    a DBAT XML script file. The result is returned in the DBAT
%    structure S. The string DOCFILE should contain the path name
%    of the XML file and is used to determine base directories.
%
%    INPUT must contain the fields 'cameras', 'images',
%    'image_pts' and may contain the fields 'Attributes','
%    'ctrl_pts', 'check_pts'.
%
%See also: PARSECAMERAS, PARSEIMAGES, GETATTRBASEDIR.

narginchk(2,2);

inputFields={'Attributes','ctrl_pts','check_pts','images', ...
             'prior_eo','image_pts','cameras','c'};
[ok,msg]=checkxmlfields(input,inputFields,[false,false,false,true, ...
                    false,true,true,false]);
if ~ok, error('DBAT XML script input error: %s',msg); end

% Check if a base directory was specified.
[baseDir,rawBaseDir]=getattrbasedir(input,docFile);

% Warn if base dir does not exist.
if ~isempty(baseDir) && ~exist(baseDir,'dir')
    warning('Base directory %s (%s) does not exist!',baseDir, ...
            rawBaseDir);
end

% Parse cameras.
cams=parsecameras(input.cameras,baseDir);

camModels=cellfun(@(x)x.Model,cams);
if ~all(isnan(camModels)) && ~isscalar(unique(camModels))
    error('Multiple camera models not supported');
end

% Parse images.
ims=parseimages(input.images,baseDir,docFile);

imDir=ims.imDir;

% Parse image points.
allPts=parseimagepts(input.image_pts,baseDir);
% Merge all measured points.
pts=allPts{1};
for i=2:length(allPts)
    pts.id=cat(2,pts.id,allPts{i}.id);
    pts.im=cat(2,pts.im,allPts{i}.im);
    pts.pos=cat(2,pts.pos,allPts{i}.pos);
    pts.std=cat(2,pts.std,allPts{i}.std);
end

% Load control and check points, if any.
ctrlPts=struct('id',zeros(1,0),'fileName','');
if isfield(input,'ctrl_pts')
    ctrlPts=parsectrlpts(input.ctrl_pts,baseDir);
end
cptFile=ctrlPts.fileName;

checkPts=struct('id',zeros(1,0));
if isfield(input,'check_pts')
    checkPts=parsectrlpts(input.check_pts,baseDir);
end

if ~isempty(intersect(ctrlPts.id,checkPts.id))
    error('Point cannot be both control and check points');
end
    
% Number of images.
nImages=length(ims.id);

% Number of mark points.
nMarkPts=length(pts.id);

% Number of object points.
nOP=length(unique([pts.id,ctrlPts.id,checkPts.id]));

s=emptydbatstruct(nImages,nOP,nMarkPts);

s=setdbatcamsandimages(s,cams,ims);

s=setdbatpts(s,ctrlPts,checkPts,pts);

% Parse and load prior EO value, if any.
EOfile='';
if isfield(input,'prior_eo')
    priorEO=parseprioreo(input.prior_eo,baseDir);
    EOfile=priorEO.fileName;
    s.prior.EO.val(1:3,:)=priorEO.pos;
    s.prior.EO.val(4:6,:)=priorEO.ang;
    s.prior.EO.std(1:3,:)=priorEO.std;
    s.prior.EO.std(4:6,:)=priorEO.angStd;
end

s=parseblockvariant(s);
s=buildparamtypes(s);

% s.zz.cams=cams;
% s.zz.ims=ims;
% s.zz.pts=pts;
% s.zz.ctrlPts=ctrlPts;
% s.zz.checkPts=checkPts;
