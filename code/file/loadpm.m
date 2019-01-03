function [prob,err]=loadpm(name,varargin)
%LOADPM Load photomodeler export text file.
%
%   S=LOADPM(NAME) loads the Photomodeler export text file NAME into a
%   struct S with fields listed below. The image size [width, height] in
%   pixels is extracted from the image files, if present, or from the
%   extended header (last 2 elements on row 2 in NAME).
%
%   S=LOADPM(NAME,IMSZ) with IMSZ=[w,h] overrides any information stored
%   in NAME or the image files.
%
%   ...=LOADPM(...,FALSE) also loads the feature information (slower).
%
%   [S,ERR]=LOADPM(...) will return an error string ERR on error instead of
%   throwing an error.
%
%   The struct S has the following fields:
%       job      - struct with fields
%                  fileName  - name of the loaded file, i.e. NAME.
%                  title     - job title.
%                  tol       - required solution tolerance.
%                  maxIter   - maximum number of iterations.
%                  defStd    - vector of default standard deviations 
%                              [photo-coordinates (mm),control points (m),
%                               object points (m),camera position x,y,z (m),
%                               camera angles kappa,phi,omega (degrees)].
%                  defCam    - vector of default calibrated camera
%                              [calibrated focal length (mm), 
%                               principal point xp,yp (mm), 
%                               format size xs,ys (mm),
%                               lens distortion parameters K1,K2,K3,P1,P2].
%                  defCamStd - standard deviation for defCam parameters.
%                  imSz      - [width, height] in pixels, or [NaN, NaN].
%       images   - array of structs with fields
%                  imName   - image name.
%                  outer    - outer orientation parameters 
%                             [x,y,z (m), kappa,phi,omega (degrees)]
%                  outerStd - standard deviation of outer parameters.
%                  outerCov - camera position covariances (xy,xz,yz) (m)
%                  inner    - [focal length, principal point xp,yp, 
%                              format size xs,ys,
%                              lens distortion K1,K2,K3,P1,P2].
%                  innerStd - standard deviation of inner parameters.
%                  imSz     - [width, height] in pixels, or [NaN, NaN].
%                  id       - same as image number.
%                  label    - image name stripped of longest common path.
%       ctrlPts  - array of control points [id,pos x,y,z (m),stdev x,y,z (m)].
%       checkPts - array of check points [id,pos x,y,z (m),stdev x,y,z (m)]. 
%                  Will always be empty on loading, but can be set later.
%       objPts   - array of object points [id,pos x,y,z (m),stdev x,y,z (m)].
%       priorCamPos - array of prior camera positions
%                     [id,pos x,y,z (m),stdev x,y,z (m)].
%       markPts  - array of marked points [photo#,id,pos x,y(px),std x,y(mm)].
%       features - cell array of vectors with points numbers in each
%                  feature (optional).
%       featVis  - array with rows [i,j] indicating that feature j is 
%                  visible in photo i (optional).
%
%   Note: LOADPM does little validation of the file content and e.g. does
%   not verify that a reconstructed object point has measured points in at
%   least two images.


skipFeatures=true;
globalImSz=nan(1,2);

if ~isempty(varargin) && islogical(varargin{end})
    skipFeatures=varargin{end};
    varargin(end)=[];
end

if ~isempty(varargin) && isnumeric(varargin{1}) && length(varargin{1})==2
    globalImSz=reshape(varargin{end},1,2);
    varargin(end)=[];
end

if ~isempty(varargin)
    error('LOADPM: Too mary parameters');
end
    
err=[];
prob=[];

% Open file.
fid=fopen(name,'rt');
if fid<0
	err=['Failed to open file ',name];
	if (nargout<2), error(err); else return; end
end

% Get file size for progress bar.
if fseek(fid,0,'eof')~=0
	err='Failed to get file size';
	if nargout<2, error(err); else return; end
end
sz=ftell(fid);
% Rewind.
fseek(fid,0,'bof');

% Initialize progress bar.
h=waitbar(0,'Loading Photomodeler export file');


% Get project global information from extended header.

% Line 1: Title
title=fgetl(fid);
% Line 2: Tolerence, maxIter, and optionally width, height of image.
s=fgetl(fid); tol=sscanf(s,'%g');
% Line 3: Default point stdev.
s=fgetl(fid); defStd=sscanf(s,'%g');
% Line 4: Default camera parameters [
s=fgetl(fid); defCam=sscanf(s,'%g');
% Default stdev of camera params.
s=fgetl(fid); defCamStd=sscanf(s,'%g');

% Extract the image size if it was present at the end of row 2.
if length(tol)>2 && any(isnan(globalImSz))
    globalImSz=tol(3:4);
end

% Package global information.
job=struct('fileName',name,'title',title,'tol',tol(1),'maxIter',tol(2),...
		   'defStd',defStd,'defCam',defCam,'defCamStd',defCamStd,...
           'imSz',globalImSz);

% Created images struct array.
images=struct('imName',cell(0,1),... % Image file names.
              'outer',zeros(0,6),... % Outer orientation parameter values...
                                 ... % [X,Y,Z (m),kappa,phi,omega (degrees)]
              'outerStd',zeros(0,6),...  % ...standard deviations...
              'outerCov',zeros(0,3),...  % ...and XYZ covariances.
              'inner',zeros(0,10),... % Inner orientation parameter values...
                                ...  % [c,xp,yp,xs,ys,K1,K2,K3,P1,P2]
              'innerStd',zeros(0,10),...  % ...and standard deviations.
              'imSz',zeros(0,2),...    % Image size [w,h].
              'label',cell(0,1),... % image id
              'id',nan); % Image label
% Scan the file line by line.
while ~feof(fid)
    s=fgetl(fid);
    
    if ~ishandle(h)
        % Waitbar closed, abort.
        fclose(fid);
        err='Aborted by user';
        if nargout<2, error(err); else return; end
    elseif rem(length(images)+1,10)==0
        % Update progressbar every 10 input lines.
        waitbar(ftell(fid)/sz,h);
    end
    
    % Here we expect sequence of photo blocks, each on the format
    % (N=zero-based image number)

    % N FILE_NAME
    % N X    Y    Z    KAPPA PHI OMEGA
    % N STDX STDY etc.
    % N C    XP    YP XS YS K1 K2 K3 P1 P2
    % N STDC STDYP etc.

    % Get photo name.
    [photo,count,err,next]=sscanf(s,'%d'); %#ok<ASGLU>
    if isempty(photo)
        % Photo block sequence terminated by blank line.
        break;
    end
    imName=s(next:end);

    % Determine image size.
    imSz=globalImSz;
    if any(isnan(imSz)) && exist(imName,'file')
        info=imfinfo(imName);
        imSz=[info.Width, info.Height];
    end
    
    % Get outer parameters.
    s=fgetl(fid);
    outer=sscanf(s,'%g')';
    outer(1)=[];
    
    % Get outer stdevs.
    s=fgetl(fid);
    outerStd=sscanf(s,'%g')';
    outerStd(1)=[];
    
    % Get outer covariances (is this really used?).
    s=fgetl(fid);
    outerCov=sscanf(s,'%g')';
    if isempty(outerCov)
        outerCov=nan(1,3);
    else
        outerCov(1)=[];
    end
    
    % Get inner parameters.
    s=fgetl(fid);
    inner=sscanf(s,'%g')';
    inner(1)=[];
    
    % Get inner stdevs.
    s=fgetl(fid);
    innerStd=sscanf(s,'%g')';
    innerStd(1)=[];
    
    id=length(images)+1;
    images(end+1)=struct('imName',imName,'outer',outer,...
                         'outerStd',outerStd,'outerCov',outerCov,...
                         'inner',inner,'innerStd',innerStd,...
                         'imSz',imSz,'label',imName,'id',id); %#ok<AGROW>
end
if ishandle(h), waitbar(ftell(fid)/sz,h); end

% Find longest common path of image names.
imPaths={images.imName};
while true
    % Strip last part of path name.
    imDirs=cellfun(@fileparts,imPaths,'uniformoutput',false);
    commonDir=unique(imDirs);
    % Stop when we have a unique common path
    if length(commonDir)<2
        break
    end
    % Otherwise, strip another level.
    imPaths=imDirs;
end

% Use image names sans common path as image labels.
if ~isempty(commonDir) && ~isempty(commonDir{1})
    c=commonDir{1};
    for i=1:length(images)
        images(i).label(1:length(c)+1)='';
    end
end

% Next expected block is control points.

ctrlPts=zeros(0,7);
nCtrlPts=0;
while ~feof(fid)
	s=fgetl(fid);
	if ~ishandle(h)
		% Waitbar closed, abort.
		fclose(fid);
		err='Aborted by user';
		if nargout<2, error(err); else return; end
	elseif rem(nCtrlPts+1,10)==0
        % Update progressbar every 10 control points.
		waitbar(ftell(fid)/sz,h);
	end
	% Get control point data [id,x,y,z,sx,sy,sz].
	cp=sscanf(s,'%g')';
	if isempty(cp)
        % End of control point list.
		break;
	end
    nCtrlPts=nCtrlPts+1;
    if size(ctrlPts,1)<nCtrlPts
        % Expand by 1000 points at a time to avoid memory fragmentation.
        ctrlPts(end+1000,1)=0; %#ok<AGROW>
    end
    ctrlPts(nCtrlPts,:)=cp;
end
if ishandle(h), waitbar(ftell(fid)/sz,h); end
% Trim unused memory.
ctrlPts=ctrlPts(1:nCtrlPts,:);

% Just set check points to empty.
checkPts=zeros(0,7);

% Next expected block is object points.
objPts=zeros(0,7);
nObjPts=0;
while ~feof(fid)
    s=fgetl(fid);
    if ~ishandle(h)
        % Waitbar closed, abort.
        fclose(fid);
        err='Aborted by user';
        if nargout<2, error(err); else return; end
    elseif rem(nObjPts+1,100)==0
        % Update progressbar every 100 object points.
        waitbar(ftell(fid)/sz,h);
    end
    % Get object point data [id,x,y,z,sx,sy,sz].
    op=sscanf(s,'%g')';
    if isempty(op)
        % End of object point list.
        break;
    end
    nObjPts=nObjPts+1;
    if size(objPts,1)<nObjPts
        % Expand by 10000 points at a time to avoid memory fragmentation.
        objPts(end+10000,1)=0; %#ok<AGROW>
    end
    objPts(nObjPts,:)=op;
end
if ishandle(h), waitbar(ftell(fid)/sz,h); end
% Trim unused memory.
objPts=objPts(1:nObjPts,:);

isCtrl=ismember(objPts(:,1),ctrlPts(:,1));

% Next block is mark points.

markPts=zeros(0,6);
nMarkPts=0;
while ~feof(fid)
    s=fgetl(fid);
    if ~ishandle(h)
        % Waitbar closed, abort.
        fclose(fid);
        err='Aborted by user';
        if nargout<2, error(err); else return; end
    elseif rem(nMarkPts+1,1000)==0
        % Update progressbar every 1000 mark points.
        waitbar(ftell(fid)/sz,h);
    end
    % Get marked point info [ph,mp,x,y,sx,sy].
    mp=sscanf(s,'%g')';
    if isempty(mp)
        % End of mark points list.
        break;
    end
    nMarkPts=nMarkPts+1;
    if size(markPts,1)<nMarkPts
        % Expand by 10000 points at a time to avoid memory fragmentation.
        markPts(end+10000,1)=0; %#ok<AGROW>
    end
    markPts(nMarkPts,:)=mp;
end
if ishandle(h), waitbar(ftell(fid)/sz,h); end
markPts=markPts(1:nMarkPts,:);

features=cell(0,0);
if ~skipFeatures
    while ~feof(fid)
        s=fgetl(fid);
        if ~ishandle(h)
            % Waitbar closed, abort.
            fclose(fid);
            err='Aborted by user';
            if nargout<2, error(err); else return; end
        else
            waitbar(ftell(fid)/sz,h);
        end
        % Get pts in each feature.
        ff=sscanf(s,'%d')';
        if isempty(ff)
            break;
        end
        features{ff(1)}=ff(2+(1:ff(2)));
    end
end

featVis=[];
if ~skipFeatures
    while (~feof(fid))
        s=fgetl(fid);
        if ~ishandle(h)
            % Waitbar closed, abort.
            fclose(fid);
            err='Aborted by user';
            if nargout<2, error(err); else return; end
        else
            waitbar(ftell(fid)/sz,h);
        end
        % Get photo-feature pair.
        fp=sscanf(s,'%d')';
        if isempty(fp)
            break;
        end
        featVis=[featVis;fp]; %#ok<AGROW>
    end
end

% Record raw OP ids.
rawOPids=objPts(:,1);

% Create labels from CP ids.
OPlabels=repmat({''},size(rawOPids));
OPlabels(isCtrl)=arrayfun(@int2str,rawOPids(isCtrl),'uniformoutput',false);

% Check if we have both normal and smart points.
isSmartMarkPt=all(markPts(:,5:6)==0,2);
normMarkId=unique(markPts(~isSmartMarkPt,2));
smartMarkId=unique(markPts(isSmartMarkPt,2));

% Are all object point ids increasing?
split=find(diff(objPts(:,1))<0);

if ~isempty(split) && ~isempty(normMarkId) && ~isempty(smartMarkId)
    warning('Found smart points and normal points, renumbering assumed smartpoints');

    maxNormId=max([-inf;normMarkId]);
    minSmartId=min([inf;smartMarkId]);

    % Shift all smart point ids to fall above normal object ids.
    shift=maxNormId+1-minSmartId;
    % Adjust to start on next power of 10.
    %shift=10^ceil(log10(minSmartId+shift)+eps)-minSmartId+1;
    
    % Adjust mark ids
    markPts(isSmartMarkPt,2)=markPts(isSmartMarkPt,2)+shift;
    
    % Adjust obj ids
    isSmartObjPt=ismember(objPts(:,1),smartMarkId);
    isSmartObjPt(1:split)=false;
    objPts(isSmartObjPt,1)=objPts(isSmartObjPt,1)+shift;
end

priorCamPos=zeros(0,7);

prob=struct('job',job,'images',images,'ctrlPts',ctrlPts,'checkPts',checkPts, ...
            'objPts',objPts,'priorCamPos',priorCamPos,'rawOPids',rawOPids, ...
            'OPlabels',{OPlabels}, 'markPts',markPts,'features',{features}, ...
            'featVis',featVis);

fclose(fid);
if ishandle(h), close(h); end
