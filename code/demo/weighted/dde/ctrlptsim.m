function [X,C,Xi,Ci]=ctrlptsim(varargin)
%CTRLPTSIM Simulate control point pertubation on a PhotoModeler project.
%
%   X=CTRLPTSIM(PROJ,PMVER,K,CPT) uses the PhotoModeler DDE interface to
%   perform K bundle runs with different control point values.  PROJ is a
%   struct returned by LOADPM.  PMVER is the PhotoModeler version number (5
%   or 6).  CPT is an N-by-3 PTGaussian with the control points to simulate,
%   in the order given by PROJ.ctrlPts.  X is an M-by-3 PTGaussian with
%   the collected statistical properties of the result, in the same order
%   as PROJ.objPts.
%
%   All info for the processing is taken from PROJ.  All images are assumed
%   to use the same camera (the default camera).  All image marks without a
%   corresponding object point are silently ignored.
%
%   PhotoModeler should be started with no project loaded.
%
%   Use K=0 to do one bundle run without modifying the control points.  This
%   can be used test that the PhotoModeler bundle function can handle the
%   project.
%
%   [X,C]=... also returns the camera center positions as a PTGaussian.
%   [X,C,Xi,Ci]=... also returns the individual samples.
%
%   ...=CTRLPTSIM(PROJ,PMVER,K,CPT,CAM,CAMIX) uses the camera information
%   with columns CAM=[f sw sh iw ih px py K1 K2 K3 P1 P2]' instead of the
%   info in PROJ.  CAMIX(i) gives the camera number for image i.  Use an
%   explicit CAM to give more precise values than what is exported to the
%   PhotoModeler text file. The elements of CAM are:
%     f       - focal length in mm.
%     sw, sh  - sensor width, height in mm.
%     iw, ih  - image width, height in pixels.
%     px, py  - principal point in mm.
%     K1, ... - lens distortion parameters.
%
%   ...=CTRLPTSIM(PROJ,PMVER,K,CPT,CAM,VERB) gives control of the
%   printout during processing.
%     VERB==0     - no printout.
%     VERB==1     - some printout (default).
%     VERB==2     - lot's of printout.
%
%   ...=CTRLPTSIM(...,TRUE) plots the result before and after the bundle.
%   Only makes sense for small K.

% $Id$

% Default values.
cam=[];
verb=1;
doPlot=false;
Xi=[];
Ci=[];
C=[];
X=[];

% Handle trailing doPlot argument.
if (~isempty(varargin) && islogical(varargin{end}))
    doPlot=varargin{end};
    varargin(end)=[];
end

switch (length(varargin))
case 4
    [proj,pmVer,k,cpt]=deal(varargin{:});
    if (k==0), verb=2; end
case 6
    [proj,pmVer,k,cpt,cam,camIx]=deal(varargin{:});
    if (k==0), verb=2; end
case 7
    [proj,pmVer,k,cpt,cam,camIx,verb]=deal(varargin{:});
otherwise
    error('Illegal number of arguments.');
end

% Handle camera defaults.
if (size(cam,1)~=12), error('CAM must contain 12 elements per column'); end

if (~ismember(pmVer,5:6)), error ('Illegal PhotoModeler version.'); end

[m,n]=size(cpt);
if (~isa(cpt,'PTGaussian') || (m~=3) || (n~=size(proj.ctrlPts,1)))
    if (~isempty(cpt))
        error('Illegal CPT size.');
    end
end

if (k==0)
    % Use the existing control points.
    CPTi=proj.ctrlPts(:,2:7)';
else
    CPTi=sample(cpt,k,3);
end

nImages=length(proj.cameras);

% Ensure that all image files exist.
for i=1:nImages
    if (~exist(proj.cameras(i).imName,'file'))
        error('Image file #%d (%s) does not exist.',i,proj.cameras(i).imName);
    end
end

% Open communication channel.
if (verb>0), disp('Opening a DDE channel to PhotoModeler...'); end
ch=ddeinit('PhotoModeler','Data');
if (ch==0)
    error('Failed to connect to PhotoModeler');
end

for kk=1:max(k,1),
    if (verb>0), fprintf('Loop #%d of %d.\n',kk,max(k,1)); end
    % Create a new project with [m] as units.
    cmd='NewProject 1';
    [ok,num,str]=ddecmd(ch,cmd); %#ok<NASGU>
    if (~ok), ddeterm(ch); error('Failed to create new project.'); end

    if (0)
        % Project file name.
        saveFile=[tempdir,'dummyproject.pmr'];
    end
    
    % Add cameras
    cName=cell(1,size(cam,2));
    for i=1:size(cam,2)
        cName{i}=sprintf('DummyCam%d',i);
        camStr=sprintf('%g ',cam(:,i).mean);

        if (verb>1), disp(sprintf('Adding camera %d...',i)); end
        [ok,num,str]=ddecmd(ch,['AddCamera ',cName{i},' ',camStr]); %#ok<NASGU>
        if (~ok), ddeterm(ch); error('Failed to add camera.'); end
    end
    
    if (verb>0), fprintf('Adding photos: '); end
    for i=1:nImages
        if (verb>1), fprintf(' %d',i); end
        cmd=sprintf('AddPhoto %d %s %s',i,cName{camIx(i)},...
                    proj.cameras(i).imName);
        [ok,num,str]=ddecmd(ch,cmd); %#ok<NASGU>
        if (~ok), ddeterm(ch); error('Failed to add photo.'); end
    end
    if (verb>0), fprintf('.\n'); end

    if (verb>0), fprintf('Adding photo stations: '); end
    for i=1:nImages
        if (verb>1), fprintf(' %d',i); end
        cmd=sprintf('SetPhotoStation %d %g %g %g 0 %g %g %g',i,...
                    proj.cameras(i).outer(1:3),...
                    proj.cameras(i).outer([6,5,4])/180*pi);
        [ok,num,str]=ddecmd(ch,cmd); %#ok<NASGU>
        if (~ok), ddeterm(ch); error('Failed to set photo station.'); end
    end
    if (verb>0), fprintf('.\n'); end

    if (verb>0), fprintf('Setting photo processing: '); end
    for i=1:nImages
        if (verb>1), fprintf(' %d',i); end
        cmd=sprintf('SetPhotoProcessing %d 0 1',i);
        [ok,num,str]=ddecmd(ch,cmd); %#ok<NASGU>
        if (~ok), ddeterm(ch); error('Failed to set photo processing.'); end
    end
    if (verb>0), fprintf('.\n'); end
    
    if (verb>0), fprintf('Defining control points: '); end
    for i=1:size(proj.ctrlPts,1)
        if (verb>1), fprintf(' %d',proj.ctrlPts(i,1)); end
        if (pmVer==5)
            nameStr=sprintf('c%d',proj.ctrlPts(i,1));
        else
            nameStr=sprintf('imp%d c%d',proj.ctrlPts(i,1),proj.ctrlPts(i,1));
        end
        % Set control points.
        if (size(CPTi,1)==3)
            cmd=sprintf('DefineControlPoint %s %g %g %g',nameStr,...
                        CPTi(:,i,kk));
        else
            cmd=sprintf('DefineControlPoint %s %g %g %g %g %g %g',nameStr,...
                        CPTi(:,i,kk));
        end
        [ok,num,str]=ddecmd(ch,cmd); %#ok<NASGU>
        if (~ok), ddeterm(ch); error('Failed to define control point.'); end
    end
    if (verb>0), fprintf('.\n'); end
    
    if (verb>0), fprintf('Setting image marks: '); end
    
    newIds=sparse(0,1);
    oldIds=sparse(0,1);
    
    nextPointID=1;
    
    % Reorder marks to mark all control points first.
    markPts=proj.markPts;
    ctrlMarkIx=ismember(markPts(:,2),proj.ctrlPts(:,1));
    markPts=markPts([find(ctrlMarkIx);find(~ctrlMarkIx)],:);
    
    for i=1:size(markPts,1)
        imNo=markPts(i,1)+1;
        ptId=markPts(i,2);
        if (ismember(ptId,proj.objPts(:,1)))
            % This mark corresponds to an object point.
            
            % Have we seen it before?
            if (ptId<=length(newIds) && newIds(ptId)~=0)
                % Yes.
                newId=newIds(ptId);
            else
                % No, allocate new id and remember mapping.
                if (1)
                    cmd='GetNextPointID';
                    [ok,num,str]=ddecmd(ch,cmd); %#ok<NASGU>
                    if (~ok), ddeterm(ch); error('Failed to get next point ID.'); end
                    newId=num{2};
                else
                    newId=nextPointID;
                    nextPointID=nextPointID+1;
                end
                newIds(ptId)=newId;
                oldIds(newId)=ptId;
            end
            
            % Data for this mark.
            xy=markPts(i,3:4);
            % Is it a control point?
            isCtrl=any(ptId==proj.ctrlPts(:,1));
            if (isCtrl)
                if (pmVer==5)
                    ctrlStr=sprintf('c%d',ptId);
                else
                    ctrlStr=sprintf('imp%d c%d',ptId,ptId);
                end
            else
                ctrlStr='';
            end
            
            if (verb>1), fprintf(' %d (%d)',newId,ptId); end
            if (verb>1 && rem(i,10)==0), fprintf('\n'); end
            if (isempty(ctrlStr))
                cmd=sprintf('MP %d %d %g %g',imNo,newId,xy);
            else
                cmd=sprintf('MP %d %d %g %g %s',imNo,newId,xy,ctrlStr);
            end
            [ok,num,str]=ddecmd(ch,cmd); %#ok<NASGU>
            if (~ok), ddeterm(ch); error('Failed to set mark.'); end
        end
    end
    if (verb>0), fprintf('.\n'); end

    if (verb>0), disp('Running bundle...'); end
    cmd='Process 4';
    [ok,num,str]=ddecmd(ch,cmd); %#ok<NASGU>
    if (~ok), ddeterm(ch); error('Failed to run bundle.'); end
    if (verb>0)
        fprintf('Bundle iters: %d, initial error: %g, final error: %g.\n',...
            num{2:end});
    end
    
    pts=ddegetallpts(ch,1:max(newIds));
    cams=ddegetallcams(ch,1:nImages);

    if (isempty(Xi))
        Xi=zeros(3,size(pts,1),max(k,1));
        Xi(:,:,1)=pts(:,2:4)';
    else
        if (size(pts,1)~=size(Xi,2))
            ddeterm(ch);
            error('The number of object points must remain the same.');
        end
        Xi(:,:,kk)=pts(:,2:4)'; %#ok<AGROW>
    end
    
    if (isempty(Ci))
        Ci=zeros(6,size(cams,1),max(k,1));
        Ci(:,:,1)=cams(:,2:7)';
    else
        if (size(cams,1)~=size(Ci,2))
            ddeterm(ch);
            error('The number of camera stations must remain the same.');
        end
        Ci(:,:,kk)=cams(:,2:7)'; %#ok<AGROW>
    end
    
    if (0)
        disp(['Saving to ',saveFile,'...']);
        cmd=['SaveProject ',saveFile];
        [ok,num,str]=ddecmd(ch,cmd); %#ok<NASGU>
        if (~ok), ddeterm(ch); error('Failed to save.'); end
        disp([saveFile,' saved.']);
    end
end

ddeterm(ch);

% Object point statistics.
Xm=mean(Xi,3);
XC=cov(reshape(Xi,[],size(Xi,3))');
X=PTGaussian(Xm,XC);

Cm=mean(Ci,3);
CC=cov(reshape(Ci,[],size(Ci,3))');
C=PTGaussian(Cm,CC);

camSize=10;

if (doPlot)
    isCtrl=pts(:,end)~=0;
    plot3(pts(~isCtrl,2),pts(~isCtrl,3),pts(~isCtrl,4),'rx')
    hold on
    plot3(pts(isCtrl,2),pts(isCtrl,3),pts(isCtrl,4),'r^')
    hold off
    axis equal

    for i=1:nImages
        % Get camera icon.
        [cam,camCol]=cameraicon(camSize,true);
        [m,n,p]=size(cam);
		
        % Camera center.
        CC=cams(i,2:4)';
		
        % Camera orientation.
        ang=cams(i,5:7);
		
        RR=pm_eulerrotmat(ang);
		
        % Apply transformation.
        T=RR*[eye(3),-CC];
        T(4,4)=1;
        cam1=reshape(applyhomoxform(inv(T),reshape(cam,m*n,p)')',m,n,p);
		
        hold on
        surf(cam1(:,:,1),cam1(:,:,2),cam1(:,:,3),camCol,... 
             'tag',sprintf('cam%d',i));
        hold off
    end
    colormap([1,0,0;0,1,0;0.5,0.5,1]);
    axis equal
end
