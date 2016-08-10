if (~exist('proj','var'))
    error('The struct imported from photomodeler should be in the proj variable');
end

ctrlNames=[201,202,211,212,101:103,113,112,111,121:123,302,303,321:323,301,311,312];

useImages=1:2;

% Open communication channel.
disp('Opening a DDE channel to PhotoModeler...');
ch=ddeinit('PhotoModeler','Data');
if (ch==0)
    error('Failed to connect to PhotoModeler');
else
    disp('DDE channel opened.');
end

pmVer=5;

ctrlPt101=[103.55;4.35;0.245]+randn(3,1000)*3;

if (1)
    % Create a new project with [m] as units.
    cmd='NewProject 1';
    [ok,num,str]=ddecmd(ch,cmd);
    if (~ok), ddeterm(ch); error('Failed to create new project.'); end

    % Save file name.
    if (pmVer==5)
        saveFile='C:\pm\ddetest\newProj5.pmr';
    else
        saveFile='C:\pm\ddetest\newProj6.pmr';
    end
    
    % Camera string
    cName='D1400L';
    camVal=[8.9996 8.2446 6.6 1280 1024 3.9322 3.1916 3.288e-3 -2.158e-5 ...
            0 1.829e-4 1.883e-4];
    camStr=[sprintf('%.4f ',camVal(1:3)),...
            sprintf('%d ',camVal(4:5)),...
            sprintf('%.4f ',camVal(6:7)),...
            sprintf('%.8f ',camVal(8:end))];

    disp('Adding camera...');
    [ok,num,str]=ddecmd(ch,['AddCamera ',cName,' ',camStr]);
    if (~ok), ddeterm(ch); error('Failed to add camera.'); end
else
    if (pmVer==5)
        loadFile='C:\pm\ddetest\cam5.pmr';
        saveFile='C:\pm\ddetest\cam5s.pmr';
    else        
        loadFile='C:\pm\ddetest\cam6.pmr';
        saveFile='C:\pm\ddetest\cam6s.pmr';
    end
    
    disp(['Opening project ',loadFile,'...']);
    cmd=['OpenProject ',loadFile];
    [ok,num,str]=ddecmd(ch,cmd);
    if (~ok), ddeterm(ch); error('Failed to open project.'); end

    cName='D1400L';
end

for i=useImages
    disp(sprintf('Adding photo #%d',i));
    cmd=sprintf('AddPhoto %d %s %s',i,cName,proj.cameras(i).imName);
    [ok,num,str]=ddecmd(ch,cmd);
    if (~ok), ddeterm(ch); error('Failed to add photo.'); end
end

for i=useImages
    disp(sprintf('Adding photo station #%d',i));
    cmd=sprintf('SetPhotoStation %d %g %g %g 0 %g %g %g',i,...
                proj.cameras(i).outer(1:3),...
                proj.cameras(i).outer([6,5,4])/180*pi);
    [ok,num,str]=ddecmd(ch,cmd);
    if (~ok), ddeterm(ch); error('Failed to set photo station.'); end
end

for i=useImages
    disp(sprintf('Setting photo processing #%d',i));
    cmd=sprintf('SetPhotoProcessing %d 0 1',i);
    [ok,num,str]=ddecmd(ch,cmd);
    if (~ok), ddeterm(ch); error('Failed to set photo processing.'); end
end

for i=1:size(proj.ctrlPts,1)
    disp(sprintf('Defining ctrl pt #%d',i));
    if (pmVer==5)
        nameStr=sprintf('c%d',ctrlNames(i));
    else
        nameStr=sprintf('imp%d c%d',ctrlNames(i),ctrlNames(i));
    end
    if (all(proj.ctrlPts(i,5:7)==0))
        cmd=sprintf('DefineControlPoint %s %g %g %g',nameStr,...
                    proj.ctrlPts(i,2:4));
    else
        cmd=sprintf('DefineControlPoint %s %g %g %g %g %g %g',nameStr,...
                    proj.ctrlPts(i,2:7));
    end
    [ok,num,str]=ddecmd(ch,cmd);
    if (~ok), ddeterm(ch); error('Failed to define control point.'); end
end

newIds=sparse(0,1);
oldIds=sparse(0,1);

nextPointID=1;

% Mapping between ctrl IDs and names.
ctrlIdsMap=sparse(ctrlNames,1,proj.ctrlPts(:,1));
ctrlNamesMap=sparse(proj.ctrlPts(:,1),1,ctrlNames);

% First mark all control points.
markPts=proj.markPts;
ctrlMarkIx=ismember(markPts(:,2),proj.ctrlPts(:,1));
ctrlMarks=markPts(ctrlMarkIx,:);
markPts(ctrlMarkIx,:)=[];

% Mark control points in ascending name order.
[dummy,i]=sort(ctrlNamesMap(ctrlMarks(:,2)));
ctrlMarks=ctrlMarks(i,:);

markPts=[ctrlMarks;markPts];

for i=1:size(markPts,1)
    imNo=markPts(i,1)+1;
    ptId=markPts(i,2);
    if (ismember(ptId,proj.objPts(:,1)) && ismember(imNo,useImages))
        % This mark corresponds to an object point.

        % Have we seen it before?
        if (ptId<=length(newIds) && newIds(ptId)~=0)
            % Yes.
            newId=newIds(ptId);
        else
            % No, allocate new id and remember mapping.
            if (1)
                cmd='GetNextPointID';
                [ok,num,str]=ddecmd(ch,cmd);
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
        if (any(ptId==proj.ctrlPts(:,1)))
            if (pmVer==5)
                ctrlStr=sprintf('c%d',ctrlNames(ptId==proj.ctrlPts(:,1)));
            else
                ctrlStr=sprintf('imp%d c%d',...
                                ctrlNames(ptId==proj.ctrlPts(:,1)),...
                                ctrlNames(ptId==proj.ctrlPts(:,1)));
            end
        else
            ctrlStr='';
        end

        disp(sprintf('Setting mark for point %d',newId));
        if (isempty(ctrlStr))
            cmd=sprintf('MP %d %d %g %g',imNo,newId,xy);
        else
            cmd=sprintf('MP %d %d %g %g %s',imNo,newId,xy,ctrlStr);
        end
        [ok,num,str]=ddecmd(ch,cmd);
        if (~ok), ddeterm(ch); error('Failed to set mark.'); end
    end
end

nImages=length(proj.cameras);

pts0=ddegetallpts(ch,1:max(newIds));
cams0=ddegetallcams(ch,useImages);

disp('Running bundle...');
cmd='Process 4';
[ok,num,str]=ddecmd(ch,cmd);
if (~ok), ddeterm(ch); error('Failed to run bundle.'); end

pts1=ddegetallpts(ch,1:max(newIds));
cams1=ddegetallcams(ch,useImages);

disp(['Saving to ',saveFile,'...']);
cmd=['SaveProject ',saveFile];
[ok,num,str]=ddecmd(ch,cmd);
if (~ok), ddeterm(ch); error('Failed to save.'); end
disp([saveFile,' saved.']);

ddeterm(ch);

camSize=10;

isCtrl=pts0(:,end)~=0;
plot3(pts0(~isCtrl,2),pts0(~isCtrl,3),pts0(~isCtrl,4),'bx')
hold on
plot3(pts0(isCtrl,2),pts0(isCtrl,3),pts0(isCtrl,4),'b^')
plot3(pts1(~isCtrl,2),pts1(~isCtrl,3),pts1(~isCtrl,4),'rx')
plot3(pts1(isCtrl,2),pts1(isCtrl,3),pts1(isCtrl,4),'r^')
hold off
axis equal
line([pts0(:,2),pts1(:,2),nan(size(pts0(:,2)))]',...
     [pts0(:,3),pts1(:,3),nan(size(pts0(:,2)))]',...
     [pts0(:,4),pts1(:,4),nan(size(pts0(:,2)))]','marker','none',...
     'color',[0,0.5,0]);

for i=useImages
    % Get camera icon.
    [cam,camCol]=cameraicon(camSize,true);
    [m,n,p]=size(cam);
		
    % Camera center.
    CC=cams0(i,2:4)';
		
    % Camera orientation.
    ang=cams0(i,5:7);
		
    RR=pm_eulerrotmat(ang);
		
    % Apply transformation.
    T=RR*[eye(3),-CC];
    T(4,4)=1;
    cam1=reshape(applyhomoxform(inv(T),reshape(cam,m*n,p)')',m,n,p);
		
    hold on
    surf(cam1(:,:,1),cam1(:,:,2),cam1(:,:,3),camCol,... 
         'tag',sprintf('cam%d',i),'facealpha',0.25);
    hold off
end

for i=useImages
    % Get camera icon.
    [cam,camCol]=cameraicon(camSize,true);
    [m,n,p]=size(cam);
		
    % Camera center.
    CC=cams1(i,2:4)';
		
    % Camera orientation.
    ang=cams1(i,5:7);
		
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
