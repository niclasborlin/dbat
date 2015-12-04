ddepmterm(ch)

addpath(fullfile(pwd,'dde'),'-end');
ch=ddepminit

cmd='NewProject 1';
[ok,num,str]=ddecmd(ch,cmd);
if (~ok), ddeterm(ch); error('Failed to create new project.'); end

% Camera string
cName='C4040Z';
camVal=[IO(3) IO(end-5:end-4)' IO(end-3:end-2)' abs(IO(1:2))' zeros(1,5)];
camStr=[sprintf('%.4f ',camVal(1:3)),...
        sprintf('%d ',camVal(4:5)),...
        sprintf('%.4f ',camVal(6:7)),...
        sprintf('%.8f ',camVal(8:end))];

disp('Adding camera...');
[ok,num,str]=ddecmd(ch,['AddCamera ',cName,' ',camStr]);
if (~ok), ddeterm(ch); error('Failed to add camera.'); end

ctrlId=OPid(OPctrl);
nonCtrlId=setdiff(OPid,ctrlId);
ctrlStd=[];

useImages=images;

fprintf('Defining control points...');
for i=1:length(ctrlId)
    fprintf(' %d',ctrlId(i));
    nameStr=sprintf('imported c%d',ctrlId(i));
    % Set control points.
    if isempty(ctrlStd)
        cmd=sprintf('DefineControlPoint %s %.10g %.10g %.10g',...
                    nameStr,OP(:,OPid==ctrlId(i)));
    else
        cmd=sprintf('DefineControlPoint %s %.10g %.10g %.10g %.10g %.10g %.10g',...
                    nameStr,simCtrlPts(:,find(ctrlId==idsToSimulate(i))));
    end
    [ok,num,str]=ddecmd(ch,cmd); %#ok<NASGU>
    if (~ok), ddepmterm(ch); error('Failed to define control point.'); end
end
fprintf(' done.\n');

fprintf('Adding photos...');
for i=useImages
    fprintf(' %d',i);
    cmd=sprintf('AddPhoto %d %s',i,cName);
    [ok,num,str]=ddecmd(ch,cmd);
    if (~ok), ddeterm(ch); error('Failed to add photo.'); end
end
fprintf(' done.\n');

fprintf('Adding mark pts for control points: ');

for ii=1:length(ctrlId)
    cmd='GetNextPointID';
    [ok,num,str]=ddecmd(ch,cmd); %#ok<NASGU>
    if (~ok), ddepmterm(ch); error('Failed to get next point ID.'); end
    markId=num{2};

    for jj=find(vis(OPid==ctrlId(ii),:))
        % Data for this mark.
        xy=mVis{jj}(:,OPid(vis(:,jj))==ctrlId(ii));
        imNo=images(jj);
        ctrlStr=sprintf('imported c%d',ctrlId(ii));
        
        fprintf(' %d %d (%d)',imNo,markId,ctrlId(ii));
        cmd=sprintf('MP %d %d %.10g %.10g %s',imNo,markId,xy,ctrlStr);
        [ok,num,str]=ddecmd(ch,cmd); %#ok<NASGU>
        if (~ok), ddepmterm(ch); cmd, error('Failed to set mark.'); end
    end
end
fprintf(' done.\n');

fprintf('Adding mark pts for non-control points: ');

for ii=1:length(nonCtrlId)
    cmd='GetNextPointID';
    [ok,num,str]=ddecmd(ch,cmd); %#ok<NASGU>
    if (~ok), ddepmterm(ch); error('Failed to get next point ID.'); end
    markId=num{2};

    for jj=find(vis(OPid==nonCtrlId(ii),:))
        % Data for this mark.
        xy=mVis{jj}(:,OPid(vis(:,jj))==nonCtrlId(ii));
        imNo=images(jj);
        
        fprintf(' %d %d',imNo,markId);
        cmd=sprintf('MP %d %d %.10g %.10g',imNo,markId,xy);
        [ok,num,str]=ddecmd(ch,cmd); %#ok<NASGU>
        if (~ok), ddepmterm(ch); cmd, error('Failed to set mark.'); end
    end
end
fprintf(' done.\n');

if 1
    fprintf('Adding photo stations...');
    for i=useImages
        fprintf(' %d',i);
        cmd=sprintf('SetPhotoStation %d %g %g %g 0 %g %g %g',i,...
                    EO(1:3,i),EO(4:6,i));
        [ok,num,str]=ddecmd(ch,cmd);
        if (~ok), ddeterm(ch); error('Failed to set photo station.'); end
    end
    fprintf(' done.\n');
end

ddepmterm(ch)
