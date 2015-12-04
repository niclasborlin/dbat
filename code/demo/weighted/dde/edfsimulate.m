function edfsimulate(year,pmVer,loops)
%EDFSIMULATE Simulate error of CP uncertainty on PM bundle adjustment.
%
%   EDFSIMULATE(YEAR, PM, N) runs a maximum of N iterations of the
%   simulation for the specified YEAR (1990 or 2007) using Photomodeler
%   version PM (5 or 6).  If any checkpoint file(s) is detected, the
%   simulation is loaded and continued after the last stored
%   iteration.
%
%   The correct Photomodeler version must have been started cleanly before
%   calling this function.
%
%   Intermediate results are stored in checkpoint files
%   tmpdata_YYYY_nofw_odd.mat and tmpdata_YYYY_nofw_even.mat in the current
%   folder, where YYYY is the year in question. The final result is
%   stored in the file done_YYYY_nofw.mat in the current folder. The
%   final result file acts as a third checkpoint file.
%
%   The simulation may be interrupted at any time, including by a power
%   outage. After restarting Windows, the simulation may be restarted and
%   will continue after the last iteration stored in a checkpoint file.

% $Id$

% Ignore points calculated by forward intersection.
noFWpts=true;

% Do we want to run with orientation of all images?
orientAll=true;

% Do we want to clear PM by pre-loading a blank project.
preLoad=true;

% Minutes between saving of data. 0=save every iteration. -1=no saving.
if noFWpts
    saveInterval=0.25;
else
    saveInterval=2;
end

% Verbosity level
verb=1;

% Set up random number generator to generate the same random number every
% time we run this file. Modify the SEED value below to get other random
% sequences.
reset(RandStream.getDefaultStream);

if noFWpts
    pmrFile=fullfile(pwd,'edfdata','orig',[int2str(year),'_nocpt_nofw.pmr']);
else
    pmrFile=fullfile(pwd,'edfdata','orig',[int2str(year),'_nocpt.pmr']);
end
pmExportFile=fullfile(pwd,'edfdata','orig',[int2str(year),'_pmr_export.txt']);

[idsToSimulate,idsToRecord,ctrlId,pt,dummy,PMid]=load_cpt(year,noFWpts);

if year==2007
    [dummy,idsToRecord]=load_cpt(1990,noFWpts);
    % Remove points not appearing in 2007 project.
    idsToRemove=[247,489,492,493,494];
    idsToRecord=setdiff(idsToRecord,idsToRemove);
end

if (~exist(pmrFile,'file'))
    pmrFile
    error('File does not exist');
end

if (~exist(pmExportFile,'file'))
    pmExportFile
    error('File does not exist');
end

if (~exist('proj','var'))
    disp('Loading PM export file.');
    proj=loadpm(pmExportFile);
    drawnow
else
    disp('Using loaded PM export file.');
end

proj.ctrlPts(:,5:7)=0;

if (~exist(pmrFile,'file'))
    pmrFile
    error('File does not exist');
end

ch=ddepminit;

[ptsToSimulate,ptsToRecord,preLoop,startDate,seed,lastSaveWasOdd,files]=...
    loadsimdata(year,noFWpts,'a',idsToSimulate,idsToRecord);

[oddFile,evenFile,doneFile]=deal(files{:});

if (size(ptsToRecord,2)<loops)
    % Expand matrix to avoid memory fragmentation.
    ptsToRecord(1,loops)=0;
end

if (size(ptsToSimulate,2)<loops)
    % Expand matrix to avoid memory fragmentation.
    ptsToSimulate(1,loops)=0;
end

% Remove n values from the random stream to get reproducible but different
% random number sequences.
dummy=randn(1,seed);

start=now;
lastSaveDate=start;

for loop=1:loops
    % Pick one sample from the X distribution.
    simCtrlPts=sample(pt,1);

    if (loop<=preLoop)
        % This simulation has already been performed.
        start=now;
        continue;
    end
    
    %ddepmterm(ch);
    %ch=ddepminit;
    
    % Store simualated control points.
    ptsToSimulate(:,loop)=...
        reshape(simCtrlPts(:,ismember(ctrlId,idsToSimulate)),[],1);

    if (verb>0), disp(' '); end
    
    if (preLoad)
        % First clear data by creating a new project.
        if (verb>0), fprintf('Clearing project: '); end
        cmd='NewProject dummy';
        [ok,num,str]=ddecmd(ch,cmd); %#ok<NASGU>
        if (~ok), ddepmterm(ch); error('Failed to clear project.'); end
        if (verb>0), fprintf('ok.\n'); end
    end
    
    if (verb>0), fprintf('Loading project: '); end
    cmd=['OpenProject ',pmrFile];
    [ok,num,str]=ddecmd(ch,cmd); %#ok<NASGU>
    if (~ok), ddepmterm(ch); error('Failed to open project.'); end

    if (verb>0), fprintf('ok.\n'); end
    
    if (verb>0), fprintf('Defining control points: '); end
    for i=1:length(idsToSimulate)
        if (verb>1 && rem(i,10)==0), fprintf('\n'); end
        if (verb>1), fprintf(' %d',idsToSimulate(i)); end
        if (pmVer==5)
            nameStr=sprintf('c%d',idsToSimulate(i));
        else
            nameStr=sprintf('imp%d c%d',idsToSimulate(i),idsToSimulate(i));
        end
        % Set control points.
        if (size(simCtrlPts,1)==3)
            cmd=sprintf('DefineControlPoint %s %.10g %.10g %.10g',...
                        nameStr,simCtrlPts(:,find(ctrlId==idsToSimulate(i))));
        else
            cmd=sprintf('DefineControlPoint %s %.10g %.10g %.10g %.10g %.10g %.10g',...
                        nameStr,simCtrlPts(:,find(ctrlId==idsToSimulate(i))));
        end
        [ok,num,str]=ddecmd(ch,cmd); %#ok<NASGU>
        if (~ok), ddepmterm(ch); error('Failed to define control point.'); end
    end
    if (verb>0), fprintf('.\n'); end
    
    if (verb>0), fprintf('Adding mark pts for control points: '); end

    for i=1:length(idsToSimulate)
        if (verb>1 && rem(i,5)==0), fprintf('\n'); end
        % Find mark points for this control point.
        mix=find(proj.markPts(:,2)==PMid(idsToSimulate(i),1));

        cmd='GetNextPointID';
        [ok,num,str]=ddecmd(ch,cmd); %#ok<NASGU>
        if (~ok), ddepmterm(ch); error('Failed to get next point ID.'); end
        markId=num{2};

        for j=mix(:)'
            % Data for this mark.
            xy=proj.markPts(j,3:4);
            ptId=idsToSimulate(i);
            imNo=proj.markPts(j,1)+1;
        
            if (pmVer==5)
                ctrlStr=sprintf('c%d',ptId);
            else
                ctrlStr=sprintf('imp%d c%d',ptId,ptId);
            end
            
            if (verb>1), fprintf(' %d %d (%d)',imNo,markId,ptId); end
            if (isempty(ctrlStr))
                cmd=sprintf('MP %d %d %.10g %.10g',imNo,markId,xy);
            else
                cmd=sprintf('MP %d %d %.10g %.10g %s',imNo,markId,xy,ctrlStr);
            end
            [ok,num,str]=ddecmd(ch,cmd); %#ok<NASGU>
            if (~ok), ddepmterm(ch); cmd, error('Failed to set mark.'); end
        end
    end
    if (verb>0), fprintf('.\n'); end

    if (orientAll)
        if (verb>0), disp('Running bundle w all orient...'); end
        cmd='Process 5';
    else
        if (verb>0), disp('Running bundle w no orient...'); end
        cmd='Process 4';
    end
    [ok,num,str]=ddecmd(ch,cmd,1000); %#ok<NASGU>
    if (~ok), ddepmterm(ch); error('Failed to run bundle.'); end
    if (verb>0)
        fprintf('Bundle iters: %d, initial error: %g, final error: %g.\n',...
                num{2:end});
    end

    % Get wanted points after bundle.
    pts=ddegetallpts(ch,PMid(idsToRecord,1));
    if (size(pts,1)~=length(idsToRecord))
        error('Not all points were calculated by PM.');
    end
    ptsToRecord(:,loop)=reshape(pts(:,2:4)',[],1);

    t=now;
    etf=(t-start)/(loop-preLoop)*(loops-preLoop)+start;
    disp(sprintf('%d loops of %d done (%d%%). Avg loop time: %s. Expected finish time: %s.',...
        loop,loops,floor(loop/loops*100),datestr((t-start)/(loop-preLoop),13),datestr(etf)));

    if (saveInterval>=0)
        if (t-lastSaveDate>datenum(0,0,0,0,saveInterval,0))
            % Time to save.
            if (lastSaveWasOdd)
                saveFile=evenFile;
            else
                saveFile=oddFile;
            end
            
            if (verb>0), disp(' '); disp(['Saving data to ',saveFile]); end
            
            saveDate=t;
            save(saveFile,'saveDate','startDate','seed','ptsToRecord',...
                'idsToRecord','ptsToSimulate','idsToSimulate','pt','ctrlId');
            
            lastSaveDate=saveDate;
            lastSaveWasOdd=~lastSaveWasOdd;            
        end
    end
end

if (verb>0), disp(' '); disp(['Saving data to ',doneFile]); end

saveDate=now;
save(doneFile,'saveDate','startDate','seed','ptsToRecord',...
                'idsToRecord','ptsToSimulate','idsToSimulate','pt','ctrlId');

% Calculate mean and covariance of recorded points.
mn=reshape(mean(ptsToRecord,2),3,[]);
C=cov(ptsToRecord');
bundlePts=PTGaussian(mn,C);

ddepmterm(ch);
