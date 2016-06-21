function [ptsToSimulate,ptsToRecord,preLoop,startDate,seed,lastSaveWasOdd,files]=loadsimdata(year,noFWpts,mode,ids,idr)
%LOADSIMDATA Load simulation data from checkpoint or final files.
%
%[pts,ptr,preLoop,sDate,seed,odd,files]=loadsimdata(year,noFWpts,mode,ids,idr)
%or
%data=loadsimdata(year,noFWpts,mode);
%year    - which simulation year.
%noFWpts - if true, use files without any forward intersection points
%          (significantly faster).
%mode    - 'i' - interactive mode, i.e. ask user which file to load.
%          'n' - start new without loading anything.
%          'd' - load from 'done' file, else fail.
%          'a' - auto, load latest checkpoint file if simulation
%                parameters match, otherwise revert to interactive mode.
%ids     - vector of IDs to simulate.
%idr     - vector of IDs to record.
%pts     - 3-by-N array containing simulated control points.
%ptr     - 3-by-N array containing recorded control points.
%preLoop - the number of loops already performed.
%sDate   - the timestamp of the start of the loaded simulation.
%seed    - the random seed of the loaded simulation.
%odd     - true if final recorded save was to the 'odd' checkpoint file.
%files   - cell array with names of {odd,even,done} file.
%data    - struct with data loaded.

% $Id$

loadOnly=nargin<4;

if ~ismember(year,[1990,2007])
    error('Bad year');
end

if (~ismember(lower(mode(1)),'inad'))
    error('Bad mode');
end
  
if noFWpts
    % Checkpoint files.
    oddFile=sprintf('tmpdata_%d_nofw_odd.mat',year);
    evenFile=sprintf('tmpdata_%d_nofw_even.mat',year);

    % File for final result.

    doneFile=sprintf('done_%d_nofw.mat',year);
else
    % Checkpoint files.
    oddFile=sprintf('tmpdata_%d_odd.mat',year);
    evenFile=sprintf('tmpdata_%d_even.mat',year);

    % File for final result.
    doneFile=sprintf('done_%d.mat',year);
end

files={oddFile,evenFile,doneFile};

nOdd=0;
nEven=0;
oddData=[];
evenData=[];
doneData=[];
if ismember(lower(mode(1)),'ia')
    if loadOnly
        [oddData,nOdd,bad]=loadfile(oddFile,'odd');
        [evenData,nEven,bad]=loadfile(evenFile,'even');
    else
        [oddData,nOdd,bad]=loadfile(oddFile,'odd',ids,idr);
        if bad, mode='i'; end
        [evenData,nEven,bad]=loadfile(evenFile,'even',ids,idr);
        if bad, mode='i'; end
    end
end

if ismember(lower(mode(1)),'iad')
    if loadOnly
        [doneData,nDone,bad]=loadfile(doneFile,'done');
    else
        [doneData,nDone,bad]=loadfile(doneFile,'done',ids,idr);
        if nDone>0 && bad && lower(mode(1))=='a'
            % Revert to interactive unless nDone==0, i.e. no done file.
            mode='i';
        end
    end
end

switch lower(mode(1))
case 'n' % New
    r='n';
case 'd' % Done
    r='d';
case 'i' % Interactive
    r=' ';
    while ~ismember(r,'oden')
        disp(' ');
        r=input(['Do you want to load stored (o)dd, (e)ven, (d)one data or ' ...
                 'start (n)ew? '],'s');
    end
case 'a' % Auto
    [dummy,i]=max([nOdd,nEven,nDone]);
    ss='oed';
    r=ss(i);
end

if r=='n'
    disp('Setting up new data');
    seed=28; % Change this number to get another random sequence.
    startDate=now;
    preLoop=0;
    ptsToRecord=zeros(length(idr)*3,1);
    ptsToSimulate=zeros(length(ids)*3,1);
    lastSaveWasOdd=false;
else
    dd={oddData,evenData,doneData};
    ss='oed';
    strs={'odd','even','done'};
    disp(' ')
    disp(['Using ',strs{find(ss==r)},' data']);

    data=dd{find(ss==r)};

    if isempty(data)
        error('Failed to load')
    end
    
    if loadOnly
        %data=loadsimdata(year,noFWpts);
        ptsToSimulate=data;
    else
        seed=data.seed;
        startDate=data.startDate;
        preLoop=nnz(any(data.ptsToRecord));
        ptsToRecord=data.ptsToRecord;
        ptsToSimulate=data.ptsToSimulate;
        lastSaveWasOdd=r=='o';
    end
end

function [data,n,badId]=loadfile(fileName,str,ids,idr)
%fileName - name of file to load.
%str      - descriptive string, e.g. 'odd','even','done'.
%ids      - ids to simulate.
%idr      - ids to record.
%data     - struct with loaded data
%n        - number of stored loops.
%badId    - true if supplied and stored ids do not match.

checkIds=nargin>=3;
data=[];
badId=true;
n=0;

if (exist(fileName,'file'))
    try
        data=load(fileName);
    catch
    end
end

if (isempty(data))
    disp(' ');
    disp(['No data was found in ',str,' file.']);
else
    s=data.startDate;
    d=data.saveDate;
    n=nnz(any(data.ptsToRecord));
    disp(' ');
    disp(['Simulation in ',str,' file ''',fileName,''' started ',datestr(s)]);
    disp(['and was saved ',datestr(d),' after ',...
          int2str(n),' iterations.']);
    if checkIds
        try
            bad=~isempty(setxor(idr,data.idsToRecord)) | ...
                ~isempty(setxor(ids,data.idsToSimulate));
            if ~bad, badId=false; end
        catch
            bad=true;
        end
        if bad
            disp(['Warning: Different ids in ',str,' file!']);
        end
    else
        bad=false;
    end
end
