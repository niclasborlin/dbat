timeBench=false

sss=cell(1,2);
ttt=cell(1,2);
iii=cell(1,2);
allRes=cell(2,10);
donePc=zeros(size(allRes));

dataSetNames={'\full{}',
              '\thin{}',
              '\corner{}',
              '\cornerthin{}',
              '\halff{}',
              '\full{} without 5 deg pts',
              '\thin{} without 5 deg pts',
              '\corner{} without 5 deg pts',
              '\cornerthin{} without 5 deg pts',
              '\halff{} without 5 deg pts'};

fileNames=cell(1,5);

nCams=repmat([60,15,20,10,30],1,2);
nOPs=[26321,4154,7250,3988,12932,24406,4154,6226,3442,nan];
lengthX=6*nCams-7+3*nOPs;

vetoAll=false;
mirror=false;

largePert=false;
smallPert=true;

allResults=cell(2,8);

vetos=[false,false,true,false,true];
chgs=[false,true,true,false,false];
rems=[false,false,false,true,true];
for e=1:5
    veto=vetos(e);
    changeBad=chgs(e);
    removeBad=rems(e);
    
    row=e;
    
    n=8;
    results=cell(1,10);
    
    ii=1;
    for i=ii
        if veto
            fileName=['romabundle_set',int2str(i),'_veto'];
            if vetoAll
                fileName=[fileName,'_all'];
            end
        else
            fileName=['romabundle_set',int2str(i)];
        end
        
        if veto & mirror
            fileName=[fileName,'_mirror'];
        end
        
        if largePert
            fileName=[fileName,'_large'];
        end

        if smallPert
            fileName=[fileName,'_small'];
        end
        
        if removeBad
            fileName=[fileName,'_removebad'];
        end
        
        if changeBad
            fileName=[fileName,'_changebad'];
        end

        if timeBench
            fileName=[fileName,'_bench_trillian'];
        end
        
        fileName=[fileName,'.mat'];

        fileNames{e}=fileName;
        
        Z=load(fileName);

        [dummy,ia,ib]=intersect(Z.results.cPert,(0:4)/100);
        Z.results.cPert=Z.results.cPert(ia);
        Z.results.iters=Z.results.iters(:,:,ia,:);
        Z.results.time=Z.results.time(:,:,ia,:);
        Z.results.code=Z.results.code(:,:,ia,:);
        Z.results.diffNormAvg=Z.results.diffNormAvg(:,:,ia,:);
        Z.results.diffNormMax=Z.results.diffNormMax(:,:,ia,:);
        Z.results.vetox0=Z.results.vetox0(:,:,ia,:);
        if isfield(Z.results,'vetoActivated')
            Z.results.vetoActivated=Z.results.vetoActivated(:,:,ia,:);
        else
            Z.results.vetoActivated=false(size(Z.results.vetox0));
        end
        if isfield(Z.results,'badPts')
            Z.results.badPts=Z.results.badPts(:,:,ia,:);
        else
            Z.results.badPts=false(size(Z.results.vetox0));
        end
        % Load time adjustment.
        host=Z.results.host;
        if strcmp(host,'trillian')
            Z.results.refTime=1;
        else
            try
                fReal=strrep(fileName,'.mat',['_bench_',host,'.mat']);
                fRef=strrep(fileName,'.mat',['_bench_pixel.mat']);
                Zreal=load(fullfile('results',fReal));
                Zref=load(fullfile('results',fRef));
                Z.results.refTime=sum(Zref.results.time(:))/...
                    sum(Zreal.results.time(:));
                Z.results.time=Z.results.time*Z.results.refTime;
            catch
                Z.results.refTime=nan;
            end
        end
        
        results{i}=Z.results;
    end
    
    %results(6)=results(2);
    
    allRes(row,:)=results;
    
    done=zeros(size(results));
    nIter=zeros(size(results));
    for i=1:length(results)
        if ~isempty(results{i})
            % Calculate # of iterations.
            sz=size(results{i}.time);
            if length(sz)<4
                sz(end+1:4)=1;
            end
            t=permute(results{i}.time,[1,3,4,2]);
            nIter(i)=sz(2);
            done(i)=min([find(any(isnan(reshape(t,prod(sz([1,3,4])),[]))))-1,nIter(i)]);
            done(i)=min(find([isnan(results{i}.time(1,:,1,1)),251]))-1
            results{i}.time=results{i}.time(:,1:done(i),:,:);
            results{i}.iters=results{i}.iters(:,1:done(i),:,:);
            results{i}.code=results{i}.code(:,1:done(i),:,:);
            results{i}.diffNormAvg=results{i}.diffNormAvg(:,1:done(i),:,:);
            results{i}.diffNormMax=results{i}.diffNormMax(:,1:done(i),:,:);
            if removeBad
                results{i}.badPts=results{i}.badPts(:,1:done(i),:,:);
            end
            donePc(row,i)=done(i)/nIter(i);
        end
    end
    
    diffNormCutoff=1e-4;
    
    vetoX0only=cell(size(results));
    vetoActiveOnly=cell(size(results));
    vetoX0Active=cell(size(results));

    success=cell(size(results));
    avgIters=cell(size(results));
    avgTime=cell(size(results));
    for i=1:length(results)
        if ~isempty(results{i})
            sz=size(results{i}.time);
            if length(sz)<4
                sz(end+1:4)=1;
            end
            figure(tagfigure('cut'))
            subplot(6,n,(row-1)*n+i)
            dn=results{i}.diffNormAvg.*(results{i}.code==0);
            semilogy(sort(dn(:)),'.')
            line([0,prod(sz)],[diffNormCutoff,diffNormCutoff],'color','r')
            % Calculate maximum "ok" diffNorm below threshold.
            ok=results{i}.code==0 & dn<diffNormCutoff;
            vetox0=results{i}.vetox0;
            vetoActivated=results{i}.vetoActivated;
            max(results{i}.diffNormAvg(ok))
            s=zeros(sz([3,4,1]));
            v0only=zeros(sz([3,4,1]));
            vAct=zeros(sz([3,4,1]));
            v0Act=zeros(sz([3,4,1]));
            iters=zeros(sz([3,4,1]));
            times=zeros(sz([3,4,1]));
            if removeBad
                badPts=zeros(sz([3,4,1]));
            end
            for j=1:size(s,3)
                ok1=results{i}.code(j,:,:,:)==0;
                ok2=results{i}.diffNormAvg(j,:,:,:)<diffNormCutoff;
                ok3=results{i}.iters(j,:,:,:)<=20;
                ok=ok1 & ok2 & ok3;
                if any(size(ok)<size(vetox0(j,:,:,:)))
                    MM=vetox0(j,:,:,:);
                    ok(size(MM,1),size(MM,2),size(MM,3),size(MM,4))=false;
                end
                s(:,:,j)=squeeze(sum(ok,2));
                v0only(:,:,j)=squeeze(sum(ok & vetox0(j,:,:,:) & ~vetoActivated(j,:,:,:),2));
                vAct(:,:,j)=squeeze(sum(ok & ~vetox0(j,:,:,:) & vetoActivated(j,:,:,:),2));
                v0Act(:,:,j)=squeeze(sum(ok & vetox0(j,:,:,:) & vetoActivated(j,:,:,:),2));
                iters(:,:,j)=squeeze(mean(results{i}.iters(j,:,:,:),2));
                times(:,:,j)=squeeze(mean(results{i}.time(j,:,:,:),2));
                if removeBad
                    %badPts(:,:,j)=squeeze(mean(results{i}.badPts(j,:,:,:),2));
                    badPts(:,:,j)=squeeze(max(results{i}.badPts(j,:,:,:),[],2));
                end
            end
            success{i}=floor(s/done(i)*100);
            ddone=repmat(squeeze(sum(~isnan(results{i}.time(1,:,:,:)),2)),[1,1,4]);
            success{i}=round(s./ddone*100);
            success{i}=(s./ddone*100);
            vetoX0only{i}=floor(v0only./ddone*100);
            vetoActiveOnly{i}=floor(vAct./ddone*100);
            vetoX0Active{i}=floor(v0Act./ddone*100);
            
            %iters(s~=done(i))=nan;
            %times(s~=done(i))=nan;
            avgIters{i}=iters;
            avgTime{i}=times;
        end
    end
    
    sss{row}=success;
    ttt{row}=avgTime;
    iii{row}=avgIters;
end

percCutOff=74;
%percCutOff=50;
boldCutOff=101;

okCutOff=99
grayCutOff=okCutOff-0.5;

dataSets=repmat({'\full{}','\thin{}','\corner{}','\cornerthin{}','\halff{}'},1,2);
methods={'\algC','\algGNA','\algLM','\algLMP'};

if ~timeBench
    for e=1:5
        veto=vetos(e);
        changeBad=chgs(e);
        removeBad=rems(e);
    
        row=e;
        
        for ds=1:length(dataSets)
            
        if isempty(allRes{row,ds})
            continue
        end
        
        if ~veto
            nDone1=round(allRes{row,ds}.done/allRes{row,ds}.total* ...
                         allRes{row,ds}.nRuns);
            nDone2=round(allRes{row+1,ds}.done/allRes{row+1,ds}.total* ...
                         allRes{row+1,ds}.nRuns);
            if nDone1==nDone2
                s1=sprintf('$n=%d$ ',nDone1);
                s2='';
                s3='';
            else
                s1='';
                s2=sprintf(' ($n=%d$)',nDone1);
                s3=sprintf(' ($n=%d$)',nDone2);
            end
            % Write header for double-table.
            if false
                if largePert
                    f=fullfile('results',sprintf('roma_ds%d_cap_large.tex',ds));
                else
                    f=fullfile('results',sprintf('roma_ds%d_cap.tex',ds));
                end
                fid=fopen(f,'w');
                fprintf(fid,'\\caption{Convergence results for the ');
                fprintf(fid,'%s and %s-5 datasets',dataSets{ds},dataSets{ds});
                fprintf(fid,'.\n The percentage of %s runs that achieved',s1);
                fprintf(fid,' convergence for varying levels\n');
                fprintf(fid,'of angular and translational perturbation are given.\n');
                fprintf(fid,'Only values above %d percent are presented.',percCutOff);
                fprintf(fid,'}\n');
                fclose(fid);
            end
        end
        
        f=fullfile('results',...
                   sprintf('roma_ds%d_veto%d_small_cb%drb%d.tex',...
                           ds,double(veto),double(changeBad),double(removeBad)));
        
        
        fid=fopen(f,'w');
        
        fprintf(fid,'%s\n','\vspace*{2ex}');
        fprintf(fid,'%s\n','\hspace{1em}');
        
        lastRow=min(find(any(any(sss{row}{ds}>percCutOff,1),3),1,'last')+1,...
                    size(sss{row}{ds},2));
        for method=1:4
            [m,n,dummy]=size(sss{row}{ds});
            fprintf(fid,'\\begin{tabular}{@{}r@{\\hspace{\\tabspc}}|@{\\hspace{\\tabspc}}*%d{r@{\\hspace{\\tabspc}}}}\n',m);
            for i=1:m
                fprintf(fid,'%s',' & \color{white}100');
            end
            fprintf(fid,'%s\n','\\[-1ex]');
            
            fprintf(fid,'~\\hspace*{1em}\\begin{picture}(0,0)\\put(-1.75,0.25){\\mbox{\\rotatebox{90}{\\bf %s}}}\\put(4,2){\\mbox{\\%% of $D=22$~m}}\\end{picture} ',methods{method});
            
            for i=1:m
                s=int2str(round(allRes{row,1}.cPert(i)*100));
                fprintf(fid,' & %s',s);
            end
            fprintf(fid,'%s\n','\\[-1ex]\hline');
            
            for j=1:lastRow
                fprintf(fid,'%.1f%s',allRes{row,1}.angPert(j)*180/pi,...
                        '$^\circ$');
                for i=1:m
                    v=round(sss{row}{ds}(i,j,method));
                    if v<=percCutOff
                        fprintf(fid,' &  - ');
                    elseif v<=grayCutOff
                        fprintf(fid,' & \\color{gray}%3d',v);
                    elseif v>=boldCutOff
                        fprintf(fid,' & \\bf %3d',v);
                    else
                        fprintf(fid,' & %3d',v);
                    end
                end
                if j<n
                    fprintf(fid,'%s','\\[-1.3ex]');
                end
                fprintf(fid,'\n');
            end
            
            fprintf(fid,'%s\n','\end{tabular}');
            fprintf(fid,'%s\n','\hspace{1em}');
        end
        fclose(fid);
        end
        
    end 
end

grayValues=zeros(4,1);
for i=1:5
    grayValues=grayValues+squeeze(sum(sum(sss{i}{1}<grayCutOff & sss{i}{1}>=90-0.5,1),2));
end

grayValues
grayValues./grayValues(1)

donePc*100

relTimes=zeros(4,5);

if timeBench
fid=fopen(fullfile('results','romatimes.tex'),'w');
for r=1:5
    exps=[1,4,2,5,3];
    e=exps(r);
    yn={'no','yes'};
    OPs='Unchanged';
    if chgs(e)
        OPs='Modified';
    end
    if rems(e)
        OPs='Removed';
    end
    fprintf(fid,'%d & %-10s & %-3s ',r,OPs,yn{vetos(e)+1});
    fprintf(fid,'& %5.1f & %5.1f',reshape([iii{e}{1};ttt{e}{1}],[],1));
    fprintf(fid,'%s','\\');
    if r<5
        fprintf(fid,'%s\n','\hline');
    else
        fprintf(fid,'\n');
    end
    relTimes(:,r)=squeeze(ttt{e}{1}./repmat(ttt{e}{1}(:,:,1),[1,1,4]));
end
fclose(fid);
relTimesVetoRemoved=squeeze(ttt{5}{1}./ttt{4}{1})
relTimesVetoChanged=squeeze(ttt{3}{1}./ttt{2}{1})
relTimes
end

if ~timeBench
    round(badPts(:,:,1))
end