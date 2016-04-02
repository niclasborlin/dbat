aIx=1:3;
pIx=4:6;

dataDir=fileparts(resFile);

fprintf('Doing PM/DBAT comparison for dir %s.\n',dataDir);

pm_s0=load(fullfile(dataDir,'pm_s0.txt'))
dbat_s0=load(fullfile(dataDir,'dbat_s0.txt'))

% IO parameters
pmIO=load(fullfile(dataDir,'pm_IO.txt'));
dbatIO=load(fullfile(dataDir,'dbat_IO.txt'));

if isempty(pmIO) && isempty(dbatIO)
    fprintf('No IO values.\n');
else
    fprintf('IO values within a factor of %.4f.\n',...
            max(pow2(abs(log2(pmIO./dbatIO)))));
end

pmIOdev=load(fullfile(dataDir,'pm_IOdev.txt'));
dbatIOdev=load(fullfile(dataDir,'dbat_IOdev.txt'));

if isempty(pmIOdev) && isempty(dbatIOdev)
    fprintf('No IO deviations.\n');
else
    ioDev=pow2(abs(log2(pmIOdev./dbatIOdev)));
    i=find(ioDev==max(ioDev),1);
    fprintf('IO deviations within a factor of %.4f (%d: %g vs %g).\n',...
            max(ioDev),i,pmIOdev(i),dbatIOdev(i));
end

pmIOcorr=load(fullfile(dataDir,'pm_IOcorr.txt'));
dbatIOcorr=load(fullfile(dataDir,'dbat_IOcorr.txt'));

if isempty(pmIOcorr) && isempty(dbatIOcorr)
    fprintf('No IO correlations.\n');
else
    fprintf('%d IO correlations within %.1f%%-units.\n',...
            length(pmIOcorr),max(abs(pmIOcorr-dbatIOcorr)));
end

% EO parameters
pmEO=reshape(load(fullfile(dataDir,'pm_EO.txt')),6,[]);
dbatEO=reshape(load(fullfile(dataDir,'dbat_EO.txt')),6,[]);

fprintf('EO angles   within %.4e degrees.\n',...
        max(max(abs(rem(pmEO(aIx,:)-dbatEO(aIx,:),360)))));

pmEOdev=reshape(load(fullfile(dataDir,'pm_EOdev.txt')),6,[]);
dbatEOdev=reshape(load(fullfile(dataDir,'dbat_EOdev.txt')),6,[]);

pmEOadev=pmEOdev(aIx,:);
dbatEOadev=dbatEOdev(aIx,:);
aDev=pow2(abs(log2(pmEOadev./dbatEOadev)));
maxAdev=max(max(aDev));
if maxAdev>1.01
    minAdev=min(min(aDev));
    [i,j]=find(minAdev==aDev,1);
    fprintf('EO angles   deviations between a factor of %.4f ((%d,%d) %g vs %g)\n',...
            minAdev,i,j,pmEOadev(i,j),dbatEOadev(i,j));
    [i,j]=find(maxAdev==aDev,1);
    fprintf('                                       and %.4f ((%d,%d) %g vs %g)\n',...
            maxAdev,i,j,pmEOadev(i,j),dbatEOadev(i,j));

    fprintf('                                       avg %.4f.\n',mean(aDev(:)));

else
    [i,j]=find(maxAdev==aDev,1);
    fprintf('EO angles   deviations within a factor of %.4f ((%d,%d) %g vs %g).\n',...
            maxAdev,i,j,pmEOadev(i,j),dbatEOadev(i,j));
end

fprintf('EO positions within %.4e project units.\n',...
        max(max(abs(pmEO(pIx,:)-dbatEO(pIx,:)))));

pmEOpdev=pmEOdev(pIx,:);
dbatEOpdev=dbatEOdev(pIx,:);
pDev=pow2(abs(log2(pmEOpdev./dbatEOpdev)));
maxPdev=max(max(pDev));
if maxPdev>1.01
    minPdev=min(min(pDev));
    [i,j]=find(minPdev==pDev,1);
    fprintf('EO position deviations between a factor of %.4f ((%d,%d) %g vs %g)\n',...
            minPdev,i,j,pmEOpdev(i,j),dbatEOpdev(i,j));
    [i,j]=find(maxPdev==pDev,1);
    fprintf('                                       and %.4f ((%d,%d) %g vs %g)\n',...
            maxPdev,i,j,pmEOpdev(i,j),dbatEOpdev(i,j));

    fprintf('                                       avg %.4f.\n',mean(pDev(:)));

else
    [i,j]=find(maxPdev==pDev,1);
    fprintf('EO position deviations within a factor of %.4f ((%d,%d) %g vs %g).\n',...
            maxPdev,i,j,pmEOpdev(i,j),dbatEOpdev(i,j));
end

pmEOcorr=load(fullfile(dataDir,'pm_EOcorr.txt'));
dbatEOcorr=load(fullfile(dataDir,'dbat_EOcorr.txt'));

if isempty(pmEOcorr) && isempty(dbatEOcorr)
    fprintf('No EO correlations.\n');
else
    corrDiff=abs(pmEOcorr-dbatEOcorr);
    maxCorr=max(corrDiff);
    if maxCorr>0.1
        minCorr=min(corrDiff);
        i=find(corrDiff==minCorr,1);
        fprintf('%d EO correlations between %.1f%%-units ((%d) %g vs %g)',...
                length(pmEOcorr),minCorr,i,pmEOcorr(i),dbatEOcorr(i));
        i=find(corrDiff==maxCorr,1);
        fprintf(' and %.1f%%-units ((%d) %g vs %g).\n',...
                maxCorr,i,pmEOcorr(i),dbatEOcorr(i));
    else
        fprintf('%d EO correlations within %.1f%%-units.\n',...
                length(pmEOcorr),maxCorr);
    end
end

% OP parameters
fprintf('OP+CP positions within %.4g project units.\n',...
        max(max(abs(ss0.OP-result{1}.OP))));

dbatOPstd=reshape(full(sqrt(diag(COP))),3,[]);
pmOPstd=ss0.OPstd;

dbatOPstdOP=dbatOPstd(:,~s0.isCtrl);
pmOPstdOP=pmOPstd(:,~s0.isCtrl);
dbatOPstdCP=dbatOPstd(:,s0.isCtrl);
pmOPstdCP=pmOPstd(:,s0.isCtrl);

OPdiff=pow2(abs(log2(pmOPstdOP./dbatOPstdOP)));
CPdiff=pow2(abs(log2(pmOPstdCP./dbatOPstdCP)));

if isempty(OPdiff)
    fprintf('No OP.\n');
else
    maxOPdiff=max(OPdiff(:));
    if maxOPdiff>1.01
        minOPdiff=min(OPdiff(:));
        [i,j]=find(OPdiff==minOPdiff,1);
        fprintf('OP position deviations between a factor of %.4f ((%d,%d) %g vs %g)\n',...
                minOPdiff,i,j,pmOPstdOP(i,j),dbatOPstdOP(i,j));
        [i,j]=find(OPdiff==maxOPdiff,1);
        fprintf('                                       and %.4f ((%d,%d) %g vs %g)\n',...
                maxOPdiff,i,j,pmOPstdOP(i,j),dbatOPstdOP(i,j));
        fprintf('                                       avg %.4f.\n',mean(OPdiff(:)));
    else
        [i,j]=find(OPdiff==maxOPdiff,1);

        fprintf('OP position deviations within a factor of %.4f ((%d,%d) %g vs %g).\n',...
                maxOPdiff,i,j,pmOPstdOP(i,j),dbatOPstdOP(i,j));
    end
end

if isempty(CPdiff)
    fprintf('No OP.\n');
else
    maxCPdiff=max(CPdiff(:));
    if maxCPdiff>1.01
        minCPdiff=min(CPdiff(:));
        [i,j]=find(CPdiff==minCPdiff,1);
        fprintf('CP position deviations between a factor of %.4f ((%d,%d) %g vs %g)\n',...
                minCPdiff,i,j,pmOPstdCP(i,j),dbatOPstdCP(i,j));
        [i,j]=find(CPdiff==maxCPdiff,1);
        fprintf('                                       and %.4f ((%d,%d) %g vs %g)\n',...
                maxCPdiff,i,j,pmOPstdCP(i,j),dbatOPstdCP(i,j));
        fprintf('                                       avg %.4f.\n',mean(CPdiff(:)));
    else
        [i,j]=find(CPdiff==maxCPdiff,1);
        
        fprintf('CP position deviations within a factor of %.4f ((%d,%d) %g vs %g).\n',...
                maxCPdiff,i,j,pmOPstdCP(i,j),dbatOPstdCP(i,j));
    end
end

pmPts2d=load(fullfile(dataDir,'2dpts-cleaned.txt'));
pmPts=nan(4,size(pmPts2d,1));
for i=1:size(pmPts2d,1)
    id=pmPts2d(i,1);
    imNo=pmPts2d(i,2);
    imNo=min(imNo,size(s.colPos,2));
    pos=s.colPos(id==s.OPid,imNo);
    pmPts(:,pos)=pmPts2d(i,3:6)';
end

r=E{1}.final.weighted.r;
J=E{1}.final.weighted.J;

fprintf('markpts   within %.2g pixels.\n',...
        max(max(abs(pmPts(1:2,:)-s.markPts))));

res=reshape(r(1:end-nnz(s.useOPobs)),2,[]).*result{1}.markStd;

fig=tagfigure('markpts');
ax=gca(fig);
plot(ax,pmPts(1,:),pmPts(2,:),'rx',s.markPts(1,:),s.markPts(2,:),'bo')
axis(ax,'equal')

fig=tagfigure('residuals');
ax=gca(fig);
plot(ax,pmPts(3,:),pmPts(4,:),'rx',res(1,:),res(2,:),'bo')
axis(ax,'equal','square')
legend(ax,'PM','DBAT')

fprintf('residuals within %.2e pixels.\n',...
        max(max(abs(pmPts(3:4,:)-res))));

fprintf('\nrms mark pts res PM: %g, DBAT: %g.\n',...
        sqrt(mean(reshape(pmPts(3:4,:).^2,[],1))),sqrt(mean(res(:).^2)));

fprintf('\nsigma0 based on PM res and red=m-n: %4g.\n',...
        sqrt(r'*r/(size(J,1)-size(J,2))));
fprintf('sigma0 based on DBAT res and red=m-n: %4g.\n',...
        sqrt(r'*r/(size(J,1)-size(J,2))));
p=nnz(s.estOP(:,any(s.vis,2))==0)+nnz(s.estEO(1:6,any(s.vis,1))==0);
fprintf('sigma0 based on DBAT res and red=m-n+3CP: %4g.\n',...
        sqrt(r'*r/(size(J,1)-size(J,2)+p)));
