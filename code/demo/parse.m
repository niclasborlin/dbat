aIx=1:3;
pIx=4:6;

dataDir=fileparts(resFile);

fprintf('Doing PM/DBAT comparison for dir %s.\n',dataDir);

pm_s0=load(fullfile(dataDir,'pm_s0.txt'))
dbat_s0=load(fullfile(dataDir,'dbat_s0.txt'))

% IO parameters
pmIO=load(fullfile(dataDir,'pm_IO.txt'));
dbatIO=load(fullfile(dataDir,'dbat_IO.txt'));

fprintf('IO values within a factor of %.4f.\n',...
        max(pow2(abs(log2(pmIO./dbatIO)))));

pmIOdev=load(fullfile(dataDir,'pm_IOdev.txt'));
dbatIOdev=load(fullfile(dataDir,'dbat_IOdev.txt'));

ioDev=pow2(abs(log2(pmIOdev./dbatIOdev)));
i=find(ioDev==max(ioDev));
fprintf('IO deviations within a factor of %.4f (%d: %g vs %g).\n',...
        max(ioDev),i,pmIOdev(i),dbatIOdev(i));

pmIOcorr=load(fullfile(dataDir,'pm_IOcorr.txt'));
dbatIOcorr=load(fullfile(dataDir,'dbat_IOcorr.txt'));

if isempty(pmIOcorr) && isempty(dbatIOcorr)
    fprintf('No IO correlations.\n');
else
    fprintf('IO correlations within %.2f%%-units.\n',...
            max(abs(pmIOcorr-dbatIOcorr)));
end

% EO parameters
pmEO=reshape(load(fullfile(dataDir,'pm_EO.txt')),6,[]);
dbatEO=reshape(load(fullfile(dataDir,'dbat_EO.txt')),6,[]);

fprintf('EO angles    within %.4g degrees.\n',...
        max(max(abs(rem(pmEO(aIx,:)-dbatEO(aIx,:),360)))));
fprintf('EO positions within %.4g project units.\n',...
        max(max(abs(pmEO(pIx,:)-dbatEO(pIx,:)))));

pmEOdev=reshape(load(fullfile(dataDir,'pm_EOdev.txt')),6,[]);
dbatEOdev=reshape(load(fullfile(dataDir,'dbat_EOdev.txt')),6,[]);

pmEOadev=pmEOdev(aIx,:);
dbatEOadev=dbatEOdev(aIx,:);
aDev=pow2(abs(log2(pmEOadev./dbatEOadev)));
maxAdev=max(max(aDev));
[i,j]=find(maxAdev==aDev);

fprintf('EO angles   deviations within a factor of %.4f ((%d,%d) %g vs %g).\n',...
        maxAdev,i,j,pmEOadev(i,j),dbatEOadev(i,j));
        
pmEOpdev=pmEOdev(pIx,:);
dbatEOpdev=dbatEOdev(pIx,:);
pDev=pow2(abs(log2(pmEOpdev./dbatEOpdev)));
maxPdev=max(max(pDev));
[i,j]=find(maxPdev==pDev);

fprintf('EO position deviations within a factor of %.4f ((%d,%d) %g vs %g).\n',...
        maxPdev,i,j,pmEOpdev(i,j),dbatEOpdev(i,j));

pmEOcorr=load(fullfile(dataDir,'pm_EOcorr.txt'));
dbatEOcorr=load(fullfile(dataDir,'dbat_EOcorr.txt'));

if isempty(pmEOcorr) && isempty(dbatEOcorr)
    fprintf('No EO correlations.\n');
else
    fprintf('EO correlations within %.2f%%-units.\n',...
            max(abs(pmEOcorr-dbatEOcorr)));
end

% OP parameters
fprintf('OP positions within %.4g project units.\n',...
        max(max(abs(ss0.OP-result{1}.OP))));

dbatOPstd=reshape(full(sqrt(diag(COP))),3,[]);
pmOPstd=ss0.OPstd;

dbatOPstdOP=dbatOPstd(:,~s0.isCtrl);
pmOPstdOP=pmOPstd(:,~s0.isCtrl);
dbatOPstdCP=dbatOPstd(:,s0.isCtrl);
pmOPstdCP=pmOPstd(:,s0.isCtrl);

OPdiff=pow2(abs(log2(pmOPstdOP./dbatOPstdOP)));
CPdiff=pow2(abs(log2(pmOPstdCP./dbatOPstdCP)));

maxOPdiff=max(OPdiff(:));
[i,j]=find(OPdiff==maxOPdiff);

fprintf('OP position deviations within a factor of %.4f ((%d,%d) %g vs %g).\n',...
        maxOPdiff,i,j,pmOPstdOP(i,j),dbatOPstdOP(i,j));

maxCPdiff=max(CPdiff(:));
[i,j]=find(CPdiff==maxCPdiff);
i=i(1);
j=j(i);
fprintf('CP position deviations within a factor of %.4f ((%d,%d) %g vs %g).\n',...
        maxCPdiff,i,j,pmOPstdCP(i,j),dbatOPstdCP(i,j));

