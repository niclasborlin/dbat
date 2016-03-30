aIx=1:3;
pIx=4:6;

pmVal=reshape(load('pm_values.txt'),6,[]);
dbatVal=reshape(load('dbat_values.txt'),6,[]);

pmDev=reshape(load('pm_dev.txt'),6,[]);
dbatDev=reshape(load('dbat_dev.txt'),6,[]);

angDiff=max(max(abs(pmVal(aIx,:)-dbatVal(aIx,:))))
camPosDiff=max(max(abs(pmVal(pIx,:)-dbatVal(pIx,:))))