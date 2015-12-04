doPrint=true;

for year=[1990,2007]
    [idsToSimulate,idsToRecord,ctrlId,pt,ptPMerr,PMid,allId,allPMerr]=load_cpt(year,true);
    data=loadsimdata(year,true,'a');
    
    if year==1990
        pt0=pt(:,ctrlId==0);
        pt0(isnan(pt0.mean))=0;
        pt0=PTGaussian(pt0.mean);
        data90=data;
        ctrlId90=ctrlId;
        PMid90=PMid;
        ptPMerr90=ptPMerr;
        allId90=allId;
        allPMerr90=allPMerr;
    else
        data07=data;
        ctrlId07=ctrlId;
        PMid07=PMid;
        ptPMerr07=ptPMerr;
        allId07=allId;
        allPMerr07=allPMerr;
    end

    n=nnz(any(data.ptsToRecord));

    idsToPlot=idsToSimulate;

    % Find given control points.
    [givenSimid,ia,ib]=intersect(ctrlId,idsToPlot);
    % Given points.
    givenSim=pt(:,ia);

    % Find simulated points.
    [dummy,ia,ib]=intersect(data.idsToSimulate,idsToPlot);
    % All simulated points.
    allSim=reshape(data.ptsToSimulate,3,[],size(data.ptsToSimulate,2));
    % Extract the wanted points.
    sim=allSim(:,ia,1:n);
    % Create PTGaussians from mean and covariance of samples.
    simpt=PTGaussian(mean(sim,3),mkblkdiag(cov(reshape(sim,[],n)'),3));
    simptFull=PTGaussian(mean(sim,3),cov(reshape(sim,[],n)'));

    mag=5000;
    fig=tagfigure(sprintf('sim%d',year));
    h=plotgivensimerr(fig,year,mag,givenSim,givenSimid,simpt);
    title(sprintf(['Planimetric error ellipses for ',...
                   'input ctrl pts, year=%d, n=%d, mag=%d'],year,n,mag));
    xlim(visconfreg(givenSim(1,givenSimid==7),2*mag));
    ylim(visconfreg(givenSim(2,givenSimid==7),2*mag));

    pt1=givenSim;
    pt2full=simptFull;
    pt2=simpt;

    isExact=any(pt1.std==0 | pt2.std==0);
    pt1=pt1(:,~isExact);
    pt2=pt2(:,~isExact);
    pt2full=pt2full(:,~isExact);
    size(pt1,2)

    disp(sprintf('p value for whole vector: %.3g',probdifferent(pt1(:),pt2(:))))
    disp(sprintf('p value for whole vector: %.3g',probdifferent(pt1(:),pt2full(:))))
    
    idsToPlot=idsToRecord;

    % Find given control points.
    [givenRecid,ia,ib]=intersect(ctrlId,idsToPlot);
    % Given points.
    givenRec=pt(:,ia);
    % PM errors for given points.
    givenRecPMerr=ptPMerr(:,ia);
    
    % Find simulated points.
    [dummy,ia,ib]=intersect(data.idsToRecord,idsToPlot);
    % All simulated points.
    allRec=reshape(data.ptsToRecord,3,[],size(data.ptsToRecord,2));
    % Extract the wanted points.
    rec=allRec(:,ia,1:n);
    % Create PTGaussians from mean and covariance for wanted points.
    recpt=PTGaussian(mean(rec,3),mkblkdiag(cov(reshape(rec,[],n)'),3));
    recptFull=PTGaussian(mean(rec,3),cov(reshape(rec,[],n)'));

    pt1=givenRec;
    pt2full=recptFull+givenRecPMerr;
    pt2=recpt+givenRecPMerr;

    isExact=any(pt1.std==0 | pt2.std==0);
    pt1=pt1(:,~isExact);
    pt2=pt2(:,~isExact);
    pt2full=pt2full(:,~isExact);
    zdOut=zeros(1,size(pt1,2));
    for i=1:length(zdOut)
        zdOut(i)=probdifferent(pt1(:,i),pt2(:,i));
    end
    length(zdOut)
    for alpha=[0.95,0.99]
        disp(sprintf('Number of points with p>%.2f=%d (expected=%.1f)',...
                     alpha,nnz(zdOut>alpha),(1-alpha)*length(zdOut)));
    end

    disp(sprintf('p value for whole vector: %.3g',probdifferent(pt1(:),pt2(:))))
    disp(sprintf('p value for whole vector: %.3g',probdifferent(pt1(:),pt2full(:))))

    if 0
        if year==2007
            % Test by leaving one point out.
            for i=1:size(pt1,2)
                p1=probdifferent(reshape(pt1(:,1:size(pt1,2)~=i),[],1),...
                                 reshape(pt2(:,1:size(pt1,2)~=i),[],1));
                p2=probdifferent(reshape(pt1(:,1:size(pt1,2)~=i),[],1),...
                                 reshape(pt2full(:,1:size(pt1,2)~=i),[],1));
                disp(sprintf('p value for whole vector sans point %d: %.3g (%.3g)',...
                             i,p1,p2));
            end
        end
    end

    mag=1000;
    fig=tagfigure(sprintf('rec%d',year));
    h=[h,plotgivensimerr(fig,year,mag,givenRec,givenRecid,recpt+givenRecPMerr)];
    title(sprintf(['Planimetric error ellipses for ',...
                   'output ctrl pts, year=%d, n=%d, mag=%d'],year,n,mag));
end

if doPrint
    print(sprintf('-f%d',tagfigure(sprintf('sim%d',1990))),'-depsc',...
          sprintf('edfdata/report/err_inp07.eps',year));
end

% Find common points
simId=intersect(data90.idsToSimulate,data07.idsToSimulate);

% Extract 1990 points.
n90=nnz(any(data90.ptsToRecord));
[dummy,ia]=intersect(data90.idsToSimulate,simId);
% All simulated points.
allSim=reshape(data90.ptsToSimulate,3,[],size(data90.ptsToSimulate,2));
% Extract the wanted points.
sim=allSim(:,ia,1:n90);
% Create PTGaussians from mean and covariance of samples.
sim90=PTGaussian(mean(sim,3),mkblkdiag(cov(reshape(sim,[],n90)'),3));

% Display aggregate uncertainties
disp(' ')
disp('1990 surveying points statistics:')
disp(sprintf('Number of points: %d',size(sim90,2)))
for i=1:3
    disp(sprintf('RMV_%c  : %.1f mm',abs('X')-1+i,...
                 sqrt(mean(sim90.var(i,:)))*1000))
end
disp(sprintf('RMV_XY : %.1f mm',sqrt(mean(sum(sim90.var(1:2,:),1)))*1000));
disp(sprintf('RMV_XYZ: %.1f mm',sqrt(mean(sum(sim90.var,1)))*1000));

% Extract 2007 points.
n07=nnz(any(data07.ptsToRecord));
[dummy,ia]=intersect(data07.idsToSimulate,simId);
% All simulated points.
allSim=reshape(data07.ptsToSimulate,3,[],size(data07.ptsToSimulate,2));
% Extract the wanted points.
sim=allSim(:,ia,1:n07);
% Create PTGaussians from mean and covariance of samples.
sim07=PTGaussian(mean(sim,3),mkblkdiag(cov(reshape(sim,[],n07)'),3));

% Display aggregate uncertainties
disp(' ')
disp('2007 surveying points statistics:')
disp(sprintf('Number of points: %d',size(sim07,2)))
for i=1:3
    disp(sprintf('RMV_%c  : %.1f mm',abs('X')-1+i,...
                 sqrt(mean(sim07.var(i,:)))*1000))
end
disp(sprintf('RMV_XY : %.1f mm',sqrt(mean(sum(sim07.var(1:2,:),1)))*1000));
disp(sprintf('RMV_XYZ: %.1f mm',sqrt(mean(sum(sim07.var,1)))*1000));

% Find common points
recId=intersect(data90.idsToRecord,data07.idsToRecord);

% Extract 1990 points.
n90=nnz(any(data90.ptsToRecord));
[id90,ia]=intersect(data90.idsToRecord,recId);
% All recorded points.
allRec=reshape(data90.ptsToRecord,3,[],size(data90.ptsToRecord,2));
% Extract the wanted points.
rec=allRec(:,ia,1:n90);
% Create PTGaussians from mean and covariance of samples.
rec90=PTGaussian(mean(rec,3),mkblkdiag(cov(reshape(rec,[],n90)'),3));
% Corresponding errors from PM.
rec90PMerr=allPMerr90(:,id90);
% Add errors from PM.
rec90=rec90+rec90PMerr;

% Display aggregate uncertainties
disp(' ')
disp('1990 photogrammetric points statistics:')
disp(sprintf('Number of points: %d',size(rec90,2)))
for i=1:3
    disp(sprintf('RMV_%c  : %.1f mm',abs('X')-1+i,...
                 sqrt(mean(rec90.var(i,:)))*1000))
end
disp(sprintf('RMV_XY : %.1f mm',sqrt(mean(sum(rec90.var(1:2,:),1)))*1000));
disp(sprintf('RMV_XYZ: %.1f mm',sqrt(mean(sum(rec90.var,1)))*1000));


% Extract 2007 points.
n07=nnz(any(data07.ptsToRecord));
[id07,ia]=intersect(data07.idsToRecord,recId);
% All recorded points.
allRec=reshape(data07.ptsToRecord,3,[],size(data07.ptsToRecord,2));
% Extract the wanted points.
rec=allRec(:,ia,1:n07);
% Create PTGaussians from mean and covariance of samples.
rec07=PTGaussian(mean(rec,3),mkblkdiag(cov(reshape(rec,[],n07)'),3));
% Corresponding errors from PM.
rec07PMerr=allPMerr07(:,id07);
% Add errors from PM.
rec07=rec07+rec07PMerr;

% Display aggregate uncertainties
disp(' ')
disp('2007 photogrammetric points statistics:')
disp(sprintf('Number of points: %d',size(rec07,2)))
for i=1:3
    disp(sprintf('RMV_%c  : %.1f mm',abs('X')-1+i,...
                 sqrt(mean(rec07.var(i,:)))*1000))
end
disp(sprintf('RMV_XY : %.1f mm',sqrt(mean(sum(rec07.var(1:2,:),1)))*1000));
disp(sprintf('RMV_XYZ: %.1f mm',sqrt(mean(sum(rec07.var,1)))*1000));


% Display motion uncertainties.

if n90==n07
    nStr=sprintf('n=%d',n90);
else
    nStr=sprintf('n90=%d, n07=%d',n90,n07);
end

mag=500;

fig=tagfigure('err902d');
figure(fig)
clf
plot(gca(fig),sim90(1:2,:),mag,'color','r')
hold on
plot(gca(fig),rec90(1:2,:),mag,'color','b')
axis equal
text(sim90.mean(1,:),sim90.mean(2,:),int2str(simId'), ...
     'fontsize',6,...
     'horizontal','center','color','r','parent',gca(fig));
text(rec90.mean(1,:),rec90.mean(2,:),int2str(recId'), ...
     'fontsize',6,...
     'horizontal','center','color','b','parent',gca(fig));

if doPrint
    print(fig,'-depsc','edfdata/report/err902d.eps');
end

fig=tagfigure('err903din');
figure(fig)
clf
plot(gca(fig),sim90,mag,'facecolor','r','edgecolor','r')
axis vis3d
set(gca,'projection','perspective')

if doPrint
    print(fig,'-depsc','edfdata/report/err903din.eps');
end

fig=tagfigure('err903dout');
figure(fig)
clf
plot(gca(fig),rec90,mag,'facecolor','b','edgecolor','k')
axis vis3d
set(gca,'projection','perspective')

if doPrint
    print(fig,'-depsc','edfdata/report/err903dout.eps');
end


fig=tagfigure('err072d');
figure(fig)
clf
plot(gca(fig),sim07(1:2,:),mag,'color','r')
hold on
plot(gca(fig),rec07(1:2,:),mag,'color','b')
axis equal
text(sim07.mean(1,:),sim07.mean(2,:),int2str(simId'), ...
     'fontsize',6,...
     'horizontal','center','color','r','parent',gca(fig));
text(rec07.mean(1,:),rec07.mean(2,:),int2str(recId'), ...
     'fontsize',6,...
     'horizontal','center','color','b','parent',gca(fig));

if doPrint
    print(fig,'-depsc','edfdata/report/err072d.eps');
end

fig=tagfigure('err073din');
figure(fig)
clf
plot(gca(fig),sim07,mag,'facecolor','r','edgecolor','r')
axis vis3d
set(gca,'projection','perspective')

if doPrint
    print(fig,'-depsc','edfdata/report/err073din.eps');
end

fig=tagfigure('err073dout');
figure(fig)
clf
plot(gca(fig),rec07,mag,'facecolor','b','edgecolor','k')
axis vis3d
set(gca,'projection','perspective')

if doPrint
    print(fig,'-depsc','edfdata/report/err073dout.eps');
end

simMtn=PTGaussian(sim90.mean,sim90.cov+sim07.cov);
recMtn=PTGaussian(rec90.mean,rec90.cov+rec07.cov);
fig=tagfigure('errmtn2d');

figure(fig)
clf
plot(gca(fig),simMtn(1:2,:),mag,'color','r')
hold on
plot(gca(fig),recMtn(1:2,:),mag,'color','b')
axis equal
text(sim07.mean(1,:),sim07.mean(2,:),int2str(simId'), ...
     'fontsize',6,...
     'horizontal','center','color','r','parent',gca(fig));
text(rec07.mean(1,:),rec07.mean(2,:),int2str(recId'), ...
     'fontsize',6,...
     'horizontal','center','color','b','parent',gca(fig));

if doPrint
    print(fig,'-depsc','edfdata/report/errMtn2d.eps');
end

fig=tagfigure('errMtn3din');
figure(fig)
clf
plot(gca(fig),simMtn,mag,'facecolor','r','edgecolor','r')
axis vis3d
set(gca,'projection','perspective')

if doPrint
    print(fig,'-depsc','edfdata/report/errMtn3din.eps');
end

fig=tagfigure('errMtn3dout');
figure(fig)
clf
plot(gca(fig),recMtn,mag,'facecolor','b','edgecolor','k')
axis vis3d
set(gca,'projection','perspective')

if doPrint
    print(fig,'-depsc','edfdata/report/errMtn3dout.eps');
end

mag=500;
fig=tagfigure('arrowIn');
figure(fig)
clf
simDiff=sim07-sim90;
h=arrow(sim90.mean,sim90.mean+simDiff.mean*mag,'length',5);
view(3)
axis equal
axis vis3d
set(gca,'projection','perspective')
arrow(sim90.mean,sim90.mean+simDiff.mean*mag,'facecolor','r','edgecolor','r','length',5)
delete(h);
xlabel x
ylabel y
zlabel z
text(sim90.mean(1,:),sim90.mean(2,:),sim90.mean(3,:),int2str(simId'), ...
     'fontsize',6,...
     'horizontal','center');
title(sprintf(['Arrow plot of change of input points, magnified %d times ' ...
               '(%s)'],mag,nStr));
if doPrint
    title(gca(fig),'')
    print(fig,'-depsc','edfdata/report/arrowIn.eps')
end

disp(' ')
disp('Position difference uncertainty of input surveying points:')
disp(sprintf('Number of points: %d',size(simDiff,2)))
for i=1:3
    disp(sprintf('RMV_%c  : %.1f mm',abs('X')-1+i,...
                 sqrt(mean(simDiff.var(i,:)))*1000))
end
disp(sprintf('RMV_XY : %.1f mm',sqrt(mean(sum(simDiff.var(1:2,:),1)))*1000));
disp(sprintf('RMV_XYZ: %.1f mm',sqrt(mean(sum(simDiff.var,1)))*1000));

mag=500;
fig=tagfigure('arrowOut');
figure(fig)
clf
recDiff=rec07-rec90;
h=arrow(rec90.mean,rec90.mean+recDiff.mean*mag,'length',5);
view(3)
axis equal
axis vis3d
set(gca,'projection','perspective')
arrow(rec90.mean,rec90.mean+recDiff.mean*mag,'facecolor','b','edgecolor','b','length',5)
delete(h);
xlabel x
ylabel y
zlabel z
text(rec90.mean(1,:),rec90.mean(2,:),rec90.mean(3,:),int2str(recId'), ...
     'fontsize',6,...
     'horizontal','center');
title(sprintf(['Arrow plot of change of output points, magnified %d times ' ...
               '(%s)'],mag,nStr));
if doPrint
    title(gca(fig),'')
    print(fig,'-depsc','edfdata/report/arrowOut.eps')
end

disp(' ')
disp('Position difference uncertainty of output photogrammetric points:')
disp(sprintf('Number of points: %d',size(recDiff,2)))
for i=1:3
    disp(sprintf('RMV_%c  : %.1f mm',abs('X')-1+i,...
                 sqrt(mean(recDiff.var(i,:)))*1000))
end
disp(sprintf('RMV_XY : %.1f mm',sqrt(mean(sum(recDiff.var(1:2,:),1)))*1000));
disp(sprintf('RMV_XYZ: %.1f mm',sqrt(mean(sum(recDiff.var,1)))*1000));


fig=tagfigure('arrowInOut');
figure(fig)
clf
h1=arrow(sim90.mean,sim90.mean+simDiff.mean*mag,'length',5);
h2=arrow(rec90.mean,rec90.mean+recDiff.mean*mag,'length',5);
view(3)
axis equal
axis vis3d
set(gca,'projection','perspective')
arrow(sim90.mean,sim90.mean+simDiff.mean*mag,'facecolor','r','edgecolor','r','length',5)
arrow(rec90.mean,rec90.mean+recDiff.mean*mag,'facecolor','b','edgecolor','b','length',5)
delete(h1);
delete(h2);
xlabel x
ylabel y
zlabel z
text(sim90.mean(1,:),sim90.mean(2,:),sim90.mean(3,:),int2str(simId'), ...
     'fontsize',6,...
     'horizontal','center');
text(rec90.mean(1,:),rec90.mean(2,:),rec90.mean(3,:),int2str(recId'), ...
     'fontsize',6,...
     'horizontal','center');
title(sprintf(['Arrow plot of change of input (red) and output (blue) ' ...
               'points, magnified %d times (%s)'],mag,nStr));

if doPrint
    title(gca(fig),'')
    print(fig,'-depsc','edfdata/report/arrowBoth.eps')
end

selStrs={'XZ','YZ','XY'};
sels={[1,3],[2,3],[1,2]};
for i=1:length(sels)
    sel=sels{i};
    fig=tagfigure(['arrowIn',selStrs{i}]);
    figure(fig)
    clf
    h=arrow(sim90.mean(sel,:),sim90.mean(sel,:)+simDiff.mean(sel,:)*mag,'length',5);
    axis equal
    arrow(sim90.mean(sel,:),sim90.mean(sel,:)+simDiff.mean(sel,:)*mag,'facecolor','r','edgecolor','r','length',5)
    delete(h);
    xlabel(char(abs('X'-1+sel(1))))
    ylabel(char(abs('X'-1+sel(2))))
    text(sim90.mean(sel(1),:),sim90.mean(sel(2),:),int2str(simId'), ...
     'fontsize',6,...
         'horizontal','center');
    title(sprintf(['2D Arrow plot of change of input points, magnified %d times ' ...
                   '(%s)'],mag,nStr));

    if doPrint
        title(gca(fig),'')
        print(fig,'-depsc',sprintf('edfdata/report/arrow2Din%s.eps',selStrs{i}))
        set(findobj(fig,'type','text'),'visible','off')
        print(fig,'-depsc',sprintf('edfdata/report/arrow2Din%s_nolabel.eps',selStrs{i}))
        set(findobj(fig,'type','text'),'visible','on')
    end
    
    fig=tagfigure(['arrowOut',selStrs{i}]);
    figure(fig)
    clf
    h=arrow(rec90.mean(sel,:),rec90.mean(sel,:)+recDiff.mean(sel,:)*mag,'length',5);
    axis equal
    arrow(rec90.mean(sel,:),rec90.mean(sel,:)+recDiff.mean(sel,:)*mag,'facecolor','b','edgecolor','b','length',5)
    delete(h);
    xlabel(char(abs('X'-1+sel(1))))
    ylabel(char(abs('X'-1+sel(2))))
    text(rec90.mean(sel(1),:),rec90.mean(sel(2),:),int2str(recId'), ...
     'fontsize',6,...
         'horizontal','center');
    title(sprintf(['2D Arrow plot of change of output points, magnified %d times ' ...
                   '(%s)'],mag,nStr));
    
    if doPrint
        title(gca(fig),'')
        print(fig,'-depsc',sprintf('edfdata/report/arrow2Dout%s.eps',selStrs{i}))
        set(findobj(fig,'type','text'),'visible','off')
        print(fig,'-depsc',sprintf('edfdata/report/arrow2Dout%s_nolabel.eps',selStrs{i}))
        set(findobj(fig,'type','text'),'visible','on')
    end
    
    fig=tagfigure(['arrowBoth',selStrs{i}]);
    figure(fig)
    clf
    h1=arrow(sim90.mean(sel,:),sim90.mean(sel,:)+simDiff.mean(sel,:)*mag,'length',5);
    h2=arrow(rec90.mean(sel,:),rec90.mean(sel,:)+recDiff.mean(sel,:)*mag,'length',5);
    axis equal
    arrow(sim90.mean(sel,:),sim90.mean(sel,:)+simDiff.mean(sel,:)*mag,'facecolor','r','edgecolor','r','length',5)
    arrow(rec90.mean(sel,:),rec90.mean(sel,:)+recDiff.mean(sel,:)*mag,'facecolor','b','edgecolor','b','length',5)
    delete([h1;h2]);
    xlabel(char(abs('X'-1+sel(1))))
    ylabel(char(abs('X'-1+sel(2))))
    text(sim90.mean(sel(1),:),sim90.mean(sel(2),:),int2str(simId'), ...
     'fontsize',6,...
         'horizontal','center');
    text(rec90.mean(sel(1),:),rec90.mean(sel(2),:),int2str(recId'), ...
     'fontsize',6,...
         'horizontal','center');
    title(sprintf(['2D Arrow plot of change of both points, magnified %d times ' ...
                   '(%s)'],mag,nStr));

    if doPrint
        title(gca(fig),'')
        print(fig,'-depsc',sprintf('edfdata/report/arrow2Dboth%s.eps',selStrs{i}))
        set(findobj(fig,'type','text'),'visible','off')
        print(fig,'-depsc',sprintf('edfdata/report/arrow2Dboth%s_nolabel.eps',selStrs{i}))
        set(findobj(fig,'type','text'),'visible','on')
    end
end

% Convert input points to polar
sim90pol=phasewrap(diag([180/pi,1,1])*cart2pol(sim90-repmat(pt0,1,size(sim90,2))),-190,360);
rec90pol=phasewrap(diag([180/pi,1,1])*cart2pol(rec90-repmat(pt0,1,size(rec90,2))),-190,360);
sim07pol=phasewrap(diag([180/pi,1,1])*cart2pol(sim07-repmat(pt0,1,size(sim07,2))),-190,360);
rec07pol=phasewrap(diag([180/pi,1,1])*cart2pol(rec07-repmat(pt0,1,size(rec07,2))),-190,360);
simDiffPol=phasewrap(sim07pol-sim90pol,-180,360);
recDiffPol=phasewrap(rec07pol-rec90pol,-180,360);

fig=tagfigure('polpos');
figure(fig)
clf
text(sim90pol.mean(1,:),sim90pol.mean(3,:),int2str(simId'), ...
     'fontsize',6,...
     'horizontal','center','color','r');
text(rec90pol.mean(1,:),rec90pol.mean(3,:),int2str(recId'), ...
     'fontsize',6,...
     'horizontal','center','color','b');
axis([-180,180,min(min(sim90pol.mean(3,:)),min(rec90pol.mean(3,:))),...
      max(max(sim90pol.mean(3,:)),max(rec90pol.mean(3,:)))]+10*[-1,1,-1,1]);
xlabel('theta (degrees)')
ylabel('Z (m)')
title('Polar plot of input coordinates')
if doPrint
    title(gca(fig),'')
    print(fig,'-depsc','edfdata/report/polarlegend.eps');
end
    
selStrs={'thZ'};
sels={[1,3]};
for i=1:length(sels)
    sel=sels{i};
    fig=tagfigure(['arrowIn',selStrs{i}]);
    figure(fig)
    clf
    h=arrow(sim90pol.mean(sel,:),sim90pol.mean(sel,:)+simDiffPol.mean(sel,:)*mag,'length',5);
    axis equal
    arrow(sim90pol.mean(sel,:),sim90pol.mean(sel,:)+simDiffPol.mean(sel,:)*mag,'facecolor','r','edgecolor','r','length',5)
    delete(h);
    xlabel('theta')
    ylabel(char(abs('X'-1+sel(2))))
    text(sim90pol.mean(sel(1),:),sim90pol.mean(sel(2),:),int2str(simId'), ...
     'fontsize',6,...
         'horizontal','center');
    title(sprintf(['2D Arrow plot of change of input points, magnified %d times ' ...
                   '(%s)'],mag,nStr));

    if doPrint
        title(gca(fig),'');
        print(fig,'-depsc',sprintf('edfdata/report/arrow2Din%s.eps',selStrs{i}))
        set(findobj(fig,'type','text'),'visible','off')
        print(fig,'-depsc',sprintf('edfdata/report/arrow2Din%s_nolabel.eps',selStrs{i}))
        set(findobj(fig,'type','text'),'visible','on')
    end
    
    fig=tagfigure(['arrowOut',selStrs{i}]);
    figure(fig)
    clf
    h=arrow(rec90pol.mean(sel,:),rec90pol.mean(sel,:)+recDiffPol.mean(sel,:)*mag,'length',5);
    axis equal
    arrow(rec90pol.mean(sel,:),rec90pol.mean(sel,:)+recDiffPol.mean(sel,:)*mag,'facecolor','b','edgecolor','b','length',5)
    delete(h);
    xlabel('theta')
    ylabel(char(abs('X'-1+sel(2))))
    text(rec90pol.mean(sel(1),:),rec90pol.mean(sel(2),:),int2str(recId'), ...
     'fontsize',6,...
         'horizontal','center');
    title(sprintf(['2D Arrow plot of change of output points, magnified %d times ' ...
                   '(%s)'],mag,nStr));
    
    if doPrint
        title(gca(fig),'');
        print(fig,'-depsc',sprintf('edfdata/report/arrow2Dout%s.eps',selStrs{i}))
        set(findobj(fig,'type','text'),'visible','off')
        print(fig,'-depsc',sprintf('edfdata/report/arrow2Dout%s_nolabel.eps',selStrs{i}))
        set(findobj(fig,'type','text'),'visible','on')
    end
    
    fig=tagfigure(['arrowBoth',selStrs{i}]);
    figure(fig)
    clf
    h1=arrow(sim90pol.mean(sel,:),sim90pol.mean(sel,:)+simDiff.mean(sel,:)*mag,'length',5);
    h2=arrow(rec90pol.mean(sel,:),rec90pol.mean(sel,:)+recDiffPol.mean(sel,:)*mag,'length',5);
    axis equal
    arrow(sim90pol.mean(sel,:),sim90pol.mean(sel,:)+simDiff.mean(sel,:)*mag,'facecolor','r','edgecolor','r','length',5)
    arrow(rec90pol.mean(sel,:),rec90pol.mean(sel,:)+recDiffPol.mean(sel,:)*mag,'facecolor','b','edgecolor','b','length',5)
    delete([h1;h2]);
    xlabel('theta')
    ylabel(char(abs('X'-1+sel(2))))
    text(sim90pol.mean(sel(1),:),sim90pol.mean(sel(2),:),int2str(simId'), ...
     'fontsize',6,...
         'horizontal','center');
    text(rec90pol.mean(sel(1),:),rec90pol.mean(sel(2),:),int2str(recId'), ...
     'fontsize',6,...
         'horizontal','center');
    title(sprintf(['2D Arrow plot of change of both points, magnified %d times ' ...
                   '(%s)'],mag,nStr));

    if doPrint
        title(gca(fig),'');
        print(fig,'-depsc',sprintf('edfdata/report/arrow2Dboth%s.eps',selStrs{i}))
        set(findobj(fig,'type','text'),'visible','off')
        print(fig,'-depsc',sprintf('edfdata/report/arrow2Dboth%s_nolabel.eps',selStrs{i}))
        set(findobj(fig,'type','text'),'visible','on')
    end
end

% Resample rms displacement and rms of variance.
th=sim90pol(1,:).mean;
rms=sqrt(mean(simDiff.mean.^2));
rmsVar=sqrt(mean(simDiff.var));
rmsScaled=sqrt(mean(simDiff.mean.^2 ./ simDiff.var));
z=sim90pol(3,:).mean;
th3s=[th,th-360,th+360];
z3s=repmat(z,1,3);
r3s=repmat(rms,1,3);
r3sv=repmat(rmsVar,1,3);
r3ss=repmat(rmsScaled,1,3);
% Resampler based on simulated points.
Fs=TriScatteredInterp(th3s',z3s',r3s','natural');
Fsv=TriScatteredInterp(th3s',z3s',r3sv','natural');
Fss=TriScatteredInterp(th3s',z3s',r3ss','natural');

th=rec90pol(1,:).mean;
rms=sqrt(mean(recDiff.mean.^2));
rmsVar=sqrt(mean(recDiff.var));
rmsScaled=sqrt(mean(recDiff.mean.^2 ./ recDiff.var));
z=rec90pol(3,:).mean;
th3r=[th,th-360,th+360];
z3r=repmat(z,1,3);
r3r=repmat(rms,1,3);
r3rv=repmat(rmsVar,1,3);
r3rs=repmat(rmsScaled,1,3);
% Resampler based on recorded points.
Fr=TriScatteredInterp(th3r',z3r',r3r','natural');
Frv=TriScatteredInterp(th3r',z3r',r3rv','natural');
Frs=TriScatteredInterp(th3r',z3r',r3rs','natural');

% Resampler based on both points.
Fb=TriScatteredInterp([th3r,th3s]',[z3r,z3s]',[r3r,r3s]','natural');
Fbv=TriScatteredInterp([th3r,th3s]',[z3r,z3s]',[r3rv,r3sv]','natural');
Fbs=TriScatteredInterp([th3r,th3s]',[z3r,z3s]',[r3rs,r3ss]','natural');

% Resampling grid.
[tt,zz]=meshgrid(-200:200,min(z)-5:max(z)+5);

fig=tagfigure('contourRMSin');
figure(fig);
rr=Fs(tt,zz);
[c,h]=contour(tt,zz,rr*1000,0:5:50);
clabel(c,h);
line(th3s,z3s,'marker','.','linestyle','none','color','r')
xlabel('theta (degrees)')
ylabel('Z (m)')
title(sprintf(['Contour plot of RMS displacement [mm] based on input ',...
               'points (%s)'],nStr));
if doPrint
    title(gca(fig),'');
    print(fig,'-depsc','edfdata/report/contourRMSIn.eps');
end

fig=tagfigure('contourRMSout');
figure(fig);
rr=Fr(tt,zz);
[c,h]=contour(tt,zz,rr*1000,0:5:50);
clabel(c,h);
line(th3r,z3r,'marker','.','linestyle','none','color','b')
xlabel('theta (degrees)')
ylabel('Z (m)')
title(sprintf(['Contour plot of RMS displacement [mm] based on output ',...
               'points (%s)'],nStr));
if doPrint
    title(gca(fig),'');
    print(fig,'-depsc','edfdata/report/contourRMSOut.eps');
end

fig=tagfigure('contourRMSboth');
figure(fig);
rr=Fb(tt,zz);
[c,h]=contour(tt,zz,rr*1000,0:5:50);
clabel(c,h);
line(th3s,z3s,'marker','.','linestyle','none','color','r')
line(th3r,z3r,'marker','.','linestyle','none','color','b')
xlabel('theta (degrees)')
ylabel('Z (m)')
title(sprintf(['Contour plot of RMS displacement [mm] based on both ',...
               'points (%s)'],nStr));
if doPrint
    title(gca(fig),'');
    print(fig,'-depsc','edfdata/report/contourRMSBoth.eps');
end

fig=tagfigure('contourVARin');
figure(fig);
rr=Fsv(tt,zz);
[c,h]=contour(tt,zz,rr*1000);
clabel(c,h);
line(th3s,z3s,'marker','.','linestyle','none','color','r')
xlabel('theta (degrees)')
ylabel('Z (m)')
title(sprintf(['Contour plot of RMS of variance [mm] based on input ',...
               'points (%s)'],nStr));
if doPrint
    title(gca(fig),'');
    print(fig,'-depsc','edfdata/report/contourVARIn.eps');
end

fig=tagfigure('contourVARout');
figure(fig);
rr=Frv(tt,zz);
[c,h]=contour(tt,zz,rr*1000);
clabel(c,h);
line(th3r,z3r,'marker','.','linestyle','none','color','b')
xlabel('theta (degrees)')
ylabel('Z (m)')
title(sprintf(['Contour plot of RMS of variance [mm] based on output ',...
               'points (%s)'],nStr));
if doPrint
    title(gca(fig),'');
    print(fig,'-depsc','edfdata/report/contourVAROut.eps');
end

fig=tagfigure('contourVARboth');
figure(fig);
rr=Fbv(tt,zz);
[c,h]=contour(tt,zz,rr*1000);
clabel(c,h);
line(th3s,z3s,'marker','.','linestyle','none','color','r')
line(th3r,z3r,'marker','.','linestyle','none','color','b')
xlabel('theta (degrees)')
ylabel('Z (m)')
title(sprintf(['Contour plot of RMS of variance [mm] based on both ',...
               'points (%s)'],nStr));
if doPrint
    title(gca(fig),'');
    print(fig,'-depsc','edfdata/report/contourVARBoth.eps');
end

fig=tagfigure('contourSCALEDin');
figure(fig);
rr=Fss(tt,zz);
[c,h]=contour(tt,zz,rr);
clabel(c,h);
line(th3s,z3s,'marker','.','linestyle','none','color','r')
xlabel('theta (degrees)')
ylabel('Z (m)')
title(sprintf(['Contour plot of scaled RMS based on input ',...
               'points (%s)'],nStr));
if doPrint
    title(gca(fig),'');
    print(fig,'-depsc','edfdata/report/contourSCALEDIn.eps');
end

fig=tagfigure('contourSCALEDout');
figure(fig);
rr=Frv(tt,zz);
[c,h]=contour(tt,zz,rr);
clabel(c,h);
line(th3r,z3r,'marker','.','linestyle','none','color','b')
xlabel('theta (degrees)')
ylabel('Z (m)')
title(sprintf(['Contour plot of scaled RMS based on output ',...
               'points (%s)'],nStr));
if doPrint
    title(gca(fig),'');
    print(fig,'-depsc','edfdata/report/contourSCALEDOut.eps');
end

fig=tagfigure('contourSCALEDboth');
figure(fig);
rr=Fbv(tt,zz);
[c,h]=contour(tt,zz,rr);
clabel(c,h);
line(th3s,z3s,'marker','.','linestyle','none','color','r')
line(th3r,z3r,'marker','.','linestyle','none','color','b')
xlabel('theta (degrees)')
ylabel('Z (m)')
title(sprintf(['Contour plot of scaled RMS based on both ',...
               'points (%s)'],nStr));
if doPrint
    title(gca(fig),'');
    print(fig,'-depsc','edfdata/report/contourSCALEDBoth.eps');
end
