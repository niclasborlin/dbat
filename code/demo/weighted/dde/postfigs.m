year=1990;

file=sprintf('done_%d_nofw',year);

load(file);

% Calculate statistics for tie points.

% Find out how many samples are available, if from a non-complete run.
n=nnz(any(ptsToRecord~=0));

if (n<size(ptsToRecord,2))
    disp(sprintf('Detected %d samples of %d, trimming...',n,size(ptsToRecord,2)));
end

ptsToRecord=ptsToRecord(:,1:n);
ptsToSimulate=ptsToSimulate(:,1:n);

% Calculate statistics for output points.
mn=reshape(mean(ptsToRecord,2),3,[]);
C=cov(ptsToRecord');
bundlePts=PTGaussian(mn,C);

% Calculate statistics for input points for verification.
mn=reshape(mean(ptsToSimulate,2),3,[]);
C=cov(ptsToSimulate');
simulatedPts=PTGaussian(mn,C);

fig=0;

textPar2D={'horizontal','right','vertical','top','fontsize',8};
textPar3D={'horizontal','right','vertical','top','fontsize',8};

pt=11;

% Present results with two magnifications
for mag=1000
    fig=fig+1;
    figure(fig)
    ix1=1:size(bundlePts,2);
    plot(simulatedPts(1:2,:),mag,'color','r')
    hold on
    plot(bundlePts(1:2,ix1),mag,'color','b')
    axis equal
    hold off
    title(sprintf('Planimetric error (x%d)',mag))
    
    fig=fig+1;
    figure(fig)
    plot(simulatedPts(1:2,:),mag,'color','r')
    hold on
    plot(bundlePts(1:2,ix1),mag,'color','b')
    axis equal
    hold off
    for i=find(ix1(:)')
        color='b';
        text(bundlePts(1,i).mean,bundlePts(2,i).mean,int2str(idsToRecord(i)),...
             textPar2D{:},'color',color);
    end
    title(sprintf('Planimetric error (x%d)',mag))
    
    fig=fig+1;
    figure(fig)
    plot(simulatedPts,mag,10,'facecolor','r')
    hold on
    plot(bundlePts(:,ix1),mag,10,'facecolor','b')
    hold off
    axis vis3d
    %set(gca,'projection','perspective')
    title(sprintf('XYZ error (x%d)',mag))
    
    fig=fig+1;
    figure(fig)
    plot(simulatedPts,mag,10,'facecolor','r')
    hold on
    plot(bundlePts(:,ix1),mag,10,'facecolor','b')
    hold off
    axis vis3d
    %set(gca,'projection','perspective')
    for i=find(ix1(:)')
        color='b';
        text(bundlePts(1,i).mean,bundlePts(2,i).mean,bundlePts(3,i).mean,int2str(idsToRecord(i)),...
             textPar3D{:},'color',color);
    end
    title(sprintf('XYZ error (x%d)',mag))
    
    
end


