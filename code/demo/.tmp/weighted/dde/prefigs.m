doPrint=true;

fig=0;

for year=[1990,2007]

[idsToSimulate,idsToRecord,ctrlId,pt,ptPMerr,PMid,allId,allPMerr]=load_cpt(year,true);
idsToRemove=[247,489,492,493,494];
idsToSimulate=setdiff(idsToSimulate,[339,340,342,469,470]);
idsToSimulate=setdiff(idsToSimulate,idsToRemove);
idsToRecord=setdiff(idsToRecord,[339,340,342,469,470]);
idsToRecord=setdiff(idsToRecord,idsToRemove);

textPar2D={'horizontal','right','vertical','top','fontsize',8};
textPar3D={'horizontal','right','vertical','top','fontsize',8};

% Present results with two magnifications
for mag=5000
    fig=fig+1;
    figure(fig)
    plot(pt(1:2,ismember(ctrlId,idsToRecord)),mag,'color',[0,0.75,0])
    hold on
    plot(pt(1:2,ismember(ctrlId,idsToSimulate)),mag,'color','r')
    axis equal
    hold off
    title(sprintf('%d planimetric error (mag=%d)',year,mag))
    
    fig=fig+1;
    figure(fig)
    plot(pt(1:2,ismember(ctrlId,idsToRecord)),mag,'color',[0,0.75,0])
    hold on
    plot(pt(1:2,ismember(ctrlId,idsToSimulate)),mag,'color','r')
    axis equal
    hold off
    title(sprintf('%d planimetric error (mag=%d)',year,mag))
    for i=find(ismember(ctrlId,idsToRecord))
        color=[0,0.75,0];
        text(pt(1,i).mean,pt(2,i).mean,int2str(ctrlId(i)),...
             textPar2D{:},'color',color);
    end
    for i=find(ismember(ctrlId,idsToSimulate))
        color='r';
        text(pt(1,i).mean,pt(2,i).mean,int2str(ctrlId(i)),...
             textPar2D{:},'color',color);
    end

    if doPrint
        print(fig,'-depsc',sprintf('edfdata/report/%d_planerr_id.eps',year));
    end
        
    fig=fig+1;
    figure(fig)
    plot(pt(:,ismember(ctrlId,idsToRecord)),mag,10,'facecolor',[0,0.75,0])
    hold on
    plot(pt(:,ismember(ctrlId,idsToSimulate)),mag,10,'facecolor','r')
    hold off
    axis vis3d
    set(gca,'projection','perspective')
    %title(sprintf('%d XYZ error (mag=%d)',year,mag))
    
    fig=fig+1;
    figure(fig)
    plot(pt(:,ismember(ctrlId,idsToRecord)),mag,10,'facecolor',[0,0.75,0])
    hold on
    plot(pt(:,ismember(ctrlId,idsToSimulate)),mag,10,'facecolor','r')
    hold off
    axis vis3d
    set(gca,'projection','perspective')
    for i=find(ismember(ctrlId,idsToRecord))
        color=[0,0.75,0];
        text(pt(1,i).mean,pt(2,i).mean,pt(3,i).mean,int2str(ctrlId(i)),...
             textPar3D{:},'color',color);
    end
    for i=find(ismember(ctrlId,idsToSimulate))
        color='r';
        text(pt(1,i).mean,pt(2,i).mean,pt(3,i).mean,int2str(ctrlId(i)),...
             textPar3D{:},'color',color);
    end
    %title(sprintf('%d XYZ error (mag=%d)',year,mag))
    if doPrint
        print(fig,'-depsc',sprintf('edfdata/report/%d_XYZerr_id.eps',year));
    end
end

end
