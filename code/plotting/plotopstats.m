function hh=plotopstats(s,e,id)
%PLOTOPSTATS Plot object point statistics after a bundle.
%
%   PLOTOPSTATS(S,E), where S and E are struct given to and returned by
%   BUNDLE, respectively, produces a plot of object point statistics. The
%   plotted statistics are
%   - mark point count,
%   - RMS residual,
%   - X/Y/Z/total variance.
%
%   PLOTOPSTATS(S,E,ID) plots the stats for the OPs with ids in ID only.
%   ID='all' plots for all OPs.
%
%   H=... returns the handle to the plot figure.
%
%See also: BUNDLE, COVERAGE, PLOTIMAGESTATS.

% $Id$

if nargin<3, id='all'; end

if strcmp(id,'all'), id=s.OPid; end

ix=find(ismember(s.OPid,id));
% Any sorting would go here.
i=1:length(ix); % Identity sort.

fig=tagfigure('opstats');

axH=zeros(3,1);

clf(fig);

nPts=length(id);

% Where to put tick marks.
% Repeated.
ticks=rem(0:nPts-1,floor(sqrt(nPts)))==0;
% Where IDs jump.
jump=find(diff(id)~=1);
%ticks(jump)=true;
ticks(jump+1)=true;
ticks=find(ticks);

% Point count plot.
ax=subplot(3,1,1,'parent',fig);
n=sum(s.vis(ix,:),2);
bar(ax,1:nPts,n);
title(ax,'Point count')
set(ax,'xtick',ticks,'xgrid','on')
set(ax,'xlim',[1,nPts],'ylim',[0,max(n)+1]);
set(ax,'xticklabel',[])
axH(1)=ax;

% Residuals.
ax=subplot(3,1,2,'parent',fig);
[rms,res]=bundle_residuals(s,e);
% Compute averages per OP.
nOP=full(sum(s.vis,2));
sqSumOP=full(sum(res.^2,2));
meanOP=sqrt(sqSumOP./nOP);
bar(ax,ix,meanOP(ix));
title(ax,'RMS point residuals (-- = global RMS)')
line(get(ax,'xlim'),rms*[1,1],'marker','none','linestyle','--',...
     'color',0.5*[0,1,0],'parent',ax);
set(ax,'xtick',ticks,'xgrid','on')
set(ax,'xlim',[1,nPts]);
set(ax,'xticklabel',[])
axH(2)=ax;


% Variance
COP=bundle_cov(s,e,'COP');
var=reshape(diag(COP),3,[]);

ax=subplot(3,1,3,'parent',fig);
bar(ax,ix,sqrt([var(1:3,ix);sum(var(1:3,ix),1)]'));
ylabel(ax,'Project units')
legend(ax,'X','Y','Z','Total','location','NorthEastOutside')
title(ax,'OP standard deviations')
set(ax,'xtick',ticks,'xgrid','on')
set(ax,'xlim',[1,nPts])
set(ax,'xticklabel',' ')
axH(3)=ax;

pos=ticks;
labelIds=id(ticks);
ylim=get(ax,'ylim');
for i=1:length(pos)
    text(pos(i),ylim(1)-0.05*diff(ylim),int2str(id(ticks(i))),'parent',ax,...
         'rotation',90,'interpreter','none','horizontal','right',...
         'vertical','middle');
end

xlabel(ax,'OP id')

yLim=nan(length(ax),2);
for i=1:length(ax)
    yLim(i,:)=get(ax(i),'ylim');
end

% Scale plots to have equal width.
scalewidth(axH);

% Keep xlim equal during zoom in all plots.
h=zoom(fig);
set(h,'actionpostcallback',@postzoom);
set(h,'motion','horizontal');
set(fig,'userdata',struct('ax',axH));

% Callback function after zooming.
function postzoom(obj,evd)

% Get limit of current axes.
xLim=get(evd.Axes,'xlim');
% Set limit of all axes in figure.
fig=get(evd.Axes,'parent');
s=get(fig,'userdata');
set(s.ax,'xlim',xLim);

