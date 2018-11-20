function [hh,COP]=plotopstats(s,e,varargin)
%PLOTOPSTATS Plot object point statistics after a bundle.
%
%   PLOTOPSTATS(S,E), where S and E are struct given to and returned by
%   BUNDLE, respectively, produces a plot of object point statistics. The
%   plotted statistics are
%   - mark point count,
%   - RMS residual,
%   - X/Y/Z/total variance.
%
%   PLOTOPSTATS(...,COP) supplies a pre-computed OP covariance matrix COP.
%
%   PLOTOPSTATS(...,ID) plots the stats for the OPs with ids in ID only.
%   ID='all' plots for all OPs.
%
%   H=... returns the handle to the plot figure.
%
%   [H,COP]=... also returns the computed OP covariance matrix COP.
%
%See also: BUNDLE, COVERAGE, PLOTIMAGESTATS.


id='all';
COP=[];

for i=1:length(varargin)
    if isvector(varargin{i})
        id=varargin{i};
    else
        COP=varargin{i};
    end
end

if strcmp(id,'all'), id=s.OP.id; end

ix=find(ismember(s.OP.id,id));
% Any sorting would go here.
i=1:length(ix); % Identity sort.

fig=tagfigure('opstats');
set(fig,'name','OP stats');

axH=zeros(4,1);

clf(fig);

nPts=length(id);

% Point count plot.
ax=subplot(4,1,1,'parent',fig);
n=sum(s.IP.vis(ix,:),2);
bar(ax,1:nPts,n);
title(ax,'Point count')
set(ax,'xtick',[],'xgrid','on')
set(ax,'xlim',[1,nPts],'ylim',[0,max(n)+1]);
set(ax,'xticklabel',[])
axH(1)=ax;


% Angles..
ax=subplot(4,1,2,'parent',fig);
ang=angles(s,'Computing angles')*180/pi;
bar(ax,ix,ang(ix))
title(ax,'Maximum intersection angles')
set(ax,'xtick',[],'xgrid','on')
set(ax,'xlim',[1,nPts]);
set(ax,'ylim',[0,91])
set(ax,'xticklabel',[])
axH(2)=ax;


% Residuals.
ax=subplot(4,1,3,'parent',fig);
[rms,res]=bundle_residuals(s,e);
% Compute averages per OP.
nOP=full(sum(s.IP.vis,2));
sqSumOP=full(sum(res.^2,2));
meanOP=sqrt(sqSumOP./nOP);
bar(ax,ix,meanOP(ix));
title(ax,'RMS point residuals (-- = global RMS)')
line(get(ax,'xlim'),rms*[1,1],'marker','none','linestyle','--',...
     'color',0.5*[0,1,0],'parent',ax);
set(ax,'xtick',[],'xgrid','on')
set(ax,'xlim',[1,nPts]);
set(ax,'xticklabel',[])
axH(3)=ax;


% Variance
if isempty(COP)
    COP=bundle_cov(s,e,'COP');
end
var=reshape(diag(COP),3,[]);

ax=subplot(4,1,4,'parent',fig);
bar(ax,ix,sqrt([var(1:3,ix);sum(var(1:3,ix),1)]'));
ylabel(ax,'Project units')
legend(ax,'X','Y','Z','Total','location','NorthEastOutside')
title(ax,'OP standard deviations')
set(ax,'xtick',[],'xgrid','on')
set(ax,'xlim',[1,nPts])
set(ax,'xticklabel',' ')
axH(4)=ax;

xlabel(ax,'OP id')

% Scale plots to have equal width.
scalewidth(axH);

s=struct('ax',axH,'axL',ax,'id',id);

setticks(s);

% Keep xlim equal during zoom in all plots.
h=zoom(fig);
set(h,'actionpostcallback',@postzoom);
set(h,'motion','horizontal');
set(fig,'userdata',s);

if nargout>0, hh=fig; end

% Set tick marks such that at most 1 tick per 10 character.
function setticks(s)

% Axes handles.
ax=s.ax;
% Axes handle with id labels.
axL=s.axL;
% Number of points
nPts=length(s.id);
xLim=get(ax(1),'xlim');

% Get width in characters.
oldUnits=get(ax(1),'units');
set(ax(1),'units','characters');
pos=get(ax(1),'pos');
set(ax(1),'units',oldUnits);
w=pos(3);

% Where to put tick marks.
ix=1:nPts;
ticks=false(1,nPts);
% Where IDs jump.
jump=find(diff(s.id)~=1);
%ticks(jump)=true;
ticks(jump+1)=true;
% Which ticks are in range?
ticks(ix<xLim(1))=false;
ticks(ix>xLim(2))=false;

if nnz(ticks)>w/10
    % Too many, keep only w/10.
    ticks=find(ticks);
    ticks=ticks(round(linspace(1,end,w/10)));
else
    % Too few, let's add more.
    nToAdd=floor(w/10)-nnz(ticks);
    ticks(round(linspace(xLim(1),xLim(2),nToAdd)))=true;
    ticks=find(ticks);
end

set(ax,'xtick',ticks);

tag='oplabels';
delete(findobj(axL,'type','text','tag',tag));

pos=ticks;
ylim=get(axL,'ylim');
for i=1:length(pos)
    text(pos(i),ylim(1)-0.05*diff(ylim),int2str(s.id(ticks(i))),'parent',axL,...
         'rotation',90,'interpreter','none','horizontal','right',...
         'vertical','middle','tag',tag);
end

% Callback function after zooming.
function postzoom(obj,evd)

% Get limit of current axes.
xLim=get(evd.Axes,'xlim');
% Set limit of all axes in figure.
fig=get(evd.Axes,'parent');
s=get(fig,'userdata');
set(s.ax,'xlim',xLim);
setticks(s);
