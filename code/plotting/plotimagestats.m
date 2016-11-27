function hh=plotimagestats(s,e,ix)
%PLOTIMAGESTATS Plot image statistics after a bundle.
%
%   PLOTIMAGESTATS(S,E), where S and E are struct given to and returned
%   by BUNDLE, respectively, produces a plot of image statistics. The
%   plotted statistics are
%   - image coverage (rectangular and convex hull),
%   - point count,
%   - RMS residual,
%   - camera station X/Y/Z/total variance.
%   - camera station omega/phi/kappa/total variance.
%
%   PLOTIMAGESTATS(S,E,IX) plots the stats for the images in IX only.
%   IX='all' plots for all images.
%
%   H=... returns the handle to the plot figure.
%
%See also: BUNDLE, COVERAGE, PLOTOPSTATS.


if nargin<3, ix='all'; end

if strcmp(ix,'all'), ix=1:size(s.EO,2); end

h=tagfigure('imagestats');

clf(h);

% Axes in this plot.
axH=[];

% Coverage plots.
ax=subplot(5,1,1,'parent',h);
axH(end+1)=ax;
[c,cr,crr]=coverage(s,ix);
bar(ax,ix,[cr;c;crr]'*100,'grouped');
legend(ax,'Rectangular','Convex','Radial','Location','NorthEastOutside')
title(ax,'Image coverage (percent)')
set(ax,'xlim',[0.5,max(ix)+0.5],'ylim',[0,100]);
set(ax,'xticklabel','')

% Point count plot.
ax=subplot(5,1,2,'parent',h);
axH(end+1)=ax;
n=sum(s.vis(:,ix),1);
bar(ax,ix,n);
set(ax,'xlim',[0.5,max(ix)+0.5],'ylim',[0,max(n)*1.05]);
title(ax,'Point count')
set(ax,'xticklabel','')

% Residuals.
ax=subplot(5,1,3,'parent',h);
axH(end+1)=ax;
[rms,res]=bundle_residuals(s,e);
% Compute averages per photo.
nPhoto=full(sum(s.vis,1));
sqSumPhoto=full(sum(res.^2,1));
meanPhoto=sqrt(sqSumPhoto./nPhoto);
bar(ax,ix,meanPhoto(ix));
set(ax,'xlim',[0.5,max(ix)+0.5]);
title(ax,'RMS point residuals (-- = global RMS)')
line(get(ax,'xlim'),rms*[1,1],'marker','none','linestyle','--',...
     'color',0.5*[0,1,0],'parent',ax);
set(ax,'xticklabel','')

% Variance
CEO=bundle_cov(s,e,'CEO');
var=reshape(diag(CEO),6,[]);

ax=subplot(5,1,4,'parent',h);
axH(end+1)=ax;
bar(ax,ix,sqrt([var(1:3,ix);sum(var(1:3,ix),1)]'));
ylabel(ax,'Project units')
legend(ax,'X','Y','Z','Total','Location','NorthEastOutside')
set(ax,'xlim',[0.5,max(ix)+0.5]);
title(ax,'Spatial standard deviations (camera station)')
set(ax,'xticklabel','')

ax=subplot(5,1,5,'parent',h);
axH(end+1)=ax;
bar(ax,ix,sqrt([var(4:6,ix);sum(var(4:6,ix),1)]'*180/pi));
ylabel(ax,'Degrees')
legend(ax,'omega','phi','kappa','Total','Location','NorthEastOutside')
title(ax,'Rotational standard deviations (camera station)')
set(ax,'xlim',[0.5,max(ix)+0.5]);
xlabel(ax,'Image number')

% Scale axes to have the same width.
scalewidth(axH);

if nargout>0, hh=h; end
