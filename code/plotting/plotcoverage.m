function hh=plotcoverage(s,varargin)
%PLOTCOVERAGE Plot image coverage.
%
%   PLOTCOVERAGE(S), where S is a struct returned by PROB2DBATSTRUCT,
%   plots the image coverage for the project.
%
%   PLOTCOVERAGE(S,IX) plots the image coverage, convex hull and
%   rectangular, for the images in IX. IX=='all' includes all images in S.
%
%   PLOTCOVERAGE(...,TRUE) additionally plots the convex hull of each
%   image separately.
%
%   H=... returns the figure handle.
%
%See also: PROB2DBATSTRUCT.

% $Id$

ix='all';
plotAll=false;

for i=1:length(varargin)
    if isnumeric(varargin{i}) || strcmp(varargin{i},'all')
        ix=varargin{i};
    elseif islogical(varargin{i})
        plotAll=varargin{i};
    else
        error('Bad parameter');
    end
end

if strcmp(ix,'all'), ix=1:size(s.EO,2); end

h=tagfigure('coverage');

clf(h);

ax=gca(h);

cla(ax);
set(ax,'xlim',[0.5,s.IO(end-3,1)+0.5],'ylim',[0.5,s.IO(end-2,1)+0.5]);
axis(ax,'image')
cc=get(ax,'colororder');
% Legend strings.
lgs={};
% Legend handles.
lh=[];
cb=@highlight;

if plotAll
    for ii=1:length(ix)
        i=ix(ii);
        % Plot observations for this image.
        j=s.colPos(s.vis(:,i),i);
        line(s.markPts(1,j),s.markPts(2,j),'marker','x',...
             'linestyle','none',...
             'color',cc(rem(ii-1,size(cc,1))+1,:),...
             'tag',sprintf('obs%d',i),'userdata',i, ...
             'parent',ax,...
             'buttondownfcn',cb);
        % Plot convex hull for this image.
        [~,~,cHull]=coverage(s,i);
        lh(end+1)=line(cHull{1}(1,:),cHull{1}(2,:),'linestyle','-',...
                       'color',cc(rem(ii-1,size(cc,1))+1,:),...
                       'tag',sprintf('chull%d',i),'userdata',i, ...
                       'parent',ax,...
                       'buttondownfcn',cb);
        lgs{end+1}=sprintf('Image %d',i);
    end
    ii=length(ix+1);
else
    ii=1;
    % Plot all points.
    line(s.markPts(1,:),s.markPts(2,:),'linestyle','none','marker','x',...
         'color',cc(rem(ii-1,size(cc,1))+1,:),...
         'tag','obsall','userdata',0,...
         'parent',ax,...
         'buttondownfcn',cb);
end

% Coverage for all points.
[c,cr,cHull,cl,ch]=coverage(s,ix,true);

% Rectangular hull for all images.
clh=[cl,ch];
lh(end+1)=line(clh(1,[1,2,2,1,1]),clh(2,[1,1,2,2,1]),'linestyle','--',...
               'color',cc(rem(ii-1,size(cc,1))+1,:),...
               'tag','rectall','userdata',-1,...
               'parent',ax,...
               'buttondownfcn',cb);
lgs{end+1}='Rect hull';
ii=ii+1;

% Conver hull for all images.
lh(end+1)=line(cHull{1}(1,:),cHull{1}(2,:),'linestyle','-',...
               'color',cc(rem(ii-1,size(cc,1))+1,:),...
               'tag','rectall','userdata',-2,...
               'parent',ax,...
               'buttondownfcn',cb);
lgs{end+1}='Convex hull';
ii=ii+1;

[legh,objh,outh,outm]=legend(ax,lh,lgs{:},'Location','NorthEastOutside');

set(objh,'selectionhighlight','off');
tags=get(objh,'tag');
for i=1:length(lh)
    % Set callback on legend lines, too.
    j=find(strcmp(lgs{i},tags));
    if ~isempty(j)
        set(objh(j),'userdata',get(lh(i),'userdata'),'buttondownfcn',cb,...
                    'hittest','on');
    end
end

if length(ix)==size(s.EO,2)
    title(ax,'Image coverage for all images');
else
    ss=sprintf('%d,',ix);
    ss(end)=[];
    title(ax,['Image coverage for images ',ss]);
end
xlabel(ax,sprintf('Rectangular coverage=%d%%. Convex hull coverage=%d%%.',...
                  round(c*100),round(cr*100)));

if nargout>0, hh=h; end


% Callback function to highlight selected object(s).
function highlight(obj,event)

if nargin<1, obj=gcbo; end
fig=gcbf;

% Get object number.
num=get(obj,'userdata');

% All objects.
all=findobj(fig,'type','line');
% All matching objects.
sel=findobj(all,'flat','userdata',num);

% Clear any previous highlights.
set(all,'linewidth',0.5);
% Select and set thick lines.
set(sel,'linewidth',2);

% Move selected object to the front.

% Get all axes that are parents to selected objects.
ax=unique(cell2mat(get(sel,'parent')));
for i=1:length(ax)
    if strcmp(get(ax(i),'type'),'axes')
        ch=get(ax(i),'children');
        % Find out which of the selected objects are in this axes.
        j=ismember(ch,sel);
        [~,k]=sort(j,'descend');
        set(ax(i),'children',ch(k));
    end
end
