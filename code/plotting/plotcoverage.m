function hh=plotcoverage(s,varargin)
%PLOTCOVERAGE Plot image coverage.
%
%   PLOTCOVERAGE(S), where S is a struct returned by PROB2DBATSTRUCT, plots
%   the image coverage for the project. The rectangular, convex hull, and
%   radial coverage are plotted separately. If the images were taken by
%   different cameras, the radial coverage may be uninformative.
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


ix='all';
plotEach=false;

for i=1:length(varargin)
    if isnumeric(varargin{i}) || strcmp(varargin{i},'all')
        ix=varargin{i};
    elseif islogical(varargin{i})
        plotEach=varargin{i};
    else
        error('Bad parameter');
    end
end

if strcmp(ix,'all'), ix=1:size(s.EO.val,2); end

h=tagfigure('coverage');
set(h,'name','Image coverage');

clf(h);

ax=gca(h);

cla(ax);
set(ax,'xlim',[0.5,s.IO.sensor.imSize(1,1)+0.5],...
       'ylim',[0.5,s.IO.sensor.imSize(2,1)+0.5]);
axis(ax,'image')
cc=get(ax,'colororder');
% Legend strings.
lgs={};
% Legend handles.
lh=[];
cb=@highlight;

if plotEach
    for ii=1:length(ix)
        i=ix(ii);
        % Plot observations for this image.
        j=s.IP.ix(s.IP.vis(:,i),i);
        line(s.IP.val(1,j),s.IP.val(2,j),'marker','x',...
             'linestyle','none',...
             'color',cc(rem(ii-1,size(cc,1))+1,:),...
             'tag',sprintf('obs%d',i),'userdata',i, ...
             'parent',ax,...
             'buttondownfcn',cb);
        % Plot convex hull for this image.
        [~,~,~,cHull,cl,ch,crp]=coverage(s,i);
        lh(end+1)=line(cHull{1}(1,:),cHull{1}(2,:),'linestyle','-',...
                       'color',cc(rem(ii-1,size(cc,1))+1,:),...
                       'tag',sprintf('chull%d',i),'userdata',i, ...
                       'parent',ax,...
                       'buttondownfcn',cb); %#ok<AGROW>
        lgs{end+1}=sprintf('Image %d',i); %#ok<AGROW>
        % Plot rectangular hull for this image.
        clh=[cl,ch];
        line(clh(1,[1,2,2,1,1]),clh(2,[1,1,2,2,1]),'linestyle','-.',...
             'color',cc(rem(ii-1,size(cc,1))+1,:),...
             'tag',sprintf('recthull%d',i),'userdata',i, ...
             'parent',ax,...
             'buttondownfcn',cb);
        % Plot radial hull for this image.
        % Radial line
        line(crp(1,:),crp(2,:),'linestyle','--','marker','o',...
             'color',cc(rem(ii-1,size(cc,1))+1,:),...
             'tag',sprintf('radhull%d',i),'userdata',i, ...
             'parent',ax,...
             'buttondownfcn',cb);
        % Circle
        c=crp(:,1);
        r=norm(diff(crp,[],2));
        t=linspace(0,2*pi);
        cpt=repmat(c,1,length(t))+r*[cos(t);sin(t)];
        line(cpt(1,:),cpt(2,:),'linestyle','--',...
             'color',cc(rem(ii-1,size(cc,1))+1,:),...
             'tag',sprintf('radhull%d',i),'userdata',i, ...
             'parent',ax,...
             'buttondownfcn',cb);
    end
    ii=length(ix)+1;
else
    ii=1;
    % Plot observations for all wanted images.
    j=s.IP.ix(:,ix);
    j=j(j~=0);
    % Plot all points.
    line(s.IP.val(1,j),s.IP.val(2,j),'linestyle','none','marker','x',...
         'color',cc(rem(ii-1,size(cc,1))+1,:),...
         'tag','obsall','userdata',0,...
         'parent',ax,...
         'buttondownfcn',cb);
    ii=ii+1;
end

% Coverage for all points.
[c,cr,crr,cHull,cl,ch,crp]=coverage(s,ix,true);

% Rectangular hull for all images.
clh=[cl,ch];
lh(end+1)=line(clh(1,[1,2,2,1,1]),clh(2,[1,1,2,2,1]),'linestyle','-.',...
               'color',cc(rem(ii-1,size(cc,1))+1,:),...
               'tag','rectall','userdata',-1,...
               'parent',ax,...
               'buttondownfcn',cb);
lgs{end+1}='Rect hull';
%ii=ii+1;

% Convex hull for all images.
lh(end+1)=line(cHull{1}(1,:),cHull{1}(2,:),'linestyle','-',...
               'color',cc(rem(ii-1,size(cc,1))+1,:),...
               'tag','rectall','userdata',-1,...
               'parent',ax,...
               'buttondownfcn',cb);
lgs{end+1}='Convex hull';

% Plot radial hull for this image.

% Radial line
lh(end+1)=line(crp(1,:),crp(2,:),'linestyle','--','marker','o',...
               'color',cc(rem(ii-1,size(cc,1))+1,:),...
               'tag','radall','userdata',-1,...
               'parent',ax,...
               'buttondownfcn',cb);
lgs{end+1}='Radial hull';
% Circle
c0=crp(:,1);
r=norm(diff(crp,[],2));
t=linspace(0,2*pi);
cpt=repmat(c0,1,length(t))+r*[cos(t);sin(t)];
line(cpt(1,:),cpt(2,:),'linestyle','--',...
     'color',cc(rem(ii-1,size(cc,1))+1,:),...
     'tag','radall','userdata',-1,...
     'parent',ax,...
     'buttondownfcn',cb);
ii=ii+1; %#ok<NASGU>

[legh,objh,outh,outm]=legend(ax,lh,lgs{:},'Location','NorthEastOutside'); %#ok<ASGLU>

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

if length(ix)==size(s.EO.val,2)
    title(ax,'Coverage for all images');
else
    ss=sprintf('%d,',ix);
    ss(end)=[];
    title(ax,['Coverage for images ',ss]);
end
xlabel(ax,sprintf('Rectangular/convex hull/radial coverage = %d%% / %d%% / %d%%.',...
                  round(cr*100),round(c*100),round(crr*100)));

if nargout>0, hh=h; end


% Callback function to highlight selected object(s).
function highlight(obj,event) %#ok<INUSD>

if nargin<1, obj=gcbo; end
fig=gcbf;

% Get object number.
num=get(obj,'userdata');

% All objects.
all=findobj(fig,'type','line');
% All matching objects.
sel=findobj(all,'flat','userdata',num);

% Clear any previous highlights.
set(all,'selected','off','linewidth',0.5);
% Select and set thick lines.
set(sel,'selected','on','linewidth',2);

% Move selected object to the front.

% Get all axes that are parents to selected objects.
axesWorkList=get(sel,'parent');
while ~isempty(axesWorkList)
    ax=axesWorkList{1};
    if strcmp(get(ax,'type'),'axes')
        ch=get(ax,'children');
        % Find out which of the selected objects are in this axes.
        j=ismember(ch,sel);
        [~,k]=sort(j,'descend');
        set(ax,'children',ch(k));
    end
    % Remove all references to this axes.
    eq=cellfun(@(x)isequaln(x,axesWorkList{1}),axesWorkList);
    axesWorkList=axesWorkList(~eq);
end
