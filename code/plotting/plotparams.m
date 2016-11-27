function hh=plotparams(s,e,varargin)
%PLOTPARAMS Plot bundle iteration parameters.
%
%   PLOTPARAMS(S,E), where S is a struct returned by PROB2DBATSTRUCT and E
%   is a struct returned by BUNDLE, plots the iteration trace of the
%   parameters estimated by BUNDLE.
%
%   PLOTPARAMS(...,'noop') does not plot the OP plot.
%   PLOTPARAMS(...,'noeo') does not plot the EO plot.
%   PLOTPARAMS(...,'noio') does not plot the IO plot.
%   PLOTPARAMS(...,'noparams') does not plot the params plot.
%
%   H=... also returns the figure handles for the IO/EO/OP/damping plots.
%
%See also: BUNDLE, PROB2DBATSTRUCT.


plotIO=true;
plotEO=true;
plotOP=true;
plotParams=true;

for i=1:length(varargin)
    if ~ischar(varargin{i})
        error('Bad argument');
    end
    switch lower(varargin{i})
      case 'noio'
        plotIO=false;
      case 'noeo'
        plotEO=false;
      case 'noop'
        plotOP=false;
      case 'noparams'
        plotParams=false;
      otherwise
        error('Bad string argument');
    end
end
    
h=nan(4,1);

[ixIO,ixEO,ixOP]=indvec([nnz(s.estIO),nnz(s.estEO),nnz(s.estOP)]);

if plotIO && any(s.estIO(:))
    ixPos=double(s.estIO);
    ixPos(s.estIO)=ixIO;
    
    % IO parameter plot.
    fig=tagfigure(sprintf('paramplot_io_%s',e.damping.name));
    h(1)=fig;
    clf(fig);

    % Axes in this plot.
    axH=[];
    ax=subplot(3,1,1,'parent',fig);
    axH(end+1)=ax;
    cla(ax);
    cc=get(ax,'colororder');
    % Legend strings.
    lgs={};
    % Line styles.
    ls={'-','--','-.'};
    
    % For each camera.
    for ci=1:size(s.IO,2)
        % Create array with fixed focal length, principal point.
        fp=repmat(s.IO(1:3,ci),1,size(e.trace,2));
        % Update with estimated values.
        ixp=ixPos(1:3,ci);
        fp(s.estIO(1:3,ci),:)=e.trace(ixp(s.estIO(1:3,ci)),:);
        % Flip y coordinate.
        v=diag([1,-1,1])*fp;
        % Change order to be f, px, py.
        v=v([3,1,2],:);
        % Line style and legend strings.
        ls={'-','--','-.'};
        fps={'f','px','py'};
        for i=1:size(v,1)
            if size(s.IO,2)==1
                color=cc(i,:);
                lgs{end+1}=fps{i};
            else
                color=cc(rem(ci-1,size(cc,1))+1,:);
                lgs{end+1}=sprintf('%s-%d',fps{i},ci);
            end
            line(0:size(e.trace,2)-1,v(i,:),'parent',ax,'linestyle',ls{i},...
                 'marker','x','color',color);
        end
    end
    legh=legend(lgs,'location','NorthEastOutside');
    title(ax,sprintf('Focal length, principal point (%s)',e.damping.name));
    set(ax,'xtick',0:size(e.trace,2)-1);
    set(ax,'xlim',[0,size(e.trace,2)-1]);
    
    ax=subplot(3,1,2,'parent',fig);
    axH(end+1)=ax;
    % Legend strings.
    lgs={};
    % For each camera.
    for ci=1:size(s.IO,2)
        % Create array with K1-K3 parameters.
        K=repmat(s.IO(4:6,ci),1,size(e.trace,2));
        % Update with estimated values.
        ixp=ixPos(4:6,ci);
        K(s.estIO(4:6,ci),:)=e.trace(ixp(s.estIO(4:6,ci)),:);
        % Scale K values.
        avgScale=floor(median(log10(abs(K)),2));
        % Use scaling of 1 for all-zero values.
        avgScale(avgScale==-inf)=0;
        v=repmat(10.^(-avgScale),1,size(K,2)).*K;
        for i=1:size(v,1)
            if avgScale(i)==0
                prefix='';
            else
                prefix=sprintf('10^{%d}',-avgScale(i));
            end
            if size(s.IO,2)==1
                color=cc(i,:);
                lgs{end+1}=sprintf('%sK%d',prefix,i);
            else
                color=cc(rem(ci-1,size(cc,1))+1,:);
                lgs{end+1}=sprintf('%sK%d-%d',prefix,i,ci);
            end
            
            line(0:size(e.trace,2)-1,v(i,:),'parent',ax,'linestyle',ls{i},...
                 'marker','x','color',color);
        end
        legend(lgs,'location','NorthEastOutside');
    end
    title(ax,'Radial distortion');
    set(ax,'xtick',0:size(e.trace,2)-1);
    
    ax=subplot(3,1,3,'parent',fig);
    axH(end+1)=ax;
    % Legend strings.
    lgs={};
    % For each camera.
    for ci=1:size(s.IO,2)
        % Create array with P1-P2 parameters.
        P=repmat(s.IO(7:8,ci),1,size(e.trace,2));
        % Update with estimated values.
        ixp=ixPos(7:8,ci);
        P(s.estIO(7:8,ci),:)=e.trace(ixp(s.estIO(7:8,ci)),:);
        avgScale=floor(median(log10(abs(P)),2));
        % Use scaling of 1 for all-zero values.
        avgScale(avgScale==-inf)=0;
        v=repmat(10.^(-avgScale),1,size(P,2)).*P;
        for i=1:size(P,1)
            if avgScale(i)==0
                prefix='';
            else
                prefix=sprintf('10^{%d}',-avgScale(i));
            end
            if size(s.IO,2)==1
                color=cc(i,:);
                lgs{end+1}=sprintf('%sP%d',prefix,i);
            else
                color=cc(rem(ci-1,size(cc,1))+1,:);
                lgs{end+1}=sprintf('%sP%d-%d',prefix,i,ci);
            end

            line(0:size(e.trace,2)-1,v(i,:),'parent',ax,'linestyle',ls{i},...
                 'marker','x','color',color);
        end
        legend(lgs,'location','NorthEastOutside');
    end
    title(ax,'Tangential distortion');
    set(ax,'xtick',0:size(e.trace,2)-1);
    xlabel(ax,'Iteration count')

    % Scale axes to have same width.
    scalewidth(axH);
end

if plotEO && any(s.estEO(:))
    ixPos=double(s.estEO);
    ixPos(s.estEO)=ixEO;

    % EO parameter plot.
    fig=tagfigure(sprintf('paramplot_eo_%s',e.damping.name));
    h(2)=fig;
    clf(fig);

    % Axes in this plot.
    axH=[];
    ax=subplot(2,1,1,'parent',fig);
    axH(end+1)=ax;
    cla(ax);
    cc=get(ax,'colororder');
    % Legend strings.
    lgs={};
    % Line styles.
    ls={'-','--','-.'};

    % Callback to clear all highlights in figure and highlight lines
    % corresponding to clicked line.
    cb=@highlight;
    
    % Plot each coordinate as the outer loop to get a better legend.
    % (There is a better way, but now this works.)
    for i=1:3
        % For each camera.
        for ci=1:size(s.EO,2)
            % Create array with camera centers.
            c=repmat(s.EO(1:3,ci),1,size(e.trace,2));
            % Update with estimated values.
            ixp=ixPos(1:3,ci);
            c(s.estEO(1:3,ci),:)=e.trace(ixp(s.estEO(1:3,ci)),:);
            % Line style and legend strings.
            ls={'-','--','-.'};
            fps={'X0','Y0','Z0'};
            color=cc(rem(ci-1,size(cc,1))+1,:);
            if i==1
                lgs{end+1}=sprintf('C%d',ci);
            end
            line(0:size(e.trace,2)-1,c(i,:),'parent',ax,'linestyle',ls{i},...
                 'marker','x','color',color,...
                 'tag',sprintf('%c0-%d',abs('X')-1+i,ci),...
                 'userdata',ci,'buttondownfcn',cb);
        end
        if i==1
            [legh,objh,outh,outm]=legend(lgs,'location','NorthEastOutside');
            % First comes text handles, then line handles.
            lineH=reshape(objh(size(s.EO,2)+1:end),2,[]);
            % Set lines to highlight when selected.
            set(lineH','selectionhighlight','on');
            for j=1:size(s.EO,2)
                set(lineH(:,j),'userdata',j,'buttondownfcn',cb,'hittest','on');
            end
        end
    end
    set(ax,'xtick',0:size(e.trace,2)-1);
    title(ax,sprintf('Camera center (%s)',e.damping.name));
    
    ax=subplot(2,1,2,'parent',fig);
    axH(end+1)=ax;
    % Angle strings.
    aStrs={'\omega','\phi','\kappa'};
    
    % For each camera.
    for ci=1:size(s.EO,2)
        % Create array with K1-K3 parameters.
        angles=repmat(s.EO(4:6,ci),1,size(e.trace,2));
        % Update with estimated values.
        ixp=ixPos(4:6,ci);
        angles(s.estEO(4:6,ci),:)=e.trace(ixp(s.estEO(4:6,ci)),:);
        % Convert to degrees.
        angles=angles*180/pi;
        for i=1:size(angles,1)
            color=cc(rem(ci-1,size(cc,1))+1,:);
            
            line(0:size(e.trace,2)-1,angles(i,:),'parent',ax,...
                 'linestyle',ls{i}, ... 
                 'marker','x','color',color,...
                 'tag',sprintf('%s-%d',aStrs{i},ci),'userdata',ci,...
                 'buttondownfcn',cb);

        end
    end
    title(ax,'Euler angles [degrees]');
    set(ax,'xtick',0:size(e.trace,2)-1);
    xlabel(ax,'Iteration count')

    % Scale axes to have same width.
    scalewidth(axH);
end

if plotOP && any(s.estOP(:))
    ixPos=double(s.estOP);
    ixPos(s.estOP)=ixOP;

    % OP parameter plot.
    fig=tagfigure(sprintf('paramplot_op_%s',e.damping.name));
    h(3)=fig;
    clf(fig);

    ax=subplot(1,1,1,'parent',fig);
    cla(ax);
    cc=get(ax,'colororder');
    % Legend strings.
    lgs={};
    % Line styles.
    ls={'-','--','-.'};

    % Callback to clear all highlights in figure and highlight lines
    % corresponding to clicked line.
    cb=@highlight;
    
    % Plot each coordinate as the outer loop to get a better legend.
    for i=1:3
        % For each point.
        for ci=1:size(s.OP,2)
            % Create array with OP coordinates.
            c=repmat(s.OP(1:3,ci),1,size(e.trace,2));
            % Update with estimated values.
            ixp=ixPos(1:3,ci);
            c(s.estOP(1:3,ci),:)=e.trace(ixp(s.estOP(1:3,ci)),:);
            % Line style and legend strings.
            ls={'-','--','-.'};
            fps={'X','Y','Z'};
            color=cc(rem(ci-1,size(cc,1))+1,:);
            if i==1
                lgs{end+1}=sprintf('P%d',ci);
            end
            line(0:size(e.trace,2)-1,c(i,:),'parent',ax,'linestyle',ls{i},...
                 'marker','x','color',color,...
                 'tag',sprintf('%c0-%d',abs('X')-1+i,ci),...
                 'userdata',ci,'buttondownfcn',cb);
        end
        if i==1
            [legh,objh,outh,outm]=legend(lgs,'location','NorthEastOutside');
            % First comes text handles, then line handles.
            lineH=reshape(objh(size(s.OP,2)+1:end),2,[]);
            % Set lines to highlight when selected.
            set(lineH','selectionhighlight','on');
            for j=1:size(s.OP,2)
                set(lineH(:,j),'userdata',j,'buttondownfcn',cb,'hittest','on');
            end
        end
    end
    title(ax,sprintf('Object points (%s)',e.damping.name));
    set(ax,'xtick',0:size(e.trace,2)-1);
    xlabel(ax,'Iteration count')

    % Scale axes to have same width.
    scalewidth(ax);
end    

if plotParams
    % Always generate iteration method plot.
    
    % IO parameter plot.
    fig=tagfigure(sprintf('paramplot_damping_%s',e.damping.name));
    h(4)=fig;
    clf(fig);

    switch e.damping.name
      case {'gm','gna'}
        if strcmp(e.damping.name,'gna')
            alpha=e.damping.alpha;
        else
            alpha=ones(1,length(e.res)-1);
        end
        ax=subplot(2,1,1,'parent',fig);
        semilogy(ax,0:length(e.res)-1,e.res,'x-');
        set(ax,'xtick',0:size(e.trace,2)-1);
        title(ax,sprintf('Residual norm (%s)',e.damping.name));
        
        ax=subplot(2,1,2,'parent',fig);
        plot(ax,0:length(e.res)-1,[nan,alpha],'x-');
        set(ax,'ylim',[0,1.1])
        set(ax,'xtick',0:size(e.trace,2)-1);
        xlabel(ax,'Iteration count')
        title(ax,'Step length (alpha)');
      case 'lm'
        % Replace zero lambda by non-zero to appear in log plot.
        lambda=e.damping.lambda;
        lambda(lambda==0)=e.damping.lambdaMin/100;
        ax=subplot(2,1,1,'parent',fig);
        semilogy(ax,0:length(e.res)-1,e.res,'x-');
        set(ax,'xtick',0:size(e.trace,2)-1);
        title(ax,sprintf('Residual norm (%s)',e.damping.name));
        ax1=ax;
        
        % Use plotyy to be able to set 0 label in semilogy plot.
        ax=subplot(2,1,2,'parent',fig);
        [h,l1,l2]=plotyy(ax,0:length(e.res)-1,lambda,[0,length(e.res)-1],e.damping.lambda0*[1,1],'semilogy','semilogy');
        set(l1,'marker','x');
        set(l2,'linestyle','--');
        set(h,'xtick',0:size(e.trace,2)-1);
        % Same y limits.
        set(h(2),'ylim',get(h(1),'ylim'));
        % Remove all primary yticks below lambdaMin.
        yTick=get(h(1),'ytick');
        set(h(1),'ytick',yTick(yTick>=e.damping.lambdaMin));
        % Set secondary yticks and labels.
        set(h(2),'ytick',e.damping.lambdaMin*[1/100,1]);
        set(h(2),'yticklabel',{'0','lambdaMin'})
        % Move y axis to the left and set default colors.
        set(h,'yaxislocation','left','ycolor',get(ax1,'ycolor'));
        xlabel(ax,'Iteration count')
        title(ax,'Damping (lambda)')
      case 'lmp'
        ax=subplot(4,1,1,'parent',fig);
        semilogy(ax,0:length(e.res)-1,e.res,'x-');
        set(ax,'xtick',0:size(e.trace,2)-1);
        set(ax,'xlim',[0,size(e.trace,2)]);
        title(ax,sprintf('Residual norm (%s)',e.damping.name));
        
        ax=subplot(4,1,2,'parent',fig);
        semilogy(ax,0:length(e.res)-1,e.damping.delta,'x-');
        set(ax,'xtick',0:size(e.trace,2)-1);
        set(ax,'xlim',[0,size(e.trace,2)]);
        title(ax,'Damping (delta)');
    
        ax=subplot(4,1,3,'parent',fig);
        plot(ax,0:length(e.res)-1,e.damping.step,'x-');
        set(ax,'xtick',0:size(e.trace,2)-1);
        set(ax,'xlim',[0,size(e.trace,2)]);
        set(ax,'ytick',0:2,'ylim',[-0.1,2.1],...
               'yticklabel',{'Gauss-Newton','Interpolated','Cauchy'});
        
        title(ax,'Step type');
        
        ax=subplot(4,1,4,'parent',fig);
        plot(ax,0:length(e.damping.rho)-1,e.damping.rho,'x-');
        line([0,length(e.damping.rho)-1],e.damping.rhoBad*[1,1],...
             'linestyle','--','parent',ax);
        line([0,length(e.damping.rho)-1],e.damping.rhoGood*[1,1],...
             'linestyle','--','parent',ax);
        set(ax,'xtick',0:size(e.trace,2)-1);
        set(ax,'xlim',[0,size(e.trace,2)]);
        set(ax,'ylim',[-0.1,1.1]);
        title(ax,'Gain ratio (rho)');
        xlabel(ax,'Iteration count')
    end
end

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
