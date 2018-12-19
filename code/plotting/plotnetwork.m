function hh=plotnetwork(s,varargin)
%PLOTNETWORK Plot DBAT camera network.
%
%  PLOTNETWORK(S), where S is a struct returned by PROB2DBATSTRUCT, plots
%  the camera network and object points in S.
%
%  PLOTNETWORK(S,E), where E is a struct returned by BUNDLE, plots a trace
%  of bundle iterations. The camera position estimates are shown as trace
%  plots of camera icons connected by lines. The object points are
%  re-plotted for every iteration.
%
%  H=PLOTNETWORK(...) returns the axes handle as H.
%
%  PLOTNETWORK(...,'trans',T0), where T0 is a 4x4 homogeneous transformation
%  matrix, applies the transformation T0 to all points and cameras before
%  plotting. Supplying T0='up' uses T0=blkdiag(1,[0,-1;1,0],1) that rotates
%  +Z from begin 'forward' to being 'up'.
%
%  PLOTNETWORK(...,'align',N), where N=I is an integer, aligns the network
%  such that camera I defines the origin and coordinate axes. If N=[I,RA],
%  camera I is assumed to have roll angle RA radians. Any transformation
%  T is applied after the alignment.
%
%  PLOTNETWORK(...,'lines',L), where L is a cell array of vectors with
%  object point indices, connects the object points specified in L by lines.
%
%  PLOTNETWORK(...,'camsize',CSZ), where CSZ is a 1-by-3 vector
%  [W,H,D], scales the camera icons to have width W, height H, and depth D.
%  CSZ defaults to [1,0.73,0.36]. If CSZ is a scalar, the default values are
%  scaled by CSZ instead.
%
%  PLOTNETWORK(...,'axes',AX), where AX is an axes or figure handle,
%  plots the network in AX instead of in GCA. If AX is a figure handle,
%  gca(AX) is used.
%
%  PLOTNETWORK(...,'plotx0pts',B), plots the initial object point
%  estimates separately if B is true.
%
%  PLOTNETWORK(...,'title',STR), uses STR as the title of the axes. If
%  STR contains a '%d', it is replaced by the iteration number. If STR
%  contains a second '%d', it is replaced by the total iteration count.
%
%  PLOTNETWORK(...,'name',STR), uses STR as the name of the figure.
%
%  PLOTNETWORK(...,'pause',P), where P is numeric scalar, pauses between
%  each iteration for P seconds. If P is 'on', waits for keypress between
%  iterations instead.
%
%  PLOTNETWORK(...,'EOplot',I), plots all cameras in the vector I. I
%  defaults to all.
%
%  PLOTNETWORK(...,'iters',I), plots only the iterations listed in the
%  vector I. The values in I are 0-based, i.e. the value of 0 means the
%  initial values in the first trace column of E. The value 'inf' may be
%  used to mean the last iteration. I defaults to all.
%
%  PLOTNETWORK(...,'camiters',I), plots separate camera icons for the
%  iterations listed in the vector I. The default is to plot a single
%  camera per station and move it during the iteration. Use I=nan to get
%  the default behaviour.
%
%See also: PLOTIMAGES, PROB2DBATSTRUCT, CAMERAICON, PM_MULTIALIGN.


% Defaults.

% Iteration trace info.
E=[];
% Default transformation.
T0=eye(4);
% Lines.
L={};
% Default camera size.
camSize=[1,0.73,0.36];
% Default axes.
ax=[];
% Should we plot X0 points separately?
plotX0pts=false;
% Axes title string.
titleStr='';
titleStrNums=0;
% Figure name string.
nameStr='';
% Which cameras to plot.
EOplot=[];
% Pause or not.
pauseMode=[];
% Which camera to align to.
align=[];
% Which iters to plot.
iters=[];
% Which iters to leave separate camera icons for.
camIters=nan;

while ~isempty(varargin)
    if isstruct(varargin{1})
        % PLOTNETWORK(S,E)
        E=varargin{1};
        varargin(1)=[];
        if ~isfield(E,{'trace','damping'})
            error('Bad E struct.');
        end
    elseif ischar(varargin{1})
        if length(varargin)<2
            error('DBAT:plotnetwork:badInput','Missing argument');
        end
        arg=[varargin{1},repmat(' ',1,4)];
        switch lower(arg(1:4))
          case 'tran' % 'trans'
            T0=varargin{2};
            if ischar(T0) && strcmpi(T0,'up')
                T0=blkdiag(1,[0,-1;1,0],1);
            end
            % T0 should be numeric 4x4
            if ~isnumeric(T0) || ~ismatrix(T0) || any(size(T0)~=4)
                error('DBAT:plotnetwork:badInput','T0 should be numeric 4x4');
            end
          case 'line' % 'lines'
            L=varargin{2};
            % L should be cell array of vectors of indices, or empty.
            if ~isempty(L) && ~iscell(L)
                error('DBAT:plotnetwork:badInput',...
                      'L should be cell array of vectors of OP indices');
            end
          case 'cams' % 'camsize'
            v=varargin{2};
            % CA should be scalar or 1-by-3.
            if ~isnumeric(v) || all([1,3]~=length(v))
                error('DBAT:plotnetwork:badInput',...
                      'CSZ should be scalar or 3-vector'); 
            end
            if isscalar(v)
                camSize=v*camSize;
            else
                camSize=v;
            end
          case 'axes' % 'axes'
            ax=varargin{2};
            % AX should be an axes or figure handle.
            if ~ishandle(ax) || ~ismember(get(ax,'type'),{'axes','figure'})
                error('DBAT:plotnetwork:badInput',...
                      'AX should be an axes or figure handle'); 
            end
            if strcmp(get(ax,'type'),'figure')
                ax=gca(ax);
            end
          case 'plot'
            plotX0pts=varargin{2};
            % plotx0pts should be scalar logical.
            if ~isscalar(plotX0pts) || ~islogical(plotX0pts)
                error('DBAT:plotnetwork:badInput',...
                      'B should be scalar boolean'); 
            end
          case 'titl'
            titleStr=varargin{2};
            % title string should be a string
            if ~ischar(titleStr)
                error('DBAT:plotnetwork:badInput',...
                      'STR should be a string'); 
            end
            % How many %d does the title string have?
            titleStrNums=length(strfind(titleStr,'%d'));
          case 'name'
            nameStr=varargin{2};
            % title string should be a string
            if ~ischar(nameStr)
                error('DBAT:plotnetwork:badInput',...
                      'STR should be a string'); 
            end
          case 'alig' % 'align'
            align=varargin{2};
            if ~isnumeric(align) || length(align)>2
                error('DBAT:plotnetwork:badInput',...
                      'P should be numeric scalar or 2-vector'); 
            end
          case 'iter' % 'iterations'
            iters=varargin{2};
          case 'cami' % 'camiterations'
            camIters=varargin{2};
          case 'paus' % 'pause'
            pauseMode=varargin{2};
          case 'eopl'
            EOplot=varargin{2};
          otherwise
            error('DBAT:plotnetwork:badInput','Bad attribute string');
        end
        % Remove processed arguments.
        varargin(1:2)=[];
    else
        error('DBAT:plotnetwork:badInput','Expected char parameter');
    end
end

if isempty(ax), ax=gca; end
if isempty(EOplot), EOplot=1:size(s.EO.val,2); end
if isempty(iters)
    if isempty(E)
        iters=0;
    else
        iters=0:size(E.trace,2)-1;
    end
end

camIters(any(isnan(camIters)))=[];

if ~isempty(E)
    camIters=unique(min(camIters,size(E.trace,2)-1));
    iters=unique(min(iters,size(E.trace,2)-1));
end

% Activate camera toolbar.
cameratoolbar(get(ax,'parent'),'show');

%[ixIO,ixEO,ixOP]=indvec([nnz(s.estIO),nnz(s.estEO),nnz(s.estOP)]);

if ~isempty(E)
    % Number of iterations.
    nIters=size(E.trace,2)-1;
else
    nIters=0;
end

% Camera centers.
camC=cell(size(EOplot));

% Store EO for all cameras that should be plotted.
EOsave=nan(size(s.EO.val,1),size(s.EO.val,2),length(iters));

x0OP=s.OP.val;

if isempty(E)
    % No iteration performed yes, we only have the current values.
    %IOtrace=s.IO.val;
    EOtrace=s.EO.val;
    OPtrace=s.OP.val;
else
    % Extract traces of all parameters.
    %IOtrace=deserialize(s,E,'all','IO');
    EOtrace=deserialize(s,E,'all','EO');
    OPtrace=deserialize(s,E,'all','OP');
end

for iter=iters
    % Extract base parameters for this iteration
    %IO=squeeze(IOtrace(:,:,iter+1));
    EO=squeeze(EOtrace(:,:,iter+1));
    OP=squeeze(OPtrace(:,:,iter+1));

    if ~isempty(align)
        if length(align)==2
            [EO,OP]=pm_multialign(EO,OP,align(1),align(2));
        else
            [EO,OP]=pm_multialign(EO,OP,align(1));
        end
    end

    % Transform points and cameras.
    [EO,OP,fail]=pm_multixform(EO,OP,T0); %#ok<ASGLU>

    EOsave(:,:,iter==iters)=EO;
    if iter==0
        x0OP=OP;
    end
    
    % Plot points.
    isCtrl=s.prior.OP.isCtrl;
    plot3(ax,OP(1,~isCtrl),OP(2,~isCtrl),OP(3,~isCtrl),'b.','tag','OP');
    hold(ax,'on');
    plot3(ax,OP(1,isCtrl),OP(2,isCtrl),OP(3,isCtrl),'r^','tag','CP');
    hold(ax,'off');

    if plotX0pts
        hold(ax,'on');
        plot3(ax,x0OP(1,~isCtrl),x0OP(2,~isCtrl),x0OP(3,~isCtrl),'.',...
              'color',0.5*ones(1,3),'tag','x0OP');
        hold(ax,'off');
    end

    axis(ax,'equal');
    colormap(ax,[1,0,0;0,1,0;0.5,0.5,1])	

    % Plot cameras for each of these iterations.
    plotCamIters=[camIters(camIters<iter),iter];
        
    for ci=plotCamIters
        % Plot cameras.
        EO=EOsave(:,:,ci==iters);
        for i=find(EOplot)
            % Get camera icon.
            [cam,camCol]=cameraicon(camSize);
            cam(:,:,3)=-cam(:,:,3);
            [m,n,p]=size(cam);
            
            % Camera center.
            CC=EO(1:3,i);

            % Don't try to plot cameras that are far away or undefined.
            if all(isfinite(CC))

                if ci==iter
                    camC{i}(:,end+1)=CC;
                end
        
                % Camera orientation.
                ang=EO(4:6,i);
                RR=pm_eulerrotmat(ang);
            
                % Apply transformation.
                T=RR*[eye(3),-CC];
                T(4,4)=1;
                cam1=reshape(applyhomoxform(inv(T),reshape(cam,m*n,p)')',m,n,p);
            
                hold(ax,'on');
                surf(ax,cam1(:,:,1),cam1(:,:,2),cam1(:,:,3),camCol,...
                     'tag',sprintf('iter%d',ci),'userdata',i);
                line(camC{i}(1,:),camC{i}(2,:),camC{i}(3,:),'marker','x',...
                     'parent',ax);
                hold(ax,'off');
            end
        end
    end
    
    % Draw 3d lines.
    for i=1:length(L)
        f=L{i};
        if (length(f)>1)
            xyz=OPT(:,f);
            if (~any(all(xyz==0)))
                % Don't draw features with uncalculated points.
                line(xyz(1,:),xyz(2,:),xyz(3,:),'parent',ax);
            end
        end
    end
    set(ax,'projection','perspective')
    axis(ax,'tight');
    view(ax,3)

    % Set title, optionally with iteration number.
    if ~isempty(titleStr)
        switch titleStrNums
        case 0
            title(ax,titleStr,'interpreter','none');
        case 1
            title(ax,sprintf(titleStr,iter),'interpreter','none');
        otherwise
            title(ax,sprintf(titleStr,iter,nIters),'interpreter','none');
        end
    end

    % Set name of figure
    if ~isempty(nameStr)
        set(get(ax,'parent'),'name',nameStr);
    end
        
    % Pause if requested unless after showing last iteration.
    if ~isempty(pauseMode) && iter<nIters
        if ischar(pauseMode) && strcmpi(pauseMode,'on')
            pause
        else
            pause(pauseMode);
        end
    end
    
end

if nargout>0, hh=ax; end
