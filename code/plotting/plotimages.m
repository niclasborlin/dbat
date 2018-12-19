function hh=plotimages(s,ims,varargin)
%PLOTIMAGES Plot image from a DBAT camera network.
%
%  PLOTIMAGES(S,I), where S is a struct returned by PROB2DBATSTRUCT, and
%  I is a vector of image numbers, plots each image in I in a separate
%  window with measured coordinates. Use I=='all' to plot all images.
%
%  PLOTIMAGES(S,I,E), where E is a struct returned by BUNDLE, plots a trace
%  of bundle iterations showing the reprojected object point coordinates.
%
%  PLOTIMAGES(...,'op',J), plots the object points with indices in the
%  vector J only. J defaults to 'all'.
%
%  PLOTNETWORK(...,'lines',L), where L is a cell array of vectors with
%  object point indices, connects the object points specified in L by lines.
%
%  PLOTIMAGES(...,'axes',AX), where AX is an axes or figure handle or vector
%  of handles, plots the images in AX instead of in GCA. If AX is a figure
%  handle, gca(AX) is used.
%
%  PLOTIMAGES(...,'plotx0pts',B), plots the initial reprojected object point
%  estimates separately if B is true.
%
%  PLOTIMAGES(...,'title',STR), uses STR as the title of the axes. If
%  STR contains a '%d', it is replaced by the iteration number. If STR
%  contains a second '%d', it is replaced by the total iteration count.
%
%  PLOTIMAGES(...,'pause',P), where P is numeric scalar, pauses between
%  each iteration for P seconds. If P is 'on', waits for keypress between
%  iterations instead.
%
%  PLOTIMAGES(...,'iters',J), plots only the iterations listed in the
%  vector J. The values in J are 0-based, i.e. the value of 0 means the
%  initial values in the first trace column of E. The value 'inf' may be
%  used to mean the last iteration. J defaults to all.
%
%  PLOTIMAGES(...,'ptiters',J), leaves trace points for the iterations in
%  the vector J. The default is to move the current points during the
%  iterations. Use J=nan to get the default behaviour.
%
%See also: PLOTNETWORK, PROB2DBATSTRUCT.


% Iteration trace info.
E=[];
% Lines.
L={};
% Default axes.
ax=[];
% Should we plot X0 points separately?
plotX0pts=false;
% Title string.
titleStr='';
titleStrNums=0;
% Which points to plot.
ptPlot=[];
% Pause or not.
pauseMode=[];
% Which iters to plot.
iters=[];
% Which iters to leave point marks for.
ptIters=nan;

if ischar(ims) && strcmp(ims,'all')
    ims=1:length(s.EO.name);
end
    
while ~isempty(varargin)
    if isstruct(varargin{1})
        % PLOTIMAGES(S,E)
        E=varargin{1};
        varargin(1)=[];
        if ~isfield(E,{'trace','damping'})
            error('Bad E struct.');
        end
    elseif ischar(varargin{1})
        if length(varargin)<2
            error('DBAT:plotimages:badInput','Missing argument');
        end
        arg=[varargin{1},repmat(' ',1,4)];
        switch lower(arg(1:4))
          case 'op  '
            ptPlot=varargin{2};
          case 'line' % 'lines'
            L=varargin{2};
            % L should be cell array of vectors of indices, or empty.
            if ~isempty(L) && ~iscell(L)
                error('DBAT:plotimages:badInput',...
                      'L should be cell array of vectors of OP indices');
            end
          case 'axes' % 'axes'
            ax=varargin{2};
            % AX should be an axes or figure handle.
            if ~all(ishandle(ax)) || ~all(ismember(get(ax,'type'),{'axes','figure'}))
                error('DBAT:plotimages:badInput',...
                      'AX should be an axes or figure handle'); 
            end
            isFig=strcmp(get(ax,'type'),'figure');
            ax(isFig)=cell2mat(gca(ax(isFig)));
          case 'plot'
            plotX0pts=varargin{2};
            % plotx0pts should be scalar logical.
            if ~isscalar(plotX0pts) || ~islogical(plotX0pts)
                error('DBAT:plotimages:badInput',...
                      'B should be scalar boolean'); 
            end
          case 'titl'
            titleStr=varargin{2};
            % title string should be a string
            if ~ischar(titleStr)
                error('DBAT:plotimages:badInput',...
                      'STR should be a string'); 
            end
            % How many %d does the title string have?
            titleStrNums=length(strfind(titleStr,'%d'));
          case 'iter' % 'iterations'
            iters=varargin{2};
          case 'ptit' % 'ptiterations'
            camIters=varargin{2};
          case 'paus' % 'pause'
            pauseMode=varargin{2};
          otherwise
            error('DBAT:plotimages:badInput','Bad attribute string');
        end
        % Remove processed arguments.
        varargin(1:2)=[];
    else
        error('DBAT:plotimages:badInput','Expected char parameter');
    end
end

if isempty(ax), ax=gca; end
if isempty(ptPlot), ptPlot=1:size(s.OP.val,2); end
if isempty(iters)
    if isempty(E)
        iters=0;
    else
        iters=0:size(E.trace,2)-1;
    end
end

[ixIO,ixEO,ixOP]=indvec([nnz(s.bundle.est.IO),nnz(s.bundle.est.EO),nnz(s.bundle.est.OP)]);

if ~isempty(E)
    % Number of iterations.
    nIters=size(E.trace,2)-1;
else
    nIters=0;
end

% Plot measurements.
for i=1:length(ims)
    imNo=ims(i);
    % Load image.
    imName=fullfile(s.proj.imDir,s.EO.name{imNo});
    % Try upper/lower-case version of same name too.
    if ~exist(imName,'file')
        imNames={s.EO.name{imNo},lower(s.EO.name{imNo}),upper(s.EO.name{imNo})};
        imDirs={s.proj.imDir,lower(s.proj.imDir),upper(s.proj.imDir)};
        for j=1:length(imDirs)
            for k=1:length(imNames)
                name=fullfile(imDirs{j},imNames{k});
                if exist(name,'file')
                    imName=name;
                    break
                end
            end
        end
    end 
    if exist(imName,'file')
        im=imread(imName);
    else
        im=repmat(uint8(127),[s.IO.sensor.imSize(2,i), ...
                            s.IO.sensor.imSize(1,i),3]);
    end
    if i<=length(ax)
        imshow(im,'parent',ax(i));
    else
        ax(i)=gca(tagfigure(sprintf('pts%d',imNo)));
        imshow(im,'parent',ax(i));
    end

    % Points to plot.
    vis=s.IP.vis(:,imNo) & ismember((1:size(s.IP.vis,1))',ptPlot);
    pts=s.IP.val(:,s.IP.ix(vis,imNo));
    
    % Draw control points.
    line(pts(1,s.prior.OP.isCtrl(vis)),pts(2,s.prior.OP.isCtrl(vis)),'marker','^','color','r',...
         'linestyle','none','parent',ax(i),'tag','ctrlpts');
    
    % Draw object points.
    line(pts(1,~s.prior.OP.isCtrl(vis)),pts(2,~s.prior.OP.isCtrl(vis)),'marker','.',...
         'color','b','linestyle','none','parent',ax(i),'tag','objpts');
    
    % Set title, optionally with iteration number.
    if ~isempty(titleStr)
        switch titleStrNums
        case 0
            title(ax,titleStr);
        case 1
            title(ax,sprintf(titleStr,iter));
        otherwise
            title(ax,sprintf(titleStr,iter,nIters));
        end
    end

    % Pause if requested unless after showing last iteration.
    if ~isempty(pauseMode) && iter<nIters
        if ischar(pauseMode)
            pause
        else
            pause(pauseMode);
        end
    end
    
end

if nargout>0, hh=ax; end
