function hh=plotnetwork(s,varargin)
%PLOTNETWORK Plot DBAT camera network.
%
%  PLOTNETWORK(S), where S is a struct returned by PROB2DBATSTRUCT, plots
%  the camera network and object points in S.
%
%  PLOTNETWORK(S,X), where X is a MPARAMS-by-NITER numeric array with
%  estimates of the MPARAMS free parameters in S, plots the network with
%  estimates during iterations 0,...,NITER-1. Each column of X is assumed to
%  hold the estimated parameters for successive iterations. The order is
%  assumed to be consistent with the cIO, cEO, and cOP fields of S, in that
%  order. The camera position estimates are shown as trace plots of camera
%  icons connected by lines. The object points are re-plotted for every
%  iteration.
%
%  H=PLOTNETWORK(...) returns the axes handle as H.
%
%  PLOTNETWORK(...,'trans',T), where T is a 4x4 homogeneous transformation
%  matrix, applies the transformation T to all points and cameras before
%  plotting.
%
%  PLOTNETWORK(...,'align',N), where N=I is an integer, aligns the network
%  such that camera I defines the origin and coordinate axes. If N=[I,RA],
%  camera I is assumed to have roll angle RA radians. Any transformation
%  T is applied after the alignment.

%  PLOTNETWORK(...,'lines',L), where L is a cell array of vectors with
%  object point indices, connects the object points specified in L by lines.
%
%  PLOTNETWORK(...,'camerasize',CSZ), where CSZ is a 1-by-3 vector
%  [W,H,D], scales the camera icons to have width W, height H, and depth D.
%  CSZ defaults to [1,0.73,0.36]. If CSZ is a scalar, the default values are
%  scaled by CSZ instead.
%
%  PLOTNETWORK(...,'axes',AX), where AX is an axes or figure handle,
%  plots the network in AX instead of in GCA. If AX is a figure handle,
%  gca(AX) is used.
%
%  PLOTNETWORK(...,'plotx0pts',B), plots the initial object point
%  estimates separately.
%
%  PLOTNETWORK(...,'title',STR), uses STR as the title of the axes. If
%  STR contains a '%d', it is replaced by the iteration number. If STR
%  contains a second '%d', it is replaced by the total iteration count.
%
%  PLOTNETWORK(...,'pause',P), where P is numeric scalar, pauses between
%  each iteration for P seconds. If P is 'on', waits for keypress between
%  iterations instead.
%
%  PLOTNETWORK(...,'EOplot',I), plots all cameras in the vector I. I
%  defaults to all.
%
%  PLOTNETWORK(...,'ixIO',ixIO), PLOTNETWORK(...,'ixEO',ixEO),
%  PLOTNETWORK(...,'ixOP',ixOP), specifies what rows of X correspond to
%  each unknown IO/EO/OP parameter. No error checking is done on the ixXX
%  arguments.
%
%See also: PROB2DBATSTRUCT, CAMERAICON, PM_MULTIALIGN.

% $Id$

% Defaults.
X=[];
T=eye(4);
L={};
camSize=[1,0.73,0.36];
ax=[];
plotX0pts=false;
titleStr='';
titleStrNums=0;
ixIO=[];
ixEO=[];
ixOP=[];
EOplot=[];
pauseMode=[];
align=[];

while ~isempty(varargin)
    if isnumeric(varargin{1})
        % PLOTNETWORKS(S,X)
        X=varargin{1};
        varargin(1)=[];
    elseif ischar(varargin{1})
        if length(varargin)<2
            error('DBAT:plotnetwork:badInput','Missing argument');
        end
        arg=[varargin{1},repmat(' ',1,3)];
        switch lower(arg(1:3))
        case 'tra' % 'trans'
            T=varargin{2};
            % T should be numeric 4x4
            if ~isnumeric(T) || ndims(T)~=2 || any(size(T)~=4)
                error('DBAT:plotnetwork:badInput','T should be numeric 4x4');
            end
        case 'lin' % 'lines'
            L=varargin{2};
            % L should be cell array of vectors of indices, or empty.
            if ~isempty(L) && ~iscell(L)
                error('DBAT:plotnetwork:badInput',...
                      'L should be cell array of vectors of OP indices');
            end
        case 'cam' % 'camerasize'
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
        case 'axe' % 'axes'
            ax=varargin{2};
            % AX should be an axes or figure handle.
            if ~ishandle(ax) || ~ismember(get(ax,'type'),{'axes','figure'})
                error('DBAT:plotnetwork:badInput',...
                      'AX should be an axes or figure handle'); 
            end
            if strcmp(get(ax,'type'),'figure')
                ax=gca(ax);
            end
        case 'plo'
            plotX0pts=varargin{2};
            % plotx0pts should be scalar logical.
            if ~isscalar(plotX0pts) || ~islogical(plotX0pts)
                error('DBAT:plotnetwork:badInput',...
                      'B should be scalar boolean'); 
            end
        case 'tit'
            titleStr=varargin{2};
            % title string should be a string
            if ~ischar(titleStr)
                error('DBAT:plotnetwork:badInput',...
                      'STR should be a string'); 
            end
            % How many %d does the title string have?
            titleStrNums=length(strfind(titleStr,'%d'));
        case 'ali' % 'align'
            align=varargin{2};
            if ~isnumeric(align) || length(align)>2
                error('DBAT:plotnetwork:badInput',...
                      'P should be numeric scalar or 2-vector'); 
            end
        case 'pau' % 'pause'
            pauseMode=varargin{2};
        case 'ixi'
            ixIO=varargin{2};
        case 'ixe'
            ixEO=varargin{2};
        case 'ixo'
            ixOP=varargin{2};
        case 'eop'
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
if isempty(ixIO), ixIO=find(s.cIO); end
if isempty(ixEO), ixEO=nnz(ixIO)+find(s.cEO); end
if isempty(ixOP), ixOP=nnz(ixIO)+nnz(ixEO)+find(s.cOP); end
if isempty(EOplot), EOplot=1:size(s.EO,2); end

iter=0;
while true
    % Extract base parameters.
    IO=s.IO;
    EO=s.EO;
    OP=s.OP;
    if ~isempty(X)
        % Replace unknown parameters with estimated value at iteration iter.
        IO(s.cIO)=X(ixIO,iter+1);
        EO(s.cEO)=X(ixEO,iter+1);
        OP(s.cOP)=X(ixOP,iter+1);
    end

    if ~isempty(align)
        if length(align)==2
            [EO,OP]=pm_multialign(EO,OP,align(1),align(2));
        else
            [EO,OP]=pm_multialign(EO,OP,align(1));
        end
    end
    
    % Transform points and cameras.
    [EO,OP]=pm_multixform(EO,OP,T);
    
    % Plot points.
    plot3(ax,OP(1,~s.isCtrl),OP(2,~s.isCtrl),OP(3,~s.isCtrl),'b.', ...
          OP(1,s.isCtrl),OP(2,s.isCtrl),OP(3,s.isCtrl),'r^');
    axis(ax,'equal');
    colormap(ax,[1,0,0;0,1,0;0.5,0.5,1])	

    % Plot cameras.
    for i=find(EOplot)
        % Get camera icon.
        [cam,camCol]=cameraicon(camSize);
        cam(:,:,3)=-cam(:,:,3);
        [m,n,p]=size(cam);
            
        % Camera center.
        CC=EO(1:3,i);
            
        % Camera orientation.
        ang=EO(4:6,i);
        RR=pm_eulerrotmat(ang);
            
        % Apply transformation.
        T=RR*[eye(3),-CC];
        T(4,4)=1;
        cam1=reshape(applyhomoxform(inv(T),reshape(cam,m*n,p)')',m,n,p);
            
        hold(ax,'on');
        surf(ax,cam1(:,:,1),cam1(:,:,2),cam1(:,:,3),camCol);
        hold(ax,'off');
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
    
    if ~isempty(titleStr)
        switch titleStrNums
        case 0
            title(ax,titleStr);
        case 1
            title(ax,sprintf(str,iter));
        otherwise
            title(ax,sprintf(str,iter,size(X,2)-1))
        end
    end

    if ~isempty(pauseMode) && iter<size(X,2)-1
        if ischar(pauseMode)
            pause
        else
            pause(pauseMode);
        end
    end
    
    iter=iter+1;
    
    if iter>=size(X,2)
        break;
    end
end
