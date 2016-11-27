function pm_plotmulti(IO,EO,OP,isCtrl,EOplot,cIO,cEO,cOP,X,ixIO,ixEO,ixOP,l,camSize,fig,T0,str)
%PM_PLOTMULTI Plot network iterations.
%
%  pm_plotmulti(IO,EO,OP,isCtrl,EOplot,cIO,cEO,cOP,X,ixIO,ixEO,ixOP,l,camSize,fig[,T0[,str]])
%  IO      - 7+nK+nP-by-K interior orientation [pp;f;K;P;a;u].
%  EO      - 7-by-N exterior orientation [C;ang;0]
%  OP      - 3-by-M object points.
%  isCtrl  - M-by-1 logical indicating which OP are control points.
%  EOplot  - N-by-1 logical with which cameras to plot.
%  cIO     - logical(size(IO)) indicating which IO elements that were estimated.
%  cEO     - logical(size(EO)) indicating which EO elements that were estimated.
%  cOP     - logical(size(OP)) indicating which OP elements that were estimated.
%  X       - P-by-nIter array with iterates. X==[] means plot what is in
%            IO, EO, OP.
%  ixIO    - indices of rows in X that containt IO elements.
%  ixEO    - indices of rows in X that containt EO elements.
%  ixOP    - indices of rows in X that containt OP elements.
%  l       - cell array with indices of object points to connect by lines.
%  camSize - size of camera icons.
%  fig     - figure to plot in.
%  T0      - 4x4 transformation matrix to improve visualisation. Defaults
%            to eye(4).
%  str     - Title string.


if nargin<16, T0=eye(4); end
if nargin<17, str=''; end

clf(fig);
ax=gca(fig);

for iter=1:max(1,size(X,2))
    % Extract estimated values.
    if ~isempty(X)
        IO(cIO)=X(ixIO,iter);
        EO(cEO)=X(ixEO,iter);
        OP(cOP)=X(ixOP,iter);
    end

    % Plot transformed object points.
    [EOT,OPT]=pm_multixform(EO,OP,T0);
    plot3(ax,OPT(1,~isCtrl),OPT(2,~isCtrl),OPT(3,~isCtrl),'b.',...
          OPT(1,isCtrl),OPT(2,isCtrl),OPT(3,isCtrl),'r^');
    axis(ax,'equal');
    colormap(ax,[1,0,0;0,1,0;0.5,0.5,1])	
        
    for i=find(EOplot)
        % Get camera icon.
        [cam,camCol]=cameraicon(camSize*[0.11,0.08,0.04]);
        cam(:,:,3)=-cam(:,:,3);
        [m,n,p]=size(cam);
            
        % Camera center.
        CC=EOT(1:3,i);
            
        % Camera orientation.
        ang=EOT(4:6,i);
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
    for i=1:length(l)
        f=l{i};
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
    
    title(ax,[str,sprintf(' Iteration %d of %d',iter-1,size(X,2)-1)]);
    if iter<size(X,2)
        pause
    end
end