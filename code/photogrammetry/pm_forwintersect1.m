function [OP,OPc]=pm_forwintersect1(P,xy,method)
%PM_FORWINTERSECT1 Estimate world coords of one pt by forward intersection.
%
%[OP,OPc]=pm_forwintersect1(P,xy[,'euclidean'])
%or
%[OP,A]=pm_forwintersect1(P,xy,'homogeneous')
%P      - 3nx4 matrix with 3x4 camera matrices [P1;P2;...].
%xy     - 2xn or 3xn matrix with euclidean or homogeneous point coordinates.
%method - 'homogeneous' or 'euclidean' (default).
%OP     - 3-vector (euclidean input) or 4-vector (homogeneous input) with
%         estimated point. 
%OPc    - 3xn or 4xn matrix with closest point on each line (method=euclidean).
%A      - homogeneous matrix A such that A*OP approx 0 (method=homogeneous).


if nargin<3
    method='euclidean';
end

switch (lower(method(1)))
  case 'h'
    % Homogeneous.
    A=zeros(size(xy,2)*2,4);
    for i=1:size(xy,2)
        Pi=P((i-1)*3+[1:3],:);
        xyi=xy(:,i);
        if (length(xyi)==2)
            xyi=[xyi;1];
        end
        SP=skew(xyi)*Pi;
        % Add the two rows that have the largest elements.
        [dummy,j]=sort(-sum(abs(SP')));
        A((i-1)*2+[1:2],:)=SP(j(1:2),:);
    end
    [U,S,V]=svd(A);
    OP=normhomo(V(:,end));
    if (size(xy,1)<3)
        OP=OP(1:3);
    end
    if (nargout>1)
        OPc=A;
    end
  case 'e'
    A=zeros(0,3);
    b=zeros(0,0);
    for i=1:size(xy,2)
        Pi=P((i-1)*3+[1:3],:);
        xyi=xy(:,i);
        if (length(xyi)==2)
            xyi=[xyi;1];
        end
        % Get two points in 3D on line.
        E=speye(3,4);
        C=E*normhomo(null(Pi));
        Ppx=pinv(Pi)*xyi;
        if (abs(Ppx(end))<1e-8)
            % Almost at infinity. Pick a closer point.
            Ppx=Ppx+[C;1];
        end
        t=E*normhomo(Ppx)-C;
        A(end+[1:3],1:3)=speye(3);
        A(end-2:end,end+1)=t/norm(t);
        b=[b;C];
    end
    % Solve for optimal point.
    x=A\b;
    % Optimal point.
    OP=x(1:3);
    if (size(xy,1)>2)
        OP=[OP;1];
    end
    if (nargout>1)
        % Points on lines.
        OPc=b-A(:,4:end)*x(4:end);
        OPc=reshape(OPc,3,length(OPc)/3);
        if (size(xy,1)>2)
            OPc=[OPc;ones(1,size(OPc,2))];
        end
    end
  otherwise
    error('Illegal method');
end

function SS=skew(x)
%Create cross-product matrix S such that S*y = x cross y.

SS=zeros(3);
SS([6,7,2,8,3,4])=[x;-x];
