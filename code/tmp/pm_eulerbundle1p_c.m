function [c,A,AA]=pm_eulerbundle1p_c(x,pt,ptCams,IO,nK,nP,EO,cams,OP,vis,PI,constr,orthoP,cIO,cEO,cOP,cPI,cSI,oIO,oEO,oOP,vIO,vEO,vOP,SI,shiftIx,ctrlPtFix)
%PM_EULERBUNDLE1P_C Constraint fun for Euler cam bundle adjust w/ plane constr.
%
%[f,J,JJ]=pm_eulerbundle1p_c(x,pt,ptCams,IO,nK,nP,EO,cams,OP,vis,PI,constr,orthoP,cIO,cEO,cOP,cPI,cSI,oIO,oEO,oOP,vIO,vEO,vOP,shiftIx,ctrlPtFix)
%x      - vector of unknowns.
%pt     - measured projected points.
%ptCams - physical camera for each point.
%IO     - vector with camera inner orientation as columns
%         pp - principal point [xp;yp].
%         f  - focal length.
%         K  - radial lens distortion coefficients.
%         P  - tangential lens distortion coefficients.
%         a  - with affine lens distortion coefficients.
%         u  - image unit.
%nK     - number of radial koefficients.
%nP     - number of tangential koefficients.
%EO     - matrix with photo camera outer orientation as columns
%         C     - camera center [Xc,Yc,Zc].
%         ang   - Euler angles.
%         ax    - Euler axis sequence
%                 0 - 'xyz'.
%                 1 - 'zxz'.
%cams   - vector of camera numbers for each photo.
%OP     - 3 x nObj matrix with object coordinates
%vis    - sparse matrix indicating if obj point i is visible in photo j.
%PI     - 4 x nPI with plane coefficients.
%orthoP - m x 2 matrix with indices of planes that are mutually orthogonal.
%constr - sparse matrix indicating if obj point i is on plane j.
%cIO    - logical matrix indicating what elements of IO are unknown.
%cEO    - logical matrix indicating what elements of EO are unknown.
%cOP    - logical matrix indicating what elements of OP are unknown.
%cPI    - logical matrix indicating what elements of PI are unknown.
%c      - constraint function.
%A      - jacobian w.r.t. the unknowns in order [dEO,dOP,dPI]
%AA     - numerical approximation of A.

% Calculate indices for unknown parameters.
base=1;

[ixIO,base]=pindex(nnz(cIO),base);
[ixEO,base]=pindex(nnz(cEO),base);

if (prod(size(cOP))==1)
	cOP=repmat(cOP,size(OP));
elseif (size(cOP,1)<size(OP,1))
	cOP=repmat(cOP,size(OP,1),1);
end
ixOPbase=base;
[ixOP,base]=pindex(nnz(cOP),base);
[ixPI,base]=pindex(nnz(cPI),base);
[ixSI,base]=pindex(nnz(cSI),base);

if (size(cEO,1)<size(EO,1))
	cEO(end+1,1)=0;
end

% Copy the current approximations of the unknown values.
IO(cIO)=x(ixIO);
EO(cEO)=x(ixEO);
OP(cOP)=x(ixOP);
PI(cPI)=x(ixPI);
SI(cSI)=x(ixSI);

if (nargout>2)
	AA=jacapprox(mfilename,x,1e-6,{pt,ptCams,IO,nK,nP,EO,cams,OP,vis,PI,constr,orthoP,cIO,cEO,cOP,cPI,cSI,oIO,oEO,oOP,vIO,vEO,vOP,SI,shiftIx,ctrlPtFix});
end

% # of points in each plane.
nPts=sum(constr,2);
% One constraint for each point + one for each plane.
c=zeros(size(PI,2)+sum(nPts),1);
% Starting index-1 for each plane.
ix0=cumsum([0;nPts]);
for i=1:size(PI,2)
    % For each plane constraint...
    c(ix0(i)+[1:nPts(i)])=PI(:,i)'*[OP(:,constr(i,:));ones(1,nPts(i))];
end
% Starting index-1 for orthogonality constraint.
ixp=length(c);
% Plane vector constraint.
c(ix0(end)+1:end)=sum(PI.^2)-1;

% Orthogonality constraints.
for i=1:size(orthoP)
    % Plane normals.
    ni=PI(1:3,orthoP(i,1));
    nj=PI(1:3,orthoP(i,2));
    c(ixp+i)=ni'*nj;
end

if (~isempty(cSI))
    shiftFirst=length(c)+1;
    for i=1:size(shiftIx,1)
        cc=OP(:,shiftIx(i,:))-repmat(SI(:,i),1,size(ctrlPtFix{i},2))-ctrlPtFix{i};
        c=[c;cc(:)];
    end
end

if (nargout>1)
    % Column index in jacobian for each free object point coordinate.
    iOP=double(cOP);
    iOP(cOP)=ixOP;
    % Column index in jacobian for each plane coefficient.
    iPI=double(cPI);
    iPI(cPI)=ixPI;
    % We will need the homogenous points.
    OP1=[OP;ones(1,size(OP,2))];

    % Row, column, value vectors.
    rr=[];
    cc=[];
    vv=[];
    for i=1:size(PI,2)
        % dc/dOP
        
        % Which object point coordinates are not fixed?
        pFree=cOP(:,constr(i,:));
        % Corresponding columns of jacobian.
        cFree=iOP(:,constr(i,:));
        cFree=cFree(pFree);
        % Get constraint numbers for each point coordinate...
        % ...within plane block...
        [dummy,cNo]=find(pFree);
        % ...absolute constraint number.
        cNo=cNo+ix0(i);
        % Corresponding plane coefficients
        pc=repmat(PI(1:3,i),1,nPts(i));
        % Add row,column,value triplets.
        rr=[rr;cNo];
        cc=[cc;cFree];
        vv=[vv;pc(pFree)];
        
        % dc/dPI
        dPIr=repmat(ix0(i)+[1:nPts(i)],nnz(cPI(:,i)),1);
        dPIc=repmat(iPI(cPI(:,i),i),1,nPts(i));
        dPIv=OP1(cPI(:,i),constr(i,:));
        rr=[rr;dPIr(:)];
        cc=[cc;dPIc(:)];
        vv=[vv;dPIv(:)];
    end
    for i=1:size(PI,2)
        rr=[rr;repmat(ix0(end)+i,nnz(cPI(:,i)),1)];
        cc=[cc;iPI(cPI(:,i),i)];
        vv=[vv;2*PI(cPI(:,i),i)];
    end
    for i=1:size(orthoP,1)
        ii=min(orthoP(i,:));
        jj=max(orthoP(i,:));
        rr=[rr;repmat(ixp+i,nnz(cPI(1:3,[ii,jj])),1)];
        cc=[cc;iPI(cPI(1:3,ii),ii);iPI(cPI(1:3,jj),jj)];
        ni=PI(1:3,ii);
        nj=PI(1:3,jj);
        vv=[vv;nj(cPI(1:3,ii));ni(cPI(1:3,jj))];
    end
    if (nnz(cSI)>0)
        ixOP3=zeros(size(OP));
        ixOP3(:)=ixOPbase-1+cumsum(cOP(:));
        ixSI3=reshape(ixSI,3,[]);
        for i=1:size(shiftIx,1)
            n=size(ctrlPtFix{i},2);
            rr=[rr;repmat(shiftFirst-1+(1:3*n)',2,1)];
            cc=[cc;reshape(ixOP3(:,shiftIx(i,:)),[],1)];
            cc=[cc;repmat(ixSI3(:,i),n,1)];
            vv=[vv;ones(3*n,1);-ones(3*n,1)];
            shiftFirst=shiftFirst+3*n;
        end
    end
    % Preallocate jacobian.
    A=sparse(rr,cc,vv,length(c),length(x));
end
