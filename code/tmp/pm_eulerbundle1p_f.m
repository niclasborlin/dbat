function [f,J,JJ]=pm_eulerbundle1p_f(x,pt,ptCams,IO,nK,nP,EO,cams,OP,vis,PI,constr,orthoP,cIO,cEO,cOP,cPI,cSI,oIO,oEO,oOP,vIO,vEO,vOP,SI,shiftIx,ctrlPtFix)
%PM_EULERBUNDLE1P_F Residual fun for Euler cam bundle adjust w/ plane constr.
%
%[f,J,JJ]=pm_eulerbundle1p_f(x,pt,ptCams,IO,nK,nP,EO,cams,OP,vis,PI,constr,orthoP,cIO,cEO,cOP,cPI,cSI,oIO,oEO,oOP,vIO,vEO,vOP,shiftIx,ctrlPtFix)
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
%constr - sparse matrix indicating if obj point i is on plane j.
%orthoP - m x 2 matrix with indices of planes that are mutually orthogonal.
%cIO    - logical matrix indicating what elements of IO are unknown.
%cEO    - logical matrix indicating what elements of EO are unknown.
%cOP    - logical matrix indicating what elements of OP are unknown.
%cPI    - logical matrix indicating what elements of PI are unknown.
%f      - residual vector.
%J      - jacobian w.r.t. the unknowns in order [dEO,dOP,dPI]
%JJ     - numerical approximation of J.

% Calculate indices for unknown parameters.
base=1;

[ixIO,base]=pindex(nnz(cIO),base);
[ixEO,base]=pindex(nnz(cEO),base);

if (prod(size(cOP))==1)
	cOP=repmat(cOP,size(OP));
elseif (size(cOP,1)<size(OP,1))
	cOP=repmat(cOP,size(OP,1),1);
end
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
	JJ=jacapprox(mfilename,x,1e-6,{pt,ptCams,IO,nK,nP,EO,cams,OP,vis,PI,constr,orthoP,cIO,cEO,cOP,cPI,cSI,oIO,oEO,oOP,vIO,vEO,vOP,SI,shiftIx,ctrlPtFix});
end

if (nargout<2)
	% Project into pinhole camera.
	xy=pm_multieulerpinhole1(IO,nK,nP,EO,cams,OP,vis);

	% Remove lens distortion from measured points.
	ptCorr=pm_multilenscorr1(pt,IO,nK,nP,ptCams,size(IO,2));
	
	f=xy(:)-ptCorr(:);
    if (any(vIO(:)))
        f=[f;(x(ixIO(vIO~=0))-oIO(:))./sqrt(vIO(vIO~=0))];
    end
    if (any(vEO(:)))
        f=[f;(x(ixEO(vEO~=0))-oEO(:))./sqrt(vEO(vEO~=0))];
    end
    if (any(vOP(:)))
        f=[f;(x(ixOP(vOP~=0))-oOP(:))./sqrt(vOP(vOP~=0))];
    end
else
	% Project into pinhole camera.
	[xy,dIO1,dEO,dOP]=pm_multieulerpinhole1(IO,nK,nP,EO,cams,OP,vis,cIO,cEO,cOP);
	
	% Remove lens distortion from measured points.
	[ptCorr,dIO2]=pm_multilenscorr1(pt,IO,nK,nP,ptCams,size(IO,2),cIO);

	f=xy(:)-ptCorr(:);
	
    dPI=sparse(length(f),nnz(cPI));
    dSI=sparse(length(f),nnz(cSI));
	J=[dIO1-dIO2,dEO,dOP,dPI,dSI];
    if (any(vIO(:)))
        f=[f;(x(ixIO(vIO~=0))-oIO(:))./sqrt(vIO(vIO~=0))];
        J=[J;sparse(1:nnz(vIO),ixIO(vIO~=0),1./sqrt(vIO(vIO~=0)),nnz(vIO),size(J,2))];
    end
    if (any(vEO(:)))
        f=[f;(x(ixEO(vEO~=0))-oEO(:))./sqrt(vEO(vEO~=0))];
        J=[J;sparse(1:nnz(vEO),ixEO(vEO~=0),1./sqrt(vEO(vEO~=0)),nnz(vEO),size(J,2))];
    end
    if (any(vOP(:)))
        f=[f;(x(ixOP(vOP~=0))-oOP(:))./sqrt(vOP(vOP~=0))];
        J=[J;sparse(1:nnz(vOP),ixOP(vOP~=0),1./sqrt(vOP(vOP~=0)),nnz(vOP),size(J,2))];
    end
end

