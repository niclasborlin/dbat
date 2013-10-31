function [f,J,JJ]=brown_euler_cam(x,S)
%BROWN_EULER_CAM Residual function for the Brown camera model with Euler angles.
%
%   F=BROWN_EULER_CAM(X,S) returns the residual vector F of the camera
%   network defined in S evaluated with the approximate values in X.
%
%   The internal camera model is from Brown (1971) with K1, K2, K3, P1, and
%   P2. The external orientation uses Euler omega-phi-kappa angles.
%   The values of S.cIO, S.cEO, and S.cOP indicate which parameters in
%   S.IO, S.EO, and S.OP, respectively are unknown.
%
%   References: Brown (1971), "Close-range camera calibration".
%       Photogrammetric Engineering, 37(8): 855-866.
%
%   See also PROB2DBATSTRUCT, BUNDLE, GAUSS_NEWTON_ARMIJO,
%       LEVENBERG_MARQUARDT, LEVENBERG_MARQUARDT_POWELL.

% $Id$

% Create index vectors for unknown parameters.
[ixIO,ixEO,ixOP,n]=indvec([nnz(s.cIO),nnz(s.cEO),nnz(s.cOP)]);

% Copy the current approximations of the unknown values.
IO=s.IO;
IO(s.cIO)=x(ixIO);
EO=s.EO;
EO(s.cEO)=x(ixEO);
OP=s.OP;
OP(s.cOP)=x(ixOP);

if (nargout>2)
    % Numerical approximation of Jacobian (used only for debugging).
    JJ=jacapprox(mfilename,x,1e-6,{S});
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

