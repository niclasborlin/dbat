function [f,J,JJ]=brown_euler_cam(x,s)
%BROWN_EULER_CAM Residual function for the Brown camera model with Euler angles.
%
%   F=BROWN_EULER_CAM(X,S) returns the residual vector F of the camera
%   network defined in S evaluated with the approximate values in X.
%
%   The internal camera model is from Brown (1971) with K1, K2, K3,
%   P1, and P2. The external orientation uses Euler omega-phi-kappa
%   angles.  The values of S.estIO, S.estEO, and S.estOP indicate what
%   parameters in S.IO, S.EO, and S.OP, respectively are unknown.
%
%   References: Brown (1971), "Close-range camera calibration".
%       Photogrammetric Engineering, 37(8): 855-866.
%
%   See also PROB2DBATSTRUCT, BUNDLE, GAUSS_NEWTON_ARMIJO,
%       LEVENBERG_MARQUARDT, LEVENBERG_MARQUARDT_POWELL.

% $Id$

% Create index vectors for unknown parameters.
[ixIO,ixEO,ixOP,n]=indvec([nnz(s.estIO),nnz(s.estEO),nnz(s.estOP)]);

% Copy the current approximations of the unknown values.
IO=s.IO;
IO(s.estIO)=x(ixIO);
EO=s.EO;
EO(s.estEO)=x(ixEO);
OP=s.OP;
OP(s.estOP)=x(ixOP);

if (nargout>2)
    % Numerical approximation of Jacobian (for debugging only).
    JJ=jacapprox(mfilename,x,1e-6,{s});
end

if (nargout<2)
    % Only residual vector requested.
    
    % Project into pinhole camera.
    xy=pm_multieulerpinhole1(IO,s.nK,s.nP,EO,s.cams,OP,s.vis);

    % Remove lens distortion from measured points.
    ptCorr=pm_multilenscorr1(diag([1,-1])*s.markPts,IO,s.nK,s.nP,s.ptCams, ...
                             size(IO,2));
	
    f=xy(:)-ptCorr(:);
else
    % Project into pinhole camera.
    [xy,dIO1,dEO,dOP]=pm_multieulerpinhole1(IO,s.nK,s.nP,EO,s.cams,OP, ...
                                            s.vis,s.estIO,s.estEO,s.estOP);
	
    % Remove lens distortion from measured points.
    [ptCorr,dIO2]=pm_multilenscorr1(diag([1,-1])*s.markPts,IO,s.nK,s.nP, ...
                                    s.ptCams,size(IO,2),s.estIO);

    f=xy(:)-ptCorr(:);
	
    J=[dIO1-dIO2,dEO,dOP];
end

