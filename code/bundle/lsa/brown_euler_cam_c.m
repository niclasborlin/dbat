function [c,A,AA]=brown_euler_cam_c(x,s)
%BROWN_EULER_CAM Constraint function for the Brown camera model.
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
    % Numerical approximation of Jacobian (for debugging only).
    AA=jacapprox(mfilename,x,1e-6,{s});
end

if nargout<2
    % Only constraint vector requested.
    
    % Project into pinhole camera.
    c=pm_multirotconstr(IO,s.nK,s.nP,EO,s.cams,OP,s.vis,s.cIO,s.cEO,s.cOP);
	
else
    % Project into pinhole camera.
    [c,dIO,dEO,dOP]=pm_multirotconstr(IO,s.nK,s.nP,EO,s.cams,OP, ...
                              s.vis,s.cIO,s.cEO,s.cOP);
	
    A=[dIO,dEO,dOP];
end

