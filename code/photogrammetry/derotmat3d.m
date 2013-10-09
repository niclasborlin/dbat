function [omega,phi,kappa]=derotmat3d(M)
%DEROTMAT3D Decompose 3d rotation matrix from world to image coordinate systems.
%
%[omega,phi,kappa]=derotmat3d(M)
%or
%angles=rotmat3d(M)
%M      - 3d rotation matrix.
%angles - [omega,phi,kappa]
%omega  - rotation about the x-axis,   -pi  <omega<=pi.
%phi    - rotation about the y'-axis,  -pi/2<phi  <=pi/2.
%kappa  - rotation about the z''-axis, -pi  <kappa<=pi.
%
%All angles are in radians.

% v1.0  99-12-07. Niclas Borlin, niclas@cs.umu.se.

phi=asin(M(3,1));
omega=atan2(-M(3,2),M(3,3));
kappa=atan2(-M(2,1),M(1,1));

if (nargout<2)
	omega=[omega,phi,kappa];
end
