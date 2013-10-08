function [M,da1,da2,da3]=rotmat(angles,seq)
%ROTMAT Generate 3d rotation matrix from world to image coordinate systems.
%
%[M,da1,da2,da3]=rotmat(angles[,seq])
%angles - rotation angles.
%seq    - axis sequence
%         - 'xyz' or 'opk' (omega, phi, kappa)
%                 or 'rpy' (roll, pitch, yaw) (default)
%         - 'zxz' or 'ats' (azimuth, tilt, swing)
%M     - Complete rotation matrix.
%da1   - Partial derivative w.r.t. the first angle.
%da2   - Partial derivative w.r.t. the second angle.
%da3   - Partial derivative w.r.t. the third angle.
%
%All angles are in radians.

% v1.0  2003-06-23. Niclas Borlin, niclas@cs.umu.se.

if (nargin<2)
	seq='xyz';
end

switch (seq)
case {'xyz','opk','rpy'}
	M1=rotmat2d(1,-angles(1));
	M2=rotmat2d(2,angles(2));
	M3=rotmat2d(3,-angles(3));
	M=M3*M2*M1;
	if (nargout>1)
		Px=[0,0,0;0,0,1;0,-1,0];
		Py=[0,0,-1;0,0,0;1,0,0];
		Pz=[0,1,0;-1,0,0;0,0,0];
		da1=M*Px;
		da2=M3*M2*Py*M1;
		da3=Pz*M;
	end
case {'zxz','ats'}
	M1=rotmat2d(3,angles(1));
	M2=rotmat2d(1,-angles(2));
	M3=rotmat2d(-3,-angles(3));
	M=M3*M2*M1;
	if (nargout>1)
		Px=[0,0,0;0,0,1;0,-1,0];
		Py=[0,0,-1;0,0,0;1,0,0];
		Pz=[0,1,0;-1,0,0;0,0,0];
		da1=-M*Pz;
		da2=M3*M2*Px*M1;
		da3=Pz*M;
	end
otherwise
	error('Illegal axis sequence')
end

function M=rotmat2d(axis,phi)
%Create 2d rotation matrix.
%axis - axis to rotate about.
%phi  - angle to rotate.

%2d rotation matrix
R=[cos(phi), -sin(phi);sin(phi), cos(phi)];
M=eye(3);
switch (axis)
case 1
	M(2:3,2:3)=R;
case -1
	M(2:3,2:3)=-R;
case 2
	M([1,3],[1,3])=R;
case -2
	M([1,3],[1,3])=-R;
case 3
	M(1:2,1:2)=R;
case -3
	M(1:2,1:2)=-R;
otherwise
	error('Illegal axis');
end

