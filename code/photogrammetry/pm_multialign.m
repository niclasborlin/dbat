function [EO,OP,T]=pm_multialign(EO,OP,i,ra)
%PM_MULTIALIGN Align network with specified camera.
%
%   [EO,OP]=pm_multialign(EO,OP,i[,ra])
%   EO - 7-by-M array with EO camera parameters.
%   OP - 3-by-N array with object point coordinates.
%   i  - which one of the cameras in EO that we should align to.
%   ra - roll angle in radians of camera i, as percieved by the
%        photographer, i.e. +90 degrees roll is camera up pointing to the
%        right. Defaults to 0 degrees roll angle.
%
%   After alignment, the i:th column will correspond to the origin.
%
%   [EO,OP,T]=... also returns the point transformation that was applied.


if nargin<4, ra=0; end

C=EO(1:3,i);
R=pm_eulerrotmat(EO(4:6,i));
% Use negative of roll angle since rotation axis is actually pointing
% toward the back of the camera.
RA=pm_eulerrotmat([0,0,-ra]);
T=[RA'*R*[eye(3),-C];0,0,0,1];

[EO,OP]=pm_multixform(EO,OP,T);
