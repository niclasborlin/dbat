function [EO,OP]=pm_multixform(EO,OP,T)
%PM_MULTIXFORM Transform a network.
%
%   [EO,OP]=pm_multixform(EO,OP,T)
%   EO - 7-by-M array with EO camera parameters.
%   OP - 3-by-N array with object point coordinates.
%   T  - 4x4 homogenous array with point transformation.

% $Id: f4493a110d2d8be4244b6611542e8cd6104f8d13 $

% Transform points.
if ~isempty(OP)
    OP=euclidean(T*homogenous(OP));
end

% Transform cameras.
for i=1:size(EO,2)
    % Create camera matrix.
    R=pm_eulerrotmat(EO(4:6,i));
    C=EO(1:3,i);
    P=R*[eye(3),-C];
    
    % Apply point transformation to camera.
    P=(T'\P')';
    % Extract camera center.
    EO(1:3,i)=euclidean(null(P));
    % Extract angles.
    EO(4:6,i)=derotmat3d(P(:,1:3));
end
