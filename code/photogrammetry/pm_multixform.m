function [EO,OP,fail]=pm_multixform(EO,OP,T)
%PM_MULTIXFORM Transform a network.
%
%   [EO,OP,fail]=pm_multixform(EO,OP,T)
%   EO - 6-by-M array with EO camera parameters.
%   OP - 3-by-N array with object point coordinates.
%   T  - 4x4 homogeneous array with point transformation.


fail=false;

% Transform points.
if ~isempty(OP)
    OP=euclidean(T*homogeneous(OP));
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
    if any(isnan(P(:)))
        C=[];
    else
        C=null(P);
    end
    if size(C,2)~=1
        fail=true;
    else
        EO(1:3,i)=euclidean(C);
        % Extract angles.
        EO(4:6,i)=derotmat3d(P(:,1:3));
    end
end
