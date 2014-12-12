function [c,dRotP]=pm_rotconstr(rotP,rotM)
%PM_ROTCONSTR Rotational constraints.
%
%[c,dRotP]=pm_rotconstr(rotP,rotM);
%rotP  - Rotation parameters.
%rotM  - Rotation model.
%c     - Constraint vector.
%dRotP - Jacobian w.r.t. c.

if nargout==1
    % Rotate point to camera coordinate system.
    c=pm_roteuler2_c(rotP,rotM);
else
    % Rotate point to camera coordinate system.
    [c,dRotP]=pm_roteuler2_c(rotP,rotM);
end
