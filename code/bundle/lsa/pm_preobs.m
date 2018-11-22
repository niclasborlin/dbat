function [f,J,JJ]=pm_preobs(x,s)
%PM_PREOBS Prior observation residual and Jacobian for bundle adjustment.
%
%   F=PM_PREOBS(X,S) returns the residual vector F for the prior
%   observations for the camera network defined in S compared to the
%   approximate values in X. The prior observation residual
%   F=[fIO;fEO;fOP] consists of IO, EO, and OP residuals, in that
%   order.
%
%   [F,J]=... also returns the analytical Jacobian J.
%
%   [F,J,JJ]=... furthermore returns the numerical Jacobian J.

% Update DBAT structure with current estimates in x.
s=deserialize(s,x);

if nargout>2
    % Numerical approximation of Jacobian (for debugging only).
    JJ=jacapprox(mfilename,x,1e-6,{s});
end

if nargout<2
    % Only residual vector requested.
    
    fIO=s.IO.val(s.IO.prior.use)-s.IO.prior.val(s.IO.prior.use);
    fEO=s.EO.val(s.EO.prior.use)-s.EO.prior.val(s.EO.prior.use);
    fOP=s.OP.val(s.OP.prior.use)-s.OP.prior.val(s.OP.prior.use);
    
    f=[fIO(:);fEO(:);fOP(:)];
else
    fIO=s.IO.val(s.IO.prior.use)-s.IO.prior.val(s.IO.prior.use);
    JIO=sparse(1:length(fIO),ixIO(s.IO.prior.use),1,length(fIO),length(x));

    fEO=EO(s.EO.prior.use)-s.EO.prior.val(s.EO.prior.use);
    JEO=sparse(1:length(fEO),ixEO(s.EO.prior.use),1,length(fEO),length(x));

    fOP=OP(s.OP.prior.use)-s.OP.prior.val(s.OP.prior.use);
    JOP=sparse(1:length(fOP),ixOP(s.OP.prior.use),1,length(fOP),length(x));

    f=[fIO(:);fEO(:);fOP(:)];
    J=[JIO;JEO;JOP];
end


