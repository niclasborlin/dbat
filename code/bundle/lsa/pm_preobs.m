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
    
    fIO=IO(s.useIOobs)-s.prior.IO(s.useIOobs);
    fEO=EO(s.useEOobs)-s.prior.EO(s.useEOobs);
    fOP=OP(s.useOPobs)-s.prior.OP(s.useOPobs);
    
    f=[fIO(:);fEO(:);fOP(:)];
else
    fIO=IO(s.useIOobs)-s.prior.IO(s.useIOobs);
    JIO=sparse(1:length(fIO),ixIO(s.useIOobs),1,length(fIO),length(x));
    fEO=EO(s.useEOobs)-s.prior.EO(s.useEOobs);
    JEO=sparse(1:length(fEO),ixEO(s.useEOobs),1,length(fEO),length(x));
    fOP=OP(s.useOPobs)-s.prior.OP(s.useOPobs);
    JOP=sparse(1:length(fOP),ixOP(s.useOPobs),1,length(fOP),length(x));

    f=[fIO(:);fEO(:);fOP(:)];
    J=[JIO;JEO;JOP];
end


