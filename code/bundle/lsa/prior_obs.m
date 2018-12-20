function [f,J,JJ]=prior_obs(x,s,asStruct)
%PRIOR_OBS Prior observation residual and Jacobian for bundle adjustment.
%
%   F=PM_PREOBS(X,S) returns the residual vector F for the prior
%   observations for the camera network defined in S compared to the
%   approximate values in X. The prior observation residual
%   F=[fIO;fEO;fOP] consists of IO, EO, and OP residuals, in that
%   order.
%   
%   [F,J]=... also returns the analytical Jacobian J.
%
%   Use the call ...=PM_PREOBS(X,S,TRUE) to have the return
%   parameters as structs with fields IO, EO, OP instead.

if nargin<3, asStruct=false; end

% Update DBAT structure with current estimates in x.
s=deserialize(s,x);

if nargout>2
    fun=@(x)feval(mfilename,x,s);
    % Numerical approximation of Jacobian (for debugging only).
    JJ=jacapprox(fun,x,1e-6);
end

if nargout<2
    % ...IO residuals...
    fIO=x(s.bundle.serial.IO.dest(s.bundle.serial.IO.obs))-...
        s.prior.IO.val(s.bundle.serial.IO.src(s.bundle.serial.IO.obs));
            
    % ...EO residuals...
    fEO=x(s.bundle.serial.EO.dest(s.bundle.serial.EO.obs))-...
        s.prior.EO.val(s.bundle.serial.EO.src(s.bundle.serial.EO.obs));

    % ...OP residuals...
    fOP=x(s.bundle.serial.OP.dest(s.bundle.serial.OP.obs))-...
        s.prior.OP.val(s.bundle.serial.OP.src(s.bundle.serial.OP.obs));
    
    if asStruct
        f=struct('IO',fIO,'EO',fEO,'OP',fOP);
    else
        f=[fIO;fEO;fOP];
    end
else
    % ...IO residuals...
    fIO=x(s.bundle.serial.IO.dest(s.bundle.serial.IO.obs))-...
        s.prior.IO.val(s.bundle.serial.IO.src(s.bundle.serial.IO.obs));
    JIO=sparse(1:length(fIO), ...
               s.bundle.serial.IO.dest(s.bundle.serial.IO.obs),1, ...
               length(fIO),length(x));
            
    % ...EO residuals...
    fEO=x(s.bundle.serial.EO.dest(s.bundle.serial.EO.obs))-...
        s.prior.EO.val(s.bundle.serial.EO.src(s.bundle.serial.EO.obs));
    JEO=sparse(1:length(fEO), ...
               s.bundle.serial.EO.dest(s.bundle.serial.EO.obs),1, ...
               length(fEO),length(x));

    % ...OP residuals...
    fOP=x(s.bundle.serial.OP.dest(s.bundle.serial.OP.obs))-...
        s.prior.OP.val(s.bundle.serial.OP.src(s.bundle.serial.OP.obs));
    JOP=sparse(1:length(fOP), ...
               s.bundle.serial.OP.dest(s.bundle.serial.OP.obs),1, ...
               length(fOP),length(x));

    if asStruct
        f=struct('IO',fIO,'EO',fEO,'OP',fOP);
        J=struct('IO',JIO,'EO',JEO,'OP',JOP);
    else
        f=[fIO;fEO;fOP];
        J=[JIO;JEO;JOP];
    end
end
