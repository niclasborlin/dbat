function [q,ui]=multiscalepts(p,IO,nK,nP,cams,nCams)
%MULTISCALESCALEPTS Scale points in many cameras from pixels to mm.
%
%   Q=SCALEPTS(P,IO,nK,nP,CAMS,NCAMS) scales the 2-by-N points P from
%   pixels to mm using the scaling constants in IO. nK and nP
%   indicates the number of lens distortion coefficients. The N-vector
%   or scalar CAMS indicate which camera was used for each point.
%
%   [Q,UI]=... also returns the scaling constants as the 2-by-NCAMS
%   array UI.

[~,~,~,~,~,u]=unpackio(IO,nK,nP);
ui=1./u;

if isscalar(cams)
    D=diag(ui(:,cams));
    q=D*p;
else
    q=nan(size(p));

    for i=1:nCams
        ix=cams==i;
        D=diag(ui(:,i));
        q(:,ix)=D*p(:,ix);
    end
end
