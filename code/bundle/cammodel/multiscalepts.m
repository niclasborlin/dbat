function q=multiscalepts(p,ui,cams)
%MULTISCALESCALEPTS Scale points in many cameras from pixels to mm.
%
%   Q=SCALEPTS(P,UI,CAMS) scales the 2-by-N points P from pixels
%   to mm using the scaling constants in UI. The N-vector or scalar
%   CAMS indicate which camera was used for each point.

if isscalar(cams)
    D=diag(ui(:,cams));
    q=D*p;
else
    q=nan(size(p));

    for i=1:max(cams)
        ix=cams==i;
        D=diag(ui(:,i));
        q(:,ix)=D*p(:,ix);
    end
end
