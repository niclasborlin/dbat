function [f,J,JJ]=multi_res(IO,EO,OP,s,fun)
%MULTI_RES Compute residuals for multiple cameras.
%
%   MULTI_RES(IO,EO,OP,S,FUN) computes residuals for all image
%   observations. The current IO, EO, OP values are stored in the
%   respective arrays. The structure S contains the project
%   information. FUN is a function handle that will compute the image
%   residuals for one camera.
%
%   [F,J]=... also returns the requested Jacobian J.
%
%   See also: RES_EULER_BROWN_0, RES_EULER_BROWN_1, RES_EULER_BROWN_2,
%   RES_EULER_BROWN_3.

nPhotos=size(EO,2);
nObjs=size(OP,2);

% Which camera was used for which image?
cams=s.cams;
if length(cams)~=nPhotos
    error('%s: bad size',mfilename);
end

% Total number of projected points.
nProj=nnz(s.vis);

if nargout<2
    % No partial derivatives at all.
	
    % Preallocate point matrix for speed.
    xy=nan(2,nProj);

    for i=find(any(s.vis))
        % Get camera station.
        camStation=EO(:,i);
        center=camStation(1:3);
        ang=camStation(4:6);
        
        % Get inner orientation.
        camNo=cams(i);
        [pp,f,K,P,b,sz]=unpackio(IO(:,camNo),s.nK,s.nP);
        
        % Trim K and P
        kn0=find(K,1,'last');
        kn0=max([kn0;0]);
        if kn0<length(K)
            K=K(1:kn0);
        end
        pn0=find(P,1,'last');
        pn0=max([pn0;0]);
        if pn0==1
            pn0=2;
        end
        if pn0<length(P)
            P=P(1:pn0);
        end
	
        % Get object points visible in this image
        v=s.vis(:,i);
        obj=OP(:,v);

        % Get corresponding mark pts.
        cp=s.colPos(v,i);
        imPts=s.markPts(:,cp);
        
        % Compute image residuals
        sz=1./sz;
        camRes=fun(obj,center,ang,f,imPts,sz(1),pp,K,P,b);

        % Find out where to store residuals.
        xy(:,cp)=camRes;
    end
    f=xy(:);
else

end
