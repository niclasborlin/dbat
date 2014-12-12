function [c,dIO,dEO,dOP]=pm_multirotconstr(IO,nK,nP,EO,cams,OP,vis,cIO,cEO,cOP)
%PM_MULTIROTCONSTR Rotation constraints for multiple camera stations.
%
%[c,dEO]=pm_multirotconstr(IO,nK,nP,EO,cams,OP,vis[,cIO,cEO,cOP])
%IO     - matrix with camera inner orientation as columns
%         pp - principal point [xp;yp].
%         f  - focal length.
%         K  - radial lens distortion coefficients.
%         P  - tangential lens distortion coefficients.
%         a  - with affine lens distortion coefficients.
%         u  - image units in x and y direction.
%nK     - number of radial koefficients.
%nP     - number of tangential koefficients.
%EO     - matrix with photo camera outer orientation as columns
%         C     - camera center [Xc,Yc,Zc].
%         rotP  - rotation parameters.
%         rotM  - rotation model.
%                 0 - Original Euler 'xyz'.
%                 1 - New Euler 'xyz'.
%                 2 - Euler 'zxz'.
%                 3 - Rodriguez abc.
%cams   - vector of camera numbers for each photo.
%OP     - 3 x nObj matrix with object coordinates
%vis    - sparse matrix indicating if obj point i is visible in photo j.
%cIO    - should we calculate partial derivatives w.r.t. internal
%         orientation? Scalar or matrix of same size as IO.
%cEO    - should we calculate partial derivatives w.r.t. external
%         orientation? Scalar or matrix of same size as EO.
%cOP    - should we calculate partial derivatives w.r.t. object points?
%xy     - projected points [x;y](:).
%dIO    - jacobian w.r.t. internal orientation.
%dEO    - jacobian w.r.t. external orientation.
%dOP    - jacobian w.r.t. object points.

% $Id$

if (nargin<8), cIO=(nargout>1); end
if (nargin<9), cEO=(nargout>2); end
if (nargin<10), cOP=(nargout>3); end

dIO=[];
dEO=[];
dOP=[];

nCam=size(IO,2);
nPhotos=size(EO,2);
nObjs=size(OP,2);

if (length(cams)==1)
    % Same camera for all photos.
    cams=repmat(cams,nPhotos,1);
end

% Total number of projected points.
nProj=nnz(vis);

if (all(~cIO(:)) & all(~cEO(:)) & all(~cOP(:)))
    % No partial derivatives at all.

    c=[];

    for i=find(any(cEO,1))
        % Get camera station.
        camStation=EO(:,i);
        center=camStation(1:3);
        % Rotation parameters.
        rotP=camStation(4:end-1);
        % Rotation model.
        rotM=camStation(end);
        
        % Project into image.
        cc=pm_rotconstr(rotP,rotM);
	
        c=[c;cc];
    end
else
    % Which IO partial derivatives are requested?
    if length(cIO)==1
        cIO=repmat(cIO,size(IO,1),size(IO,2));
    end
    
    % Preallocate IO jacobian.
    % Number of wanted parameters
    ioCols=nnz(cIO);
    
    % Which EO partial derivatives are requested?
    if length(cEO)==1
        cEO=repmat(cEO,size(EO));
        cEO(end,:)=false;
    end
    
    % Preallocate EO jacobian.
    % Number of wanted parameters
    eoCols=nnz(cEO);
    dEO=sparse(0,eoCols);
    
    % Create arrays of columns indices for EO derivatives.
    [~,ixRotP]=createeocolumnindices(cEO);
    
    % Preallocate OP jacobian.
    % Number of wanted parameters
    opCols=nnz(cOP);
	
    c=[];

    % Jacobian.
    dEOi=[];
    dEOj=[];
    dEOv=[];
    lc=length(c);

    for i=find(any(cEO,1))
        % Get camera station.
        camStation=EO(:,i);
        center=camStation(1:3);
        % Rotation parameters.
        rotP=camStation(4:end-1);
        % Rotation model.
        rotM=camStation(end);
        
        % Which outer orientation parameters are interesting?
        cRotP=cEO(4:end-1,i);
        
        % Compute constrint for this camera.
        [cc,dRotP]=pm_rotconstr(rotP,rotM);
	
        c=[c;cc];
        
        % EO jacobians.
        %dEO(ixPt,ixRotP(cRotP,i))=dRotP(:,cRotP);
        [ii,jj,vv]=find(dRotP(:,cRotP));
        ix2=ixRotP(cRotP,i);
        dEOi=[dEOi;lc+ii];
        dEOj=[dEOj;ix2(jj)];
        dEOv=[dEOv;vv];
        lc=length(c);
    end
    
    %dOP=dOPT';
    dIO=sparse(length(c),ioCols);
    dEO=sparse(dEOi,dEOj,dEOv,length(c),eoCols);
    dOP=sparse(length(c),opCols);
    % Verify
    %disp([nnz(dIO),ioMaxNnz]);
    %disp([nnz(dOP),opMaxNnz]);
end
