function [xy,dIO,dEO,dOP]=multieulerpinhole(IO,nK,nP,EO,cams,OP,vis,cIO,cEO,cOP)
%MULTIEULERPINHOLE Pinhole camera projection into multiple camera stations.
%
%[xy,dIO,dEO,dOP]=multieulerpinhole(IO,nK,nP,EO,cams,OP,vis[,cIO,cEO,cOP])
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
%         ang   - Euler angles.
%         ax    - Euler axis sequence
%                 0 - 'xyz'.
%                 1 - 'zxz'.
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


if (nargin<8), cIO=(nargout>1); end
if (nargin<9), cEO=(nargout>2); end
if (nargin<10), cOP=(nargout>3); end

dIO=[];
dEO=[];
dOP=[];

nPhotos=size(EO,2);
nObjs=size(OP,2);

if (length(cams)==1)
    % Same camera for all photos.
    cams=repmat(cams,nPhotos,1);
end

% Total number of projected points.
nProj=nnz(vis);

if (all(~cIO(:)) && all(~cEO(:)) && all(~cOP(:)))
    % No partial derivatives at all.
	
    % Preallocate point matrix for speed.
    xy=zeros(2*nProj,1);
    % Last used row index into xy.
    xyBase=0;

    for i=find(any(vis))
        % Get camera station.
        camStation=EO(:,i);
        center=camStation(1:3);
        ang=camStation(4:6);
        seq='xyz';
        
        % Get inner orientation.
        camNo=cams(i);
        [pp,f]=unpackio(IO(:,camNo),nK,nP);
	
        % Get object points visible in this image
        v=vis(:,i);
        obj=OP(:,v);
	
        % Project into image.
        imPt=eulerpinhole(pp,f,obj,center,ang,seq);
	
        % Find out where to store points.
        [ixPt,xyBase]=indvec(nnz(v)*2,xyBase);
        xy(ixPt)=imPt(:);
    end
    xy=reshape(xy,2,[]);
else
    % Which IO partial derivatives are requested?
    if isscalar(cIO)
        cIO=repmat(cIO,size(IO));
    end
    
    % Pre-allocate full IO Jacobian since it will be dense and nnz(cIO)
    % will be small.
    dIO=zeros(nProj*2,nnz(cIO));
    
    % Create arrays of columns indices for IO derivatives.
    [ixpp,ixf]=createiocolumnindices(cIO,nK,nP);
    
    
    % Which EO partial derivatives are requested?
    if isscalar(cEO)
        cEO=repmat(cEO,6,nPhotos);
    end
    
    % Max number of non-zero elements of EO Jacobian.
    eoMaxNnz=nProj*2*max([0,sum(cEO)]);
    % Collect (i,j) indices and values of dEO during loop.
    dEOi=zeros(eoMaxNnz,1);
    dEOj=zeros(eoMaxNnz,1);
    dEOv=zeros(eoMaxNnz,1);
    dEOix=0;
    
    % Create arrays of columns indices for EO derivatives.
    [ixC,ixAng]=createeocolumnindices(cEO);
    
    
    % Which OP partial derivatives are requested?
    if isscalar(cOP)
        cOP=repmat(cOP,3,nObjs);
    end
    
    % Max number of non-zero elements.
    opMaxNnz=nProj*2*max(sum(cOP));
    % Collect (i,j) indices and values of dOP during loop.
    dOPi=zeros(opMaxNnz,1);
    dOPj=zeros(opMaxNnz,1);
    dOPv=zeros(opMaxNnz,1);
    dOPix=0;

    % Create array of columns indices for OP derivatives.
    ixOP=createopcolumnindices(cOP);
	
    % Preallocate point matrix for speed.
    xy=zeros(2*nProj,1);
    % Last used row index into xy.
    xyBase=0;

    for i=find(any(vis))
        % Get camera station.
        camStation=EO(:,i);
        center=camStation(1:3);
        ang=camStation(4:6);
        seq='xyz';
        
        % Which outer orientation parameters are interesting?
        cC=cEO(1:3,i);
        cAng=cEO(4:6,i);
        
        % Get inner orientation.
        camNo=cams(i);
        [pp,f]=unpackio(IO(:,camNo),nK,nP);
	
        % Which inner orientation parameters are interesting?
        [cpp,cf]=unpackio(cIO(:,camNo),nK,nP);
        
        % Get object points visible in this image
        v=vis(:,i);
        obj=OP(:,v);
        % ...and which partial derivatives to calculate...
        cObj=cOP(:,v);
        % ... and in which columns to store them.
        ixObj=ixOP(:,v);
        
        % Project into image.
        [imPt,dpp,df,dO,dC,dAng]=eulerpinhole(pp,f,obj,center,ang, ...
                                              seq,any(cpp),cf, ...
                                              cObj,any(cC),any(cAng));
	
        % Find out where to store points.
        [ixPt,xyBase]=indvec(nnz(v)*2,xyBase);
        xy(ixPt)=imPt(:);
        
        % IO jacobians.
        if (any(cpp))
            dIO(ixPt,ixpp(cpp,camNo))=dpp(:,cpp);
        end
        if (cf)
            dIO(ixPt,ixf(camNo))=df;
        end
        
        % EO jacobians.
        if (any(cC))
            %dEO(ixPt,ixC(cC,i))=dC(:,cC);
            [ii,jj,vv]=find(dC(:,cC));
            ix2=ixC(cC,i);
            dEOi(dEOix+(1:length(ii)))=ixPt(ii);
            dEOj(dEOix+(1:length(ii)))=ix2(jj);
            dEOv(dEOix+(1:length(ii)))=vv;
            dEOix=dEOix+length(ii);
        end
        if (any(cAng))
            %dEO(ixPt,ixAng(cAng,i))=dAng(:,cAng);
            [ii,jj,vv]=find(dAng(:,cAng));
            ix2=ixAng(cAng,i);
            dEOi(dEOix+(1:length(ii)))=ixPt(ii);
            dEOj(dEOix+(1:length(ii)))=ix2(jj);
            dEOv(dEOix+(1:length(ii)))=vv;
            dEOix=dEOix+length(ii);
        end
        
        % OP jacobians.
        if (any(cObj(:)))
            %dOP(ixPt,ixObj(cObj))=dO;
            [ii,jj,vv]=find(dO);
            ix2=ixObj(cObj);
            dOPi(dOPix+(1:length(ii)))=ixPt(ii);
            dOPj(dOPix+(1:length(ii)))=ix2(jj);
            dOPv(dOPix+(1:length(ii)))=vv;
            dOPix=dOPix+length(ii);
        end
    end
    
    if dEOix<length(dEOi)
        % Trim vectors before calling sparse.
        dEOi(dEOix+1:end)=[];
        dEOj(dEOix+1:end)=[];
        dEOv(dEOix+1:end)=[];
    end
    dEO=sparse(dEOi,dEOj,dEOv,nProj*2,nnz(cEO));

    if dOPix<length(dOPi)
        % Trim vectors before calling sparse.
        dOPi(dOPix+1:end)=[]; 
        dOPj(dOPix+1:end)=[]; 
        dOPv(dOPix+1:end)=[];
    end
    dOP=sparse(dOPi,dOPj,dOPv,nProj*2,nnz(cOP));
    xy=reshape(xy,2,[]);
end
