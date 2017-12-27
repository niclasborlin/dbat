function [ld,dIO,dp]=multilensdist(p,IO,nK,nP,cams,cIO,cp)
%MULTILENSDIST Compute lens distortion in multiple images.
%
%   MULTILENSDIST(P,IO,nK,nP) computes the lens distortion of the
%   points in the 2-by-N array P using the internal orientation
%   parameters in the 16-vector IO. The scalars nK and nP contain the
%   number of radial and tangential distortion parameters,
%   respectively, stored in IO.
%
%   MULTILENSDIST(P,IO,nK,nP,CAMS), where IO is 16-by-nCams, and
%   CAMS is as N-vector, computes the lens distortion for points in
%   multiple cameras. The vector CAMS indicate which IO column
%   holds the parameters for the camera the imaged each point.
%
%   [LD,dIO,dP]=MULTILENSDIST(...) also computes the Jacobians dIO and
%   dP with respect to the IO parameters, and the points P,
%   respectively. If a subset of the Jacobians is wanted, use
%   [LD,dIO,dP]=MULTILENSDIST(P,IO,nK,nP,CAMS,cIO,cP), where cIO
%   and cP indicate what Jacobians are requested. cIO can be scalar
%   or of the same size as IO. cP must be scalar.
%
%   Each column of the IO array holds the following parameters:
%       pp - principal point (2-vector),
%       f  - focal length (scalar),
%       K  - radial lens distortion coefficients (nK-vector),
%       P  - tangential lens distortion coefficients (nP-vector),
%       b  - aspect and skew parameters (2-vector),
%       s  - sensor width and height (2-vector),
%       i  - image width and height (2-vector),
%       r  - image resolution (2-vector)
%
%   See also: BROWNDIST.

if nargin==1 && ischar(p), selftest, return, end

if nargin<5, cams=1; end
if nargin<6, cIO=(nargout>1); end
if nargin<7, cp=(nargout>2); end

if ~(isscalar(cIO) || all(size(cIO)==size(IO)))
    error('%s: Bad size of cIO parameter',mfilename);
end

if ~isscalar(cp), error('%s: cp parameter must be scalar',mfilename); end
    
% Total number of points.
n=size(p,2);

if isscalar(cams)
    % Same camera for all points.
    cams=repmat(cams,n,1);
end

if ~(any(cIO(:)) || cp)
    % No partial derivatives.
    dIO=zeros(2*n,0);
    dp=sparse(2*n,2*n);
    
    % Preallocate correction matrix for speed.
    ld=nan(size(p));
    
    for i=1:max(cams)
        % Get inner orientation.
        [pp,~,K,P]=unpackio(IO(:,i),nK,nP);
	
        % Get points taken with this camera.
        ix=cams==i;
        q=p(:,ix);
        
        % Trim K and P
        kn0=find(K,1,'last');
        kn0=max([kn0;0]);
        if kn0<length(K)
            K=K(1:kn0);
        end
        pn0=find(P,1,'last');
        pn0=max([pn0;0]);
        if pn0<length(P) && pn0~=1
            P=P(1:pn0);
        end
            
        % Compute lens distortion for these points.
        lens=browndist(q,pp,K,P);
        
        % Store.
        ld(:,ix)=lens;
    end
else
    % Which IO partial derivatives are requested?
    if isscalar(cIO)
        cIO=repmat(cIO,size(IO));
    end
    
    % Pre-allocate full IO Jacobian since it will be dense and nnz(cIO)
    % will be small.
    dIO=zeros(2*n,nnz(cIO));
    
    % Create arrays of columns indices for IO derivatives.
    [ixpp,~,ixK,ixP]=createiocolumnindices(cIO,nK,nP);

    if cp
        % Jacobian w.r.t. points will be block diagonal with 2-by-2
        % blocks. Preallocate array to hold blocks.
        dpBlock=zeros(2,2,n);
        % Ask browndist for packed Jacobian.
        cp=-1;
    else
        dp=sparse(2*n,2*n);
    end
    
    % Preallocate correction matrix for speed.
    ld=nan(size(p));
    
    for i=1:max(cams)
        % Get inner orientation.
        [pp,~,K,P]=unpackio(IO(:,i),nK,nP);
	
        % Which inner orientation parameters are interesting?
        [cpp,~,cK,cP]=unpackio(cIO(:,i),nK,nP);
        
        % Get points taken with this camera.
        ix=find(cams==i);
        q=p(:,ix);
        
        % Trim K and P
        tixK=ixK;
        kn0=find(K,1,'last');
        ckn0=find(cK,1,'last');
        kn0=max([kn0;ckn0;0]);
        if kn0<length(K)
            K=K(1:kn0);
            cK=cK(1:kn0);
            tixK=tixK(1:kn0);
        end
        tixP=ixP;
        pn0=find(P,1,'last');
        cpn0=find(cP,1,'last');
        pn0=max([pn0;cpn0;0]);
        if pn0<length(P) && pn0~=1
            P=P(1:pn0);
            cP=cP(1:pn0);
            tixP=tixP(1:pn0);
        end
            
        % Lens distortion.
        [lens,dldq,dldpp,dldK,dldP]=browndist(q,pp,K,P,cp,cpp,cK,cP);
        
        % Correct for lens distortion.
        ld(:,ix)=lens;
        
        % Calculate row indices in Jacobians.
        ixRow=[(ix-1)*2+1;(ix-1)*2+2];
        
        % IO jacobians.
        if any(cpp)
            dIO(ixRow,ixpp(cpp,i))=dldpp(:,cpp);
        end
        if any(cK)
            dIO(ixRow,tixK(cK,i))=dldK(:,cK);
        end
        if any(cP)
            dIO(ixRow,tixP(cP,i))=dldP(:,cP);
        end
        if cp && nnz(ix)>0
            % Put diagonal blocks in stack.
            dpBlock(:,:,ix)=dldq;
        end
    end
end

if cp
    % Convert the diagonal blocks to block-diagonal matrix.
    [i,j,v]=find(reshape(dpBlock,2,2*n));
    dp=sparse(i+floor((j-1)/2)*2,j,v,2*n,2*n);
end

function selftest
% Compare the analytical and numerical Jacobians and report the
% maximum deviation. Should be below 1e-9.

vec=@(x)x(:);

rng('default');
nCams=35;
IO=rand(16,nCams);
nK=3;
nP=2;
n=15;
p=rand(2,n);
cams=floor(rand(1,n)*nCams)+1;
fIO=@(IO)vec(multilensdist(p,reshape(IO,16,[]),nK,nP,cams));
fp=@(p)vec(multilensdist(reshape(p,2,[]),IO,nK,nP,cams));
dIO2=jacapprox(fIO,IO(:));
dp2=jacapprox(fp,p(:));
[~,dIO,dp]=multilensdist(reshape(p,2,[]),IO,nK,nP,cams);
mx=-inf;
mx=max(mx,full(max(max(abs(dIO-dIO2)))));
mx=max(mx,full(max(max(abs(dp-dp2)))));
thres=1e-9;
if mx<thres
    fprintf('%s selftest: Maximum diff = %g, max expected=%g, OK.\n',mfilename,mx,...
            thres);
else
    warning('%s selftest: Maximum diff = %g, max expected=%g.\n',mfilename,mx,...
            thres);
end
