function [r,J]=multi_res(IO,EO,OP,s,fun)
%MULTI_RES Compute residuals for multiple cameras.
%
%   MULTI_RES(IO,EO,OP,S,FUN) computes residuals for all image
%   observations. The current IO, EO, OP values are stored in the
%   respective arrays. The structure S contains the project
%   information. FUN is a function handle that will compute the image
%   residuals for one camera.
%
%   [R,J]=... also returns the requested Jacobian J.
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
        K=trimkp(K,false);
        P=trimkp(P,true);
	
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
    r=xy(:);
else
    % Create indices into the vector of unknowns.
    [ixIO,ixEO,ixOP,n]=indvec([nnz(s.estIO),nnz(s.estEO),nnz(s.estOP)]);
    % Re-pack to original shapes.
    destIOcols=zeros(size(s.estIO));
    destIOcols(s.estIO)=ixIO;
    destEOcols=zeros(size(s.estEO));
    destEOcols(s.estEO)=ixEO;
    destOPcols=zeros(size(s.estOP));
    destOPcols(s.estOP)=ixOP;
    
    % Preallocate residual vector.
    xy=nan(2*nProj,1);

    % Pre-allocate JIO as full.
    JIO=nan(numel(xy),nnz(s.estIO));
    
    % Vectors for "sparse" EO: rows, cols, and values.
    EOrows=nan(numel(xy)*max(sum(s.estEO,1)),1);
    EOcols=EOrows;
    EOvals=EOrows;
    
    % Vectors for "sparse" OP: rows, cols, and values.
    OProws=nan(numel(xy)*max(sum(s.estOP,1)),1);
    OPcols=OProws;
    OPvals=OProws;

    % Last used rows in EOrows/cols/vals and OProws/cols/vals, respectively.
    EOlast=0;
    OPlast=0;
    
    % Last processed row in residual/Jacobian.
    jacLast=0;

    % For each camera with observations
    for i=find(any(s.vis))
        % Get inner orientation.
        camNo=cams(i);
        % Values
        [pp,f,K,P,b,sz]=unpackio(IO(:,camNo),s.nK,s.nP);
        % Do we need the partials?
        [cpp,cf,cK,cP,cb,csz]=unpackio(s.estIO(:,camNo),s.nK,s.nP);
        % Where should we store the partials?
        [ppIx,fIx,Kix,Pix,bIx,szIx]=unpackio(destIOcols,s.nK,s.nP);
        
        % Trim K and/or P unless we need the Jacobians.
        if ~any(cK)
            K=trimkp(K,false);
        end
        if ~any(cP)
            P=trimkp(P,true);
        end

        % Get external orientation.
        camStation=EO(:,i);
        center=camStation(1:3);
        ang=camStation(4:6);
        
        % Do we need the partials?
        cEO=s.estEO(1:6,i);
        cC=cEO(1:3);
        cA=cEO(4:6);
        % Where should we store the partials?
        EOix=removezeros(destEOcols(1:6,i));
	
        % Get object points visible in this image
        v=s.vis(:,i);
        obj=OP(:,v);
        % Do we need the partials?
        cOP=s.estOP(:,v);
        % Where should we store the partials?
        OPix=removezeros(destOPcols(:,v));

        % Get corresponding mark pts.
        cp=full(s.colPos(v,i));
        imPts=s.markPts(:,cp);
        
        % Compute image residuals
        sz=1./sz;
        
        if any(csz)
            warning('Internal error: should not estimate pixel size!');
        end
        
        [camRes,camJac]=fun(obj,center,ang,f,imPts,sz(1),pp,K,P,b, ...
                            any(cOP(:)),any(cC),any(cA),cf,false,false, ...
                            any(cpp),any(cK),any(cP),any(cb));

        % Rows within this block row in residual and Jacobian.
        blockRowIx=jacLast+(1:2*numel(cp));
        xy(blockRowIx)=camRes(:);

        % Store IO partials.
        switch nnz(cpp)
          case 2
            JIO(blockRowIx,ppIx)=camJac.dU0;
          case 0
            % Do nothing.
          otherwise
            JIO(blockRowIx,ppIx(cpp))=camJac.dU0(:,cpp);
        end
        if cf
            JIO(blockRowIx,fIx)=camJac.dF;
        end
        if any(cK)
            if all(cK)
                JIO(blockRowIx,Kix)=camJac.dK;
            else
                JIO(blockRowIx,Kix(cK))=camJac.dK(:,cK);
            end
        end
        if any(cP)
            if all(cP)
                JIO(blockRowIx,Pix)=camJac.dP;
            else
                JIO(blockRowIx,Pix(cP))=camJac.dP(:,cP);
            end
        end
        switch nnz(cb)
          case 2
            JIO(blockRowIx,bIx)=camJac.dB;
          case 0
            % Do nothing.
          otherwise
            JIO(blockRowIx,bIx(cb))=camJac.dB(:,cb);
        end

        if any(cEO)
            % Pack EO partials.
            switch nnz(cC)
              case 3
                dC=camJac.dQ0;
              case 0
                dC=zeros(size(camJac.dQ0,1),1);
              otherwise
                dC=camJac.dQ0(:,cC);
            end
            switch nnz(cA)
              case 3
                dA=camJac.dA;
              case 0
                dA=zeros(size(camJac.dA,1),1);
              otherwise
                dA=camJac.dA(:,cA);
            end
            % Dissect into row, column, value.
            [ii,jj,vv]=find([dC,dA]);

            % Where to store rows, cols, values in the respective EO vectors.
            EOblockIx=EOlast+(1:length(ii));
            % Store rows, cols, values.
            EOrows(EOblockIx)=jacLast+ii;
            EOcols(EOblockIx)=EOix(jj);
            EOvals(EOblockIx)=vv;
            % Advance past used block.
            EOlast=EOlast+length(ii);
        end

        if any(cOP(:))
            % Remove any unwanted OP partials (fixed OP).
            if all(cOP(:))
                dOP=camJac.dQ;
            else
                dOP=camJac.dQ(:,cOP(:));
            end
            % Dissect into row, column, value.
            [ii,jj,vv]=find(dOP);

            % Where to store rows, cols, values in the respective OP vectors.
            OPblockIx=OPlast+(1:length(ii));
        
            OProws(OPblockIx)=jacLast+ii;
            OPcols(OPblockIx)=OPix(jj);
            OPvals(OPblockIx)=vv;
        
            % Advance past used block.
            OPlast=OPlast+length(ii);
        end
        
        % Update Jacobian row and EO row.
        jacLast=jacLast+numel(camRes);
    end

    [IOrows,IOcols,IOvals]=find(JIO);

    % Re-sort OP values after column.
    [OPcols,i]=sort(OPcols(1:OPlast));
    OProws=OProws(i);
    OPvals=OPvals(i);
    
    ii=[IOrows;EOrows(1:EOlast);OProws];
    jj=[IOcols;EOcols(1:EOlast);OPcols]; 
    vv=[IOvals;EOvals(1:EOlast);OPvals];
    
    J=sparse(ii,jj,vv,2*nProj,n);
    r=xy(:);
end


function K=trimkp(K,firstIsPair)
% Remove trailing zeros from K or P vector. If firstIsPair is true,
% always keep or remove first two elements together.

i=find(K,1,'last');
if isempty(i)
    % Only zeros.
    K=zeros(0,size(K,2));
    return;
end

if firstIsPair && i==1
    % Cannot trim P2 if P1 is non-zero.
    i=2;
end

if i==length(K)
    % Last element is non-zero => cannot trim.
    return;
end

% Trim all zeros after last non-zero.
K=K(1:i);


function v=removezeros(v)
% Return v as a column vector with all zeros removed.

v=reshape(v(v~=0),[],1);

