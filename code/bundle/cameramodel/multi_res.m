function [r,J]=multi_res(s,fun)
%MULTI_RES Compute residuals for multiple cameras.
%
%   MULTI_RES(S,FUN) computes residuals for all image observations.
%   The structure S contains the project information. FUN is a
%   function handle that will compute the image residuals for one
%   camera.
%
%   [R,J]=... also returns the requested Jacobian J.
%
%   See also: RES_EULER_BROWN_0, RES_EULER_BROWN_1, RES_EULER_BROWN_2,
%   RES_EULER_BROWN_3.

nPhotos=size(s.EO.val,2);
nObjs=size(s.OP.val,2);

% Total number of projected points.
nProj=nnz(s.IP.vis);

if nargout<2
    % No partial derivatives at all.
	
    % Preallocate point matrix for speed.
    xy=nan(2,nProj);

    for i=find(any(s.IP.vis))
        % Get camera station.
        camStation=s.EO.val(:,i);
        center=camStation(1:3);
        ang=camStation(4:6);
        
        % Get inner orientation.
        camNo=i;
        [pp,f,K,P,b]=unpackio(s.IO.val(:,camNo),s.IO.model.nK,s.IO.model.nP);
        sz=s.IO.sensor.pxSize(:,camNo);
        
        % Trim K and P
        K=trimkp(K,false);
        P=trimkp(P,true);
	
        % Get object points visible in this image
        v=s.IP.vis(:,i);
        obj=s.OP.val(:,v);

        % Get corresponding mark pts.
        cp=s.IP.ix(v,i);
        imPts=s.IP.val(:,cp);
        
        % Compute image residuals
        camRes=fun(obj,center,ang,f,imPts,sz(1),pp,K,P,b);

        % Find out where to store residuals.
        xy(:,cp)=camRes;
    end
    r=xy(:);
else
    % Create indices into the vector of unknowns = columns of J
    destIOcols=zeros(size(s.bundle.est.IO));
    destIOcols(s.bundle.deserial.IO.dest)=s.bundle.deserial.IO.src;
    destEOcols=zeros(size(s.bundle.est.EO));
    destEOcols(s.bundle.deserial.EO.dest)=s.bundle.deserial.EO.src;
    destOPcols=zeros(size(s.bundle.est.OP));
    destOPcols(s.bundle.deserial.OP.dest)=s.bundle.deserial.OP.src;
    % Preallocate residual vector.
    xy=nan(2*nProj,1);

    % Vectors for "sparse" IO: rows, cols, and values.
    IOrows=nan(numel(xy)*max(sum(s.bundle.est.IO,1)),1);
    IOcols=IOrows;
    IOvals=IOrows;
    
    % Vectors for "sparse" EO: rows, cols, and values.
    EOrows=nan(numel(xy)*max(sum(s.bundle.est.EO,1)),1);
    EOcols=EOrows;
    EOvals=EOrows;
    
    % Vectors for "sparse" OP: rows, cols, and values.
    OProws=nan(numel(xy)*max(sum(s.bundle.est.OP,1)),1);
    OPcols=OProws;
    OPvals=OProws;

    % Last used rows in IOrows/cols/vals, EOrows/cols/vals and
    % OProws/cols/vals, respectively.
    IOlast=0;
    EOlast=0;
    OPlast=0;
    
    % Last processed row in residual/Jacobian.
    jacLast=0;

    % For each camera with observations
    for i=find(any(s.IP.vis))
        % Get inner orientation.
        camNo=i;
        % Values
        [pp,f,K,P,b]=unpackio(s.IO.val(:,camNo),s.IO.model.nK,s.IO.model.nP);
        sz=s.IO.sensor.pxSize(:,camNo);
        % Do we need the partials?
        cIO=s.bundle.est.IO(:,camNo);
        [cpp,cf,cK,cP,cb]=unpackio(cIO,s.IO.model.nK,s.IO.model.nP);
        % Where should we store the partials?
        [ppIx,fIx,Kix,Pix,bIx]=unpackio(destIOcols(:,camNo), ...
                                        s.IO.model.nK,s.IO.model.nP);
        
        % Trim K and/or P unless we need the Jacobians.
        if ~any(cK)
            K=trimkp(K,false);
        end
        if ~any(cP)
            P=trimkp(P,true);
        end

        % Get external orientation.
        camStation=s.EO.val(:,i);
        center=camStation(1:3);
        ang=camStation(4:6);
        
        % Do we need the partials?
        cEO=s.bundle.est.EO(1:6,i);
        cC=cEO(1:3);
        cA=cEO(4:6);
        % Where should we store the partials?
        EOix=removezeros(destEOcols(:,i));
	
        % Get object points visible in this image
        v=s.IP.vis(:,i);
        obj=s.OP.val(:,v);
        % Do we need the partials?
        cOP=s.bundle.est.OP(:,v);
        % Where should we store the partials?
        OPix=removezeros(destOPcols(:,v));

        % Get corresponding mark pts.
        cp=full(s.IP.ix(v,i));
        imPts=s.IP.val(:,cp);
        
        % Compute image residuals
        [camRes,camJac]=fun(obj,center,ang,f,imPts,sz(1),pp,K,P,b, ...
                            any(cOP(:)),any(cC),any(cA),cf,false,false, ...
                            any(cpp),any(cK),any(cP),any(cb));

        % Rows within this block row in residual and Jacobian.
        blockRowIx=jacLast+(1:2*numel(cp));
        xy(blockRowIx)=camRes(:);

        if any(cIO)
            % Pack IO partials
            if any(cpp)
                if all(cpp)
                    [ii,jj,vv]=find(camJac.dU0);
                else
                    [ii,jj,vv]=find(camJac.dU0(:,cpp));
                    jj=jj-1+find(cpp);
                end

                % Where to store rows, cols, values in the respective IO vectors.
                IOblockIx=IOlast+(1:length(ii));
                % Store rows, cols, values.
                IOrows(IOblockIx)=jacLast+ii;
                IOcols(IOblockIx)=ppIx(jj);
                IOvals(IOblockIx)=vv;
                % Advance past used block.
                IOlast=IOlast+length(ii);
            end
            
            if cf
                [ii,~,vv]=find(camJac.dF);
                
                % Where to store rows, cols, values in the respective IO vectors.
                IOblockIx=IOlast+(1:length(ii));
                % Store rows, cols, values.
                IOrows(IOblockIx)=jacLast+ii;
                IOcols(IOblockIx)=fIx;
                IOvals(IOblockIx)=vv;
                % Advance past used block.
                IOlast=IOlast+length(ii);
            end
            
            if any(cK)
                % Ensure that cK has one 1-block followed by at
                % most one 0-block.
                cK1=find(cK);
                cK0=find(cK==0);
                if ~isempty(cK0) && min(cK0)<max(cK1)
                    error('Illegal cK vector');
                end
                
                % Will always work since we're always asking for
                % the leading columns.
                [ii,jj,vv]=find(camJac.dK(:,cK));

                % Where to store rows, cols, values in the respective IO vectors.
                IOblockIx=IOlast+(1:length(ii));
                % Store rows, cols, values.
                IOrows(IOblockIx)=jacLast+ii;
                IOcols(IOblockIx)=Kix(jj);
                IOvals(IOblockIx)=vv;
                % Advance past used block.
                IOlast=IOlast+length(ii);
            end
            
            if any(cP)
                % Ensure that cP has one 1-block followed by at most one 0-block.
                % First two elements count as one block.
                cP1=find(cP);
                cP0=find(cP==0);
                if any(cP0<=2) || (~isempty(cP0) && min(cP0)<max(cP1))
                    error('Illegal cP vector');
                end
                
                % Will always work since we're always asking for
                % the leading columns.
                [ii,jj,vv]=find(camJac.dP(:,cP));

                % Where to store rows, cols, values in the respective IO vectors.
                IOblockIx=IOlast+(1:length(ii));
                % Store rows, cols, values.
                IOrows(IOblockIx)=jacLast+ii;
                IOcols(IOblockIx)=Pix(jj);
                IOvals(IOblockIx)=vv;
                % Advance past used block.
                IOlast=IOlast+length(ii);
            end
            
            if any(cb)
                if all(cb)
                    [ii,jj,vv]=find(camJac.dB);
                else
                    [ii,jj,vv]=find(camJac.dB(:,cb));
                    jj=jj-1+find(cb);
                end

                % Where to store rows, cols, values in the respective IO vectors.
                IOblockIx=IOlast+(1:length(ii));
                % Store rows, cols, values.
                IOrows(IOblockIx)=jacLast+ii;
                IOcols(IOblockIx)=bIx(jj);
                IOvals(IOblockIx)=vv;
                % Advance past used block.
                IOlast=IOlast+length(ii);
            end
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

    IOrows=IOrows(1:IOlast);
    IOcols=IOcols(1:IOlast);
    IOvals=IOvals(1:IOlast);

    % Re-sort OP values after column.
    [OPcols,i]=sort(OPcols(1:OPlast));
    OProws=OProws(i);
    OPvals=OPvals(i);
    
    ii=[IOrows;EOrows(1:EOlast);OProws];
    jj=[IOcols;EOcols(1:EOlast);OPcols]; 
    vv=[IOvals;EOvals(1:EOlast);OPvals];
    
    J=sparse(ii,jj,vv,2*nProj,s.bundle.serial.n);
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

