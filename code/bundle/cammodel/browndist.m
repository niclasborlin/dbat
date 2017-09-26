function [d,dw,dK,dP,dw2,dK2,dP2]=browndist(w,K,P,cw,cK,cP)
%BROWNDIST Compute the lens distortion according to the Brown'71 model.
%
%   D=BROWNDIST(W,K,P) returns a 2-by-N array D with the lens
%   distortion of the 2-by-N image point array W according to the
%   Brown (1971) lens distortion model. The vectors K and P contain
%   the radial and tangential coefficients, respectively. The
%   coordinates in the point array W is relative to the principal
%   point. The K vector may be of any length. The P vector may be of
%   any length except one.
%
%   The Brown (1971) lens distortion model splits the distortion into
%   a radial and tangential part
%
%       D = Dr + Dt.
%        
%   For each column w=[x;y] of W, the radial part is given by
%
%       Dr = w * ( K(1)r^2 + K(2)r^4 + ... ),
%
%   where r=norm(w). The tangential part is given by
%
%       Dt = [ P(1) ( r^2 + 2x^2 ) + 2P(2)xy;
%              2P(1)xy + P(2) ( r^2 + 2y^2 )  ] (1 + P(3)r^2 + ...).
%
%   [D,dW,dK,dP]=BROWNDIST(W,K,P) also computes the Jacobians of D
%   with respect to the variables in W, K, and P, respectively. If a
%   subset of the Jacobians are wanted, use BROWNDIST(W,K,P,cW,cK,cP),
%   where the logical parameters cW, cK, and cP control whether the
%   corresponding Jacobian is computed.
%
%   References: Brown (1971), "Close-range camera calibration".
%       Photogrammetric Engineering, 37(8): 855-866.

%   Undocumented: [D,dW,dK,dP,dW2,dK2,dP2]=... will also approximate
%   the numerical Jacobians dW2, dK2, dP2.

if nargin<4, cw=(nargout>1); end
if nargin<5, cK=(nargout>2); end
if nargin<6, cP=(nargout>3); end

dw=[];
dK=[];
dP=[];
dw2=[];
dK2=[];
dP2=[];
drdw=0;
dtdw=0;

% Number of points.
n=size(w,2);

% Split w into components.
x=w(1,:);
x2=x.^2;
y=w(2,:);
y2=y.^2;
xy=x.*y;
% Radial distance squared.
r2=x2+y2;

if isempty(K)
    dr=0;
    dK=zeros(2*n,0);
else
    % Create r2 exponent matrix.
    nK=length(K);
    r2e=repmat(1:nK,n,1);

    % r2.^[1,2,3,...]
    r2k=repmat(r2',1,nK).^r2e;

    % Inner product K(1)*r^2+K(2)*r^4+K(3)*r^6+...
    Kr=(r2k*K)';

    % Compute radial distortion
    dr=w.*repmat(Kr,2,1);

    if cK
        % Analytical Jacobian w.r.t. K.
    
        % Each 2-by-nK block row is w(:,i)*r2k(i,:).
        r2kr=reshape(repmat(r2k(:)',2,1),size(r2k,1)*2,[]);
        wr=repmat(w(:),1,nK);
        dK=r2kr.*wr;
    end
    
    if cw
        % Analytical Jacobian of dr w.r.t. w.
        
        r2km1=repmat(r2',1,nK).^(r2e-1).*r2e;
        Kdr=(r2km1*K)';
        
        drdw=zeros(2,2*n);
        drdw(1,:)=reshape([Kr+2*Kdr.*x2;2*Kdr.*xy],[],1)';
        drdw(2,:)=reshape([2*Kdr.*xy;Kr+2*Kdr.*y2;],[],1)';
    end
end

nP=length(P);

if nP==0
    dt=0;
    dP=zeros(2*n,0);
else

    % Construct 2-by-2 block rows of w'*w*eye(2)+2*w*w'.
    Aw=zeros(2*n,2);
    Aw(:,1)=reshape([r2+2*x.^2;2*xy],[],1);
    Aw(:,2)=reshape([2*xy;r2+2*y.^2],[],1);
    AwQ=Aw*P(1:2);
    dt=reshape(AwQ,2,[]);

    % Do we have a scaled tangential distortion?
    if nP>2
        S=P(3:end);
        nS=length(S);

        r2e=repmat(1:nS,n,1);

        % r2.^[1,2,3,...]
        r2s=repmat(r2',1,nS).^r2e;

        % Inner product K(1)*r^2+K(2)*r^4+K(3)*r^6+...
        Sr=(r2s*S)';

        % Scale the tangential distortion.
        dt=dt.*repmat(1+Sr,2,1);
    end

    if cP
        % Analytical Jacobian w.r.t. P.

        switch nP
          case 0
            dP=zeros(2*n,0);
          case 2
            dP=Aw;
          otherwise
            dP=zeros(2*n,nP);
            % Compute w.r.t. P1, P2.
            dP(:,1:2)=Aw.*repmat(reshape(repmat(1+Sr,2,1),[],1),1,2);
            dP(:,3:end)=repmat(AwQ,1,nS).*kron(r2s,ones(2,1));
        end
    end
    
    if cw
        % Analytical Jacobian of dr w.r.t. w.

        if nP==0
            dtdw=0;
        else
            dtdw=reshape(P(1)*[6*x;2*y;2*y;2*x]+P(2)*[2*y;2*x;2*x;6*y],2,[]);
            if nP>2
                dtdw=dtdw.*reshape(repmat(1+Sr,4,1),2,[]);
                
                r2sm1=repmat(r2',1,nS).^(r2e-1).*r2e;
                Sdr=(r2sm1*S)';
                
                dtdw=dtdw+reshape(repmat(reshape(AwQ,2,[]),2,1),2,[]).* ...
                     repmat(reshape(w.*repmat(2*Sdr,2,1),[],1)',2,1);
            end
        end
    end
end

d=dr+dt;

if nargout>4
    % Numerical approximation of dw.
    vec=@(x)x(:);
    f=@(w)vec(browndist(reshape(w,2,[]),K,P));
    dw2=jacapprox(f,w(:));
end

if nargout>5
    % Numerical approximation of dK.
    vec=@(x)x(:);
    f=@(K)vec(browndist(w,K,P));
    dK2=jacapprox(f,K);
end

if nargout>6
    % Numerical approximation of dP.
    vec=@(x)x(:);
    f=@(P)vec(browndist(w,K,P));
    dP2=jacapprox(f,P);
end

if cw
    % Analytical Jacobian w.r.t. w.
    
    % Add radial and tangential parts.
    dgdw=drdw+dtdw;
    
    [i,j,v]=find(dgdw);
    
    % Convert to 2-by-2 block diagonal.
    dw=sparse(i+floor((j-1)/2)*2,j,v,2*n,2*n);
end
