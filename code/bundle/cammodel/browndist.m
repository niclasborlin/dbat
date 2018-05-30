function [d,ds,dpp,dK,dP,ds2,dpp2,dK2,dP2]=browndist(s,pp,K,P,cs,cpp,cK,cP)
%BROWNDIST Compute the lens distortion according to the Brown'71 model.
%
%   D=BROWNDIST(S,PP,K,P) returns a 2-by-N array D with the lens
%   distortion of the 2-by-N image point array S according to the
%   Brown (1971) lens distortion model. The 2-vector PP contain the
%   principal point. The vectors K and P contain the radial and
%   tangential coefficients, respectively. The K vector may be of any
%   length. The P vector may be of any length except one.
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
%   where
%
%       w = s - pp,
%
%   and r=norm(w). The tangential part is given by
%
%       Dt = [ P(1) ( r^2 + 2x^2 ) + 2P(2)xy;
%              2P(1)xy + P(2) ( r^2 + 2y^2 )  ] (1 + P(3)r^2 + ...).
%
%   [D,dS,dPP,dK,dP]=BROWNDIST(S,PP,K,P) also computes the Jacobians
%   of D with respect to the variables in S, PP, K, and P,
%   respectively. If a subset of the Jacobians are wanted, use
%   BROWNDIST(S,PP,K,P,cS,cPP,cK,cP), where the logical parameters cS,
%   cPP, cK, and cP control whether the corresponding Jacobian is
%   computed. cS must be scalar, whereas the others may be scalar or
%   of the same size as their respective parameter. As a special
%   case, if cS is negative, dS is returned as a 2-by-2-by-N dense
%   array with the 2-by-2 diagonal blocks of dS.
%
%   References: Brown (1971), "Close-range camera calibration".
%       Photogrammetric Engineering, 37(8): 855-866.

%   Undocumented 1: BROWNDIST('SELFTEST') will run a self-test of
%   the analytical Jacobians.
%
%   Undocumented 2: [D,dS,dPP,dK,dP,dS2,dSPP2,dK2,dP2]=... will also
%   compute numerical approximations of the Jacobians dS2, dK2, and
%   dP2.

if nargin==1 && ischar(s), selftest, return; end

if nargin<5, cs=(nargout>1); end
if nargin<6, cpp=(nargout>2); end
if nargin<7, cK=(nargout>3); end
if nargin<8, cP=(nargout>4); end

if ~isscalar(cs), error('%s: cs parameter must be scalar',mfilename); end

% We need Jacobian w.r.t. w in two cases.
cw=cs || any(cpp);

% Number of points.
n=size(s,2);

ds=sparse(2*n,2*n);
dpp=sparse(2*n,2);
dK=sparse(2*n,length(K));
dP=sparse(2*n,length(P));
ds2=sparse(2*n,2*n);
dpp2=sparse(2*n,2);
dK2=sparse(2*n,length(K));
dP2=sparse(2*n,length(P));
drdw=sparse(2,2*n);
dtdw=sparse(2,2*n);

if nargout>5
    % Numerical approximation of dw.
    vec=@(x)x(:);
    f=@(s)vec(browndist(reshape(s,2,[]),pp,K,P));
    ds2=jacapprox(f,s(:));
end

if nargout>6
    % Numerical approximation of dw.
    vec=@(x)x(:);
    f=@(pp)vec(browndist(s,pp,K,P));
    dpp2=jacapprox(f,pp,1e-8);
end

if nargout>7
    % Numerical approximation of dK.
    vec=@(x)x(:);
    f=@(K)vec(browndist(s,pp,K,P));
    dK2=jacapprox(f,K);
end

if nargout>8
    % Numerical approximation of dP.
    vec=@(x)x(:);
    f=@(P)vec(browndist(s,pp,K,P));
    dP2=jacapprox(f,P);
end

% Subtract principal point.
w=s-repmat(pp,1,n);
    
% Split w into components.
x=w(1,:);
x2=x.^2;
y=w(2,:);
y2=y.^2;
xy=x.*y;
% Radial distance squared.
r2=x2+y2;

if isempty(s)
    % No points.
    d=zeros(size(s));
    return;
end

if isempty(K)
    % No radial distortion.
    dr=sparse(2,n);
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

    if any(cK)
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
    % No tangential distortion.
    dt=sparse(2,n);
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

    if any(cP)
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

if cw
    % Analytical Jacobian w.r.t. w.
    
    % Add radial and tangential parts.
    dgdw=drdw+dtdw;
    
    if any(cpp)
        % Stack the 2-by-2 blocks on top of each other.
        dpp=-reshape([dgdw(:,1:2:end),dgdw(:,2:2:end)],[],2);
        if ~all(cpp)
            dpp=dpp(:,cpp);
        end
    end
    
    if cs
        if cs>0
            % Convert to 2-by-2 block diagonal.
            [i,j,v]=find(dgdw);
            ds=sparse(i+floor((j-1)/2)*2,j,v,2*n,2*n);
        else
            % Return in packed format.
            ds=reshape(full(dgdw),2,2,n);
        end
    end
end


function selftest
% Compare the analytical and numerical Jacobians and report the
% maximum deviation. Should be below 1e-6.

s=rand(2,8);
pp=rand(2,1);
K=rand(5,1);
P=rand(5,1);

mx=-inf;
thres=1e-6;
for kLen=0:length(K)
    for pLen=[0,2:length(P)]
        [d,ds,dpp,dK,dP,ds2,dpp2,dK2,dP2]=browndist(s,pp,K(1:kLen),P(1:pLen));
        dif=full(max(max(abs(ds-ds2)))); mx=max([mx,dif]);
        if dif>thres, warning('%s selftest: %s diff = %g.\n',mfilename,'ds',dif); end
        dif=full(max(max(abs(dpp-dpp2)))); mx=max([mx,dif]);
        if dif>thres, warning('%s selftest: %s diff = %g.\n',mfilename,'dpp',dif); end
        dif=full(max(max(abs(dK-dK2)))); mx=max([mx,dif]);
        if dif>thres, warning('%s selftest: %s diff = %g.\n',mfilename,'dK',dif); end
        dif=full(max(max(abs(dP-dP2)))); mx=max([mx,dif]);
        if dif>thres, warning('%s selftest: %s diff = %g.\n',mfilename,'dP',dif); end
    end
end
if mx<=thres
    fprintf('%s selftest: Maximum diff = %g, max expected=%g, OK.\n',mfilename,mx,...
            thres);
end
