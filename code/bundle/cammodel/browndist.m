function [d,dw,dK,dP]=browndist(w,K,P,cw,cK,cP)
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

if nargin<4, cp=(nargout>1); end
if nargin<5, cK=(nargout>2); end
if nargin<6, cP=(nargout>3); end

% Quick return for 'no distortion'.
if isempty(K) && isempty(P)
    % Create sparse zero matrices of the correct sizes.
    d=sparse(size(w,1),size(w,2));
    dp=sparse(numel(p),numel(p));
    dK=sparse(numel(p),0);
    dP=sparse(numel(p),0);
    return;
end

dw=[];
dK=[];
dP=[];

% Number of points.
n=size(w,2);

% Radial distance squared.
r2=sum(w.^2,1);

% Create r2.^[1..nK], where nK is the number of K coefficients.

if isempty(K)
    dr=0;
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
end

if isempty(P)
    dt=0;
else
    % Split w into components.
    x=w(1,:);
    y=w(2,:);
    xy=x.*y;
    % Compute basic distortion.
    dtx=(P(1)*(r2+2*x.^2)+2*P(2)*xy);
    dty=(P(2)*(r2+2*y.^2)+2*P(1)*xy);
    dt=[dtx;dty];
    % Do we have a scaled tangential distortion?
    if length(P)>2
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
end

d=dr+dt;

% TODO Jacobians
