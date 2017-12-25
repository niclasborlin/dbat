function [Q,db,dc,dP,db2,dc2,dP2]=affine(b,c,P,cb,cc,cP)
%AFFINE Apply affine transformation to points.
%
%   Q=AFFINE(B,C,P) applies the affine transformation defined by the
%   2-vector B to the points in the 2-by-N array P. The transformation
%   is defined with respect to the principal point in the 2-vector C
%   as
%
%     Q = [1+B(1), B(2)  
%               0,    1] * (P-C) + C.
%  
%   [Q,DB,DC,DP]=... also returns the Jacobians with respect to B, C,
%   and P, respectively. Use AFFINE(B,C,P,cB,cC,cP), where cB, cC, and
%   cP are logical, to individually control which Jacobians are
%   computated.
%
%   If cP==-1, DP is returned in compact form.

%   Undocumented 1: AFFINE('SELFTEST') will run a self-test of the
%   analytical Jacobians.
%
%   Undocumented 2: [Q,DB,DC,DP,DB2,DC2,DP2]=... will also compute
%   numerical approximations of the Jacobians DB2, DC2, DP2.

if nargin==1 && ischar(b), selftest, return; end

if nargin<4, cb=nargout>1; end
if nargin<5, cc=nargout>2; end
if nargin<6, cP=nargout>3; end

% Number of points.
n=size(P,2);
C=repmat(c,1,n);
pmc=P-C;
B=[1;0]*b';
IB=eye(2)+B;
Q=IB*pmc+C;

% Dummy Jacobians.
db=sparse(2*n,2);
dc=sparse(2*n,2);
dP=sparse(2*n,2*n);

% Numerical Jacobians
if nargout>4
    vec=@(x)x(:);
    fun=@(b)vec(affine(b,c,P));
    db2=jacapprox(fun,b);
end

if nargout>5
    vec=@(x)x(:);
    fun=@(c)vec(affine(b,c,P));
    dc2=jacapprox(fun,c);
end    

if nargout>6
    vec=@(x)x(:);
    fun=@(P)vec(affine(b,c,reshape(P,2,[])));
    dP2=jacapprox(fun,P(:));
end

if cb
    db=zeros(2*n,2);
    db(1:2:end,:)=pmc';
end

if cc
    dc=repmat(-B,n,1);
end

if cP
    if cP>0
        dP=kron(speye(n),sparse(IB));
    else
        dP=IB;
    end
end

    
function selftest

% Compare the analytical and numerical Jacobians and report the
% maximum deviation. Should be below 1e-9.

n=3;
P=rand(2,n);
b=rand(2,1);
c=rand(2,1);

mx=-inf;
thres=1e-9;
[~,db,dc,dP,db2,dc2,dP2]=affine(b,c,P);

dif=full(max(max(abs(db-db2)))); mx=max([mx,dif]);
if dif>thres, warning('%s selftest: %s diff = %g.\n',mfilename,'db',dif); end

dif=full(max(max(abs(dc-dc2)))); mx=max([mx,dif]);
if dif>thres, warning('%s selftest: %s diff = %g.\n',mfilename,'dc',dif); end

dif=full(max(max(abs(dP-dP2)))); mx=max([mx,dif]);
if dif>thres, warning('%s selftest: %s diff = %g.\n',mfilename,'dP',dif); end
