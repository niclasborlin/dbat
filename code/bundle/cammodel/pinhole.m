function [q,dP,df,dpp,dfa,dfs,dP2,df2,dpp2,dfa2,dfs2]=pinhole(P,f,pp,fa,fs,cP,cf,cpp,cfa,cfs)
%PINHOLE Projection of points into a pinhole camera.
%
%   Q=PINHOLE(P,F,PP) projects the 3D points in the 3-by-N array P
%   into a pinhole camera with the scalar focal length F and 2-by-1
%   principal point PP. The 3D coordinates in P are in the 3D camera
%   coordinate system. The 2D coordinates in the 2-by-N array Q will
%   be in the same units as F and PP. The camera is assumed to have
%   unit aspect ratio and no skew.
%
%   Q=PINHOLE(P,F,PP,FA) computes the projected coordinates with an
%   aspect ratio of (F+FA)/F.
%
%   Q=PINHOLE(P,F,PP,FA,FS) computes the projected coordinates with a
%   skew parameter FS.
%
%   [Q,DP,DF,DPP,DFA,DFS]=PINHOLE(...) also computed the Jacobians of
%   Q with respect to the input parameters. If a subset of the
%   Jacobians are wanted, use PINHOLE(P,F,PP,FA,FS,cP,cF,cPP,cFA,cFS),
%   where the scalar logical parameters cP, cF, cPP, cFA, cFS control
%   whether the corresponding Jacobian is computed.

%   Undocumented 1: PINHOLE('SELFTEST') will run a self-test of the
%   analytical Jacobians.
%
%   Undocumented 2: [Q,DP,DF,DPP,DFA,DFS,DP2,DF2,DPP2,DFA2,DFS2]=...
%   will also compute numerical approximations of the Jacobians DP2,
%   DF2, DPP2, DFA2, and DFS2.

if nargin==1 && ischar(P), selftest, return; end

if nargin<4, fa=0; end
if nargin<5, fs=0; end
if nargin<6, cP=(nargout>1); end
if nargin<7, cf=(nargout>2); end
if nargin<8, cpp=(nargout>3); end
if nargin<9, cfa=(nargout>4); end
if nargin<10, cfs=(nargout>5); end

% Number of points.
n=size(P,2);

% Dummy Jacobians.
dP=sparse(2*n,3*n);
df=sparse(2*n,1);
dpp=sparse(2*n,1);
dfa=sparse(2*n,1);
dfs=sparse(2*n,1);

% Numerical Jacobians.
if nargout>6
    vec=@(x)x(:);
    fun=@(P)vec(pinhole(reshape(P,3,[]),f,pp,fa,fs));
    dP2=jacapprox(fun,P(:));
end

if nargout>7
    vec=@(x)x(:);
    fun=@(f)vec(pinhole(P,f,pp,fa,fs));
    df2=jacapprox(fun,f);
end

if nargout>8
    vec=@(x)x(:);
    fun=@(pp)vec(pinhole(P,f,pp,fa,fs));
    dpp2=jacapprox(fun,pp);
end

if nargout>9
    vec=@(x)x(:);
    fun=@(fa)vec(pinhole(P,f,pp,fa,fs));
    dfa2=jacapprox(fun,fa);
end

if nargout>10
    vec=@(x)x(:);
    fun=@(fs)vec(pinhole(P,f,pp,fa,fs));
    dfs2=jacapprox(fun,fs);
end

Px=P(1,:);
Py=P(2,:);
Pz=P(3,:);

% Pre-compute Px/Pz, Py/Pz
Pxz=Px./Pz;
Pyz=Py./Pz;

fx=f+fa;

% Pre-compute vectors that will be needed in the Jacobians.
fxPxz=fx*Pxz;
fPyz=f*Pyz;
% Only compute this one if necessary.
fsPyz=sparse(1,n);

q=[fxPxz+pp(1);
    fPyz+pp(2)];
if fs~=0
    fsPyz=fs*Pyz;
    q(1,:)=q(1,:)+fsPyz;
end

if cP
    if fs==0
        % Each block is 1/qz*[fx, 0, -fx*qx/qz]
        %                    [ 0, f, -f* qy/qz]
        v1=fx./Pz;
        v2=f./Pz;
        v3=-fxPxz./Pz;
        v4=-fPyz./Pz;

        % Put values in in-memory order.
        vv=reshape([v1;v2;v3;v4],2,[]);
        % Row indices for each block.  [1,2,1,2]      
        i0=reshape(repmat(0:2:2*n-1,2,1),1,[]);
        ii=[i0+1;i0+2];
        % Column indices for each block. [1,2,3,3]
        jj=repmat(0:3:3*n-1,4,1)+repmat([1;2;3;3],1,n);
        % Build Jacobian.
        dP=sparse(ii,jj,vv,2*n,3*n);
    else
        % Each block is 1/qz*[fx, fs, -(fx*qx+fs*qy)/qz]
        %                    [ 0,  f,     -f* qy/qz]
        v1=fx./Pz;
        v2=fs./Pz;
        v3=f./Pz;
        v4=-(fxPxz+fsPyz)./Pz;
        v5=-fPyz./Pz;

        % Put values in in-memory order.
        vv=reshape([v1;v2;v3;v4;v5],2,[]);
        % Row indices for each block. [1,1,2,1,2]
        ii=repmat(0:2:2*n-1,5,1)+repmat([1;1;2;1;2],1,n);
        % Column indices for each block. [1,2,2,3,3]
        jj=repmat(0:3:3*n-1,5,1)+repmat([1;2;2;3;3],1,n);
        % Build Jacobian.
        dP=sparse(ii,jj,vv,2*n,3*n);
    end
end

if cf
    df=reshape([Pxz;Pyz],2*n,1);
end

if cpp
    dpp=repmat(eye(2),n,1);
end

if cfa
    dfa=zeros(2*n,1);
    dfa(1:2:end)=Pxz';
end

if cfs
    dfs=zeros(2*n,1);
    dfs(1:2:end)=Pyz';
end


function selftest

% Compare the analytical and numerical Jacobians and report the
% maximum deviation. Should be below 1e-6.

n=12;
P=rand(3,n)+1;
f=rand+1;
pp=rand(2,1);
fa=rand;
fs=rand;

mx=-inf;
thres=1e-9;
% Self-test with minimum.
[~,dP,df,dpp,dfa,dfs,dP2,df2,dpp2,dfa2,dfs2]=pinhole(P,f,pp);
dif=full(max(max(abs(dP-dP2)))); mx=max([mx,dif]);
if dif>thres, warning('%s selftest: %s diff = %g.\n',mfilename,'dP',dif); end
dif=full(max(max(abs(df-df2)))); mx=max([mx,dif]);
if dif>thres, warning('%s selftest: %s diff = %g.\n',mfilename,'df',dif); end
dif=full(max(max(abs(dpp-dpp2)))); mx=max([mx,dif]);
if dif>thres, warning('%s selftest: %s diff = %g.\n',mfilename,'dpp',dif); end
dif=full(max(max(abs(dfa-dfa2)))); mx=max([mx,dif]);
if dif>thres, warning('%s selftest: %s diff = %g.\n',mfilename,'dfa',dif); end
dif=full(max(max(abs(dfs-dfs2)))); mx=max([mx,dif]);
if dif>thres, warning('%s selftest: %s diff = %g.\n',mfilename,'dfs',dif); end

% Self-test with fa.
[~,dP,df,dpp,dfa,dfs,dP2,df2,dpp2,dfa2,dfs2]=pinhole(P,f,pp,fa);
dif=full(max(max(abs(dP-dP2)))); mx=max([mx,dif]);
if dif>thres, warning('%s selftest: %s diff = %g.\n',mfilename,'dP',dif); end
dif=full(max(max(abs(df-df2)))); mx=max([mx,dif]);
if dif>thres, warning('%s selftest: %s diff = %g.\n',mfilename,'df',dif); end
dif=full(max(max(abs(dpp-dpp2)))); mx=max([mx,dif]);
if dif>thres, warning('%s selftest: %s diff = %g.\n',mfilename,'dpp',dif); end
dif=full(max(max(abs(dfa-dfa2)))); mx=max([mx,dif]);
if dif>thres, warning('%s selftest: %s diff = %g.\n',mfilename,'dfa',dif); end
dif=full(max(max(abs(dfs-dfs2)))); mx=max([mx,dif]);
if dif>thres, warning('%s selftest: %s diff = %g.\n',mfilename,'dfs',dif); end

% Self-test with all parameters.
[~,dP,df,dpp,dfa,dfs,dP2,df2,dpp2,dfa2,dfs2]=pinhole(P,f,pp,fa,fs);
dif=full(max(max(abs(dP-dP2)))); mx=max([mx,dif]);
if dif>thres, warning('%s selftest: %s diff = %g.\n',mfilename,'dP',dif); end
dif=full(max(max(abs(df-df2)))); mx=max([mx,dif]);
if dif>thres, warning('%s selftest: %s diff = %g.\n',mfilename,'df',dif); end
dif=full(max(max(abs(dpp-dpp2)))); mx=max([mx,dif]);
if dif>thres, warning('%s selftest: %s diff = %g.\n',mfilename,'dpp',dif); end
dif=full(max(max(abs(dfa-dfa2)))); mx=max([mx,dif]);
if dif>thres, warning('%s selftest: %s diff = %g.\n',mfilename,'dfa',dif); end
dif=full(max(max(abs(dfs-dfs2)))); mx=max([mx,dif]);
if dif>thres, warning('%s selftest: %s diff = %g.\n',mfilename,'dfs',dif); end

