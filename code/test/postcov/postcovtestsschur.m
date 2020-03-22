dataDir=fullfile(fileparts(dbatroot),'data','test');

files={'camcaldemo.mat',...
      };

fix=1;

f=files{fix};
Z=load(fullfile(dataDir,f));
s=Z.s;
e=Z.E;

JTJ=e.final.weighted.J'*e.final.weighted.J;
N=JTJ;

nImages=size(s.EO.val,2);
nOP=size(s.OP.val,2);
nIP=size(s.IP.val,2);
raysPerOP=nIP/nOP;
sizeN=size(N,2);
sparsity=nnz(N)/numel(N);
selfCal=nnz(s.bundle.est.IO)>0;
if selfCal
    str='Yes';
else
    str='No';
end

% IO blocks.
[i,j]=ind2sub(size(s.bundle.est.IO),s.bundle.serial.IO.src);
bixIO=full(sparse(i,j,s.bundle.serial.IO.dest,size(s.bundle.est.IO,1),size(s.bundle.est.IO,2)));
% EO blocks.
[i,j]=ind2sub(size(s.bundle.est.EO),s.bundle.serial.EO.src);
bixEO=full(sparse(i,j,s.bundle.serial.EO.dest,size(s.bundle.est.EO,1),size(s.bundle.est.EO,2)));
% OP blocks.
[i,j]=ind2sub(size(s.bundle.est.OP),s.bundle.serial.OP.src);
bixOP=full(sparse(i,j,s.bundle.serial.OP.dest,size(s.bundle.est.OP,1),size(s.bundle.est.OP,2)));

%p=blkcolperm(JTJ,bixIO,bixEO,bixOP);
p=[bixEO(:);bixOP(:)];
p=p(p~=0);

xLim=[0,nnz(bixOP)+nnz(bixEO)];

nOP=nnz(bixOP);
nEO=nnz(bixEO);
p=p(1:nOP+nEO);
Nperm=JTJ(p,p);

Ninv=inv(Nperm);
NSigmaC=Ninv(1:nEO,1:nEO);
NSigmaP=Ninv(nEO+1:end,nEO+1:end);

C=Nperm(1:nEO,1:nEO);
B=Nperm(1:nEO,nEO+1:end);
P=Nperm(nEO+1:end,nEO+1:end);

SC=C-B*(P\B');
SP=P-B'*(C\B);

SigmaC=inv(SC);
SigmaP=inv(SP);

% Woodbury

invP=inv(P);

SigmaP2=invP+invP*B'*SigmaC*B*invP;

abserr=@(A,B)norm(A-B,'fro');
relerr=@(A,B)abserr(A,B)/norm(A,'fro');

invT=chol(SigmaC,'lower');

V=invT'*B*invP;

SigmaP3=invP+V'*V;

I=speye(size(P));

% Diagonals for 3-by-3 block diagonal Sigma_P
n=size(P,1);
d0=zeros(n,1);
dp1=zeros(n,1);
dp2=zeros(n,1);
dm1=zeros(n,1);
dm2=zeros(n,1);

SigmaPblocks=spalloc(n,n,3*n);

for i=1:3:n
    E=I(:,i+(0:2));
    invPE=P\E;
    EinvPE=E'*invPE;
    
    BinvPE=B*invPE;
    VE=invT'*BinvPE;
    EVTVE=VE'*VE;
    
    SigmaPblocks(i+(0:2),i+(0:2))=EinvPE+EVTVE;
end
