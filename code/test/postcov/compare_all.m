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
p=[bixOP(:);bixEO(:);bixIO(:)];
p=p(p~=0);

xLim=[0,nnz(bixOP)+nnz(bixEO)];

% Macro to extract blocks of matrix. Blocksizes are M followed by N.
Part11=@(A,m,n)A(1:m,1:m);
Part12=@(A,m,n)A(1:m,m+1:m+n);
Part22=@(A,m,n)A(m+1:m+n,m+1:m+n);

% Partition to swap blocks of size M followed by N.
Swap=@(A,m,n)A([m+1:m+n,1:m],[m+1:m+n,1:m]);

nOP=nnz(bixOP);
nEO=nnz(bixEO);
p=p(1:nOP+nEO);
Nperm=JTJ(p,p);

N0=Swap(Nperm,nOP,nEO);

N11=Part11(N0,nEO,nOP);
N12=full(Part12(N0,nEO,nOP));
N22=Part22(N0,nEO,nOP);

A=N11;
B=N12;
D=N22;

Di=inv(D);

N0inv=inv(N0);

tmp=A-B*Di*B';
tmpF=full(tmp);
C1f=inv(tmpF);

C1L=chol(C1f,'lower');

E=chol(D);
Ei=inv(E);

% Naive
BDi=B*Di;
C2naive=Di+BDi'*C1*BDi;

L=chol(Nperm,'lower');

