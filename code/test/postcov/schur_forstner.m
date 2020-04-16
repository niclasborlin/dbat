dataDir=fullfile(fileparts(dbatroot),'data','test');

files={'camcaldemo.mat',...
       'romabundledemo.mat',...
       'romabundledemo-selfcal.mat',...
       'romabundledemo-imagevariant.mat',...
       'sxb.mat',...
       'stpierre.mat',...
       'normier.mat',...
       'sewu-filt25.mat',...
       'sewu-filt35.mat'
      };

fix=2;

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

% Macro to extract blocks of matrix. Blocksizes are M followed by N.
Blk11=@(A,m,n)A(1:m,1:m);
Blk12=@(A,m,n)A(1:m,m+1:m+n);
Blk22=@(A,m,n)A(m+1:m+n,m+1:m+n);

% Partition to swap blocks of size M followed by N.
Swap=@(A,m,n)A([m+1:m+n,1:m],[m+1:m+n,1:m]);

sparsity=@(A)nnz(A)/numel(A);

Nf=Swap(Nperm,nEO,nOP);

% Nkk correspond to points, Npp to cameras
Nkk=Blk11(Nf,nOP,nEO);
Npp=Blk22(Nf,nOP,nEO);
Nkp=Blk12(Nf,nOP,nEO);
Npk=Nkp';

% Spp correspond to cameras, Skk to points.
Spp=Npp-Npk*(Nkk\Nkp);
Skk=Nkk-Nkp*(Npp\Npk);

Nfi=inv(Nf);
Cpp=Blk22(Nfi,nOP,nEO);
Ckk=Blk11(Nfi,nOP,nEO);

abserr=@(A,B)norm(A-B,'fro');
relerr=@(A,B)abserr(A,B)/norm(A,'fro');

% Try with restricted symamd reordering.
p=symamd(Nf);
p2=[p(p<=nOP),p(p>nOP)];
Nf2=Nf(p2,p2);
Nkk2=Blk11(Nf2,nOP,nEO);
Npp2=Blk22(Nf2,nOP,nEO);
Nkp2=Blk12(Nf2,nOP,nEO);
Npk2=Nkp2';

% Spp correspond to cameras, Skk to points.
Spp2=Npp2-Npk2*(Nkk2\Nkp2);
Skk2=Nkk2-Nkp2*(Npp2\Npk2);


%abserr(inv(Spp),Cpp)
%abserr(inv(Skk),Ckk)

% Wolf, Dewitt based on EO, OP ordering

A=Blk11(Nperm,nEO,nOP);
B=Blk12(Nperm,nEO,nOP);
D=Blk22(Nperm,nEO,nOP);

DL=chol(D,'lower');
DLi=inv(DL);

Di=DLi'*DLi;

C1=inv(A-B*Di*B');

BDi=B*Di;

C2=Di+BDi'*C1*BDi;

C22=inv(D-B'*inv(A)*B);
