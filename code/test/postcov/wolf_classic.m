dataDir=fullfile(fileparts(dbatroot),'data','test');

files={'camcaldemo.mat',...
       'romabundledemo.mat',...
       'romabundledemo-selfcal.mat',...
       'romabundledemo-imagevariant.mat',...
       'stpierre.mat',...
       'normier.mat',...
       'sewu-filt25.mat',...
       'sewu-filt35.mat'
      };

%fix=1;

f=files{fix}
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

startClock=now;

A=Blk11(Nperm,nEO,nOP);
D=Blk22(Nperm,nEO,nOP);
B=Blk12(Nperm,nEO,nOP);

if sparseB
    B=sparse(B);
else
    B=full(B);
end

%sparsity(B)

prepClock=now;
prepTime=(prepClock-startClock)*86400;

DL=chol(D,'lower');
DLi=inv(DL);

Di=DLi'*DLi;

diClock=now;
diTime=(diClock-prepClock)*86400;

C1=inv(full(A-B*Di*B'));

c1Clock=now;
c1Time=(c1Clock-diClock)*86400;

BDi=B*Di;

C2f=diagblkouter(C1,BDi,3);
C2=Di+C2f;

c2Clock=now;
c2Time=(c2Clock-c1Clock)*86400;

totalTime=(c2Clock-startClock)*86400;

allTimes=[diTime,c1Time,c2Time,totalTime]

abserr=@(A,B)norm(A-B,'fro');
relerr=@(A,B)abserr(A,B)/norm(A,'fro');

if fix==1
    C2check=extractdiagblocks(Blk22(inv(Nperm),nEO,nOP),3);

    abserr(C2check,C2)
    relerr(C2check,C2)
end
