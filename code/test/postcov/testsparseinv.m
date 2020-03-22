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

nOP=nnz(bixOP);
nEO=nnz(bixEO);
p=p(1:nOP+nEO);
Nperm=JTJ(p,p);

addpath(fullfile(dbatroot,'test','sparseinv'),'-end');

A=Nperm;

[Z, Zpattern] = sparseinv (A) ;
S = inv (A) ;
err = norm (Zpattern .* (Z - S), 1) / norm (S, 1)
