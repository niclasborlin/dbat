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

if 0
for fix=1:length(files)
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
    fprintf('|%s | | %s | %d | %d | %d | %.1f | %d | %.1e |\n',f,str,...
            nImages,nOP,nIP,raysPerOP,sizeN,sparsity);

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

    Nperm=JTJ(p,p);
    spy(Nperm)
    print('-deps',strrep(f,'.mat','.eps'))
end
asdfa
end

%fix=2;

f=files{fix}

Z=load(fullfile(dataDir,f));

s=Z.s;
e=Z.E;
clear Z

% Permute OP first, then EO, then IO.
    
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

estOP=s.bundle.est.OP;

clear Z s E

profile off
profile on
tic
JTJ=e.final.weighted.J'*e.final.weighted.J;
N=JTJ;

Nperm=N(p,p);
% Perform Cholesky on permuted J'*J.
[LT,fail]=chol(Nperm);

if fail==0
    L=LT';
    % Number of OP parameters
    nOP=nnz(bixOP);
    % Extract blocks of L = [ A, 0; B, C].
    % Diagonal OP block of L
    LA=LT(1:nOP,1:nOP)';
    % Diagonal non-OP block of L
    LC=full(LT(nOP+1:end,nOP+1:end))';
    % Subdiagonal block
    LB=full(LT(1:nOP,nOP+1:end))';
else
    warning(['Posterior covariance matrix was not positive definite. ' ...
             'Results will be inaccurate.'])
    n=size(JTJ,1);
    L=sparse(1:n,1:n,nan,n,n);
    LA=[];
    LB=[];
    LC=[];
end
Lblocks=struct('LA',LA,'LB',LB,'LC',LC);
e.final.factorized=struct('p',p,'L',L,'fail',fail,'Lblocks',Lblocks);
ok=~fail;

CC=VectorizedCOP(Lblocks,estOP);
toc

profile report

asfdd

if size(N,1)<500
    Ninv=inv(full(N));
    NinvOP=zeros(size(CC));
    NinvOP(1:nnz(bixOP),1:nnz(bixOP))=Ninv(bixOP(1):end,bixOP(1):end);
    NinvOP=mkblkdiag(NinvOP,3);
    max(max(abs(NinvOP-CC)))
end

% Method in Wolfe, DeWitt
minOPix=min(bixOP(bixOP>0));
tic
Ndot=N(1:minOPix-1,1:minOPix-1);
Nbar=N(1:minOPix-1,minOPix:end);
Ndot2=N(minOPix:end,minOPix:end);
C1inv=Ndot-Nbar*(Ndot2\Nbar');
LL=chol(Ndot2,'lower');
Ndot2inv=inv(LL')*inv(LL);
C2=Ndot2inv+Ndot2inv*Nbar'*(C1inv\Nbar)*Ndot2inv;
toc
C2block3=mkblkdiag(C2,3);
C2folded=zeros(size(CC));
C2folded(1:nnz(bixOP),1:nnz(bixOP))=C2block3;
max(max(abs(C2folded-CC)))
