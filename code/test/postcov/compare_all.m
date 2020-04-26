dispProp=true;

addpath('../sparseinv','-end');

%dataDir=fullfile(fileparts(dbatroot),'data','test');
dataDir=fullfile('/scratch','niclas','dbat_data');

files={'camcaldemo.mat',...
       'stpierre.mat',...
       'romabundledemo-selfcal.mat',...
       'vexcel.mat',...
       'mit-2d-3ray.mat',...
       'sewu-filt35.mat',...
       'mit-3d-xhatch2-3ray.mat',...
       'mit-3d-xhatch-3ray.mat',...
       'mit-2lyr-3ray.mat',...
       'mit-3lyr-3ray.mat'};

%'mit-2d-strip-and-border-3ray.mat',...
      %  'sewu-filt25.mat',...
      %  'mit-2d-2ray.mat',...
      %  'mit-2d-strip-and-border-2ray.mat',...
      %  'mit-3d-xhatch-2ray.mat',...
      %  'mit-3lyr-2ray.mat',...
      % };

timeTable=zeros(length(files),5);
timeTables=cell(length(files),5);

abserr=@(A,B)norm(A-B,'fro');
relerr=@(A,B)abserr(A,B)/norm(A,'fro');

sparsity=@(A)nnz(A)/numel(A);

for fix=1:length(files)

f=files{fix}
Z=load(fullfile(dataDir,f));
s=Z.s;
e=Z.E;

spVis=sparsity(s.IP.vis);

JTJ=e.final.weighted.J'*e.final.weighted.J;
N=JTJ;

clear Z e

% IO blocks.
[i,j]=ind2sub(size(s.bundle.est.IO),s.bundle.serial.IO.src);
bixIO=full(sparse(i,j,s.bundle.serial.IO.dest,size(s.bundle.est.IO,1),size(s.bundle.est.IO,2)));
bixIO=bixIO(bixIO~=0);

% EO blocks.
[i,j]=ind2sub(size(s.bundle.est.EO),s.bundle.serial.EO.src);
bixEO=full(sparse(i,j,s.bundle.serial.EO.dest,size(s.bundle.est.EO,1),size(s.bundle.est.EO,2)));
bixEO=bixEO(bixEO~=0);

% OP blocks.
[i,j]=ind2sub(size(s.bundle.est.OP),s.bundle.serial.OP.src);
bixOP=full(sparse(i,j,s.bundle.serial.OP.dest,size(s.bundle.est.OP,1),size(s.bundle.est.OP,2)));
bixOP=bixOP(bixOP~=0);

% Macro to extract blocks of matrix. Blocksizes are M followed by N.
Part11=@(A,m,n)A(1:m,1:m);
Part12=@(A,m,n)A(1:m,m+1:m+n);
Part22=@(A,m,n)A(m+1:m+n,m+1:m+n);

% Re-partition to swap blocks of size M followed by N.
Swap=@(A,m,n)A([m+1:m+n,1:m],[m+1:m+n,1:m]);

nIO=nnz(bixIO);
nOP=nnz(bixOP);
nEO=nnz(bixEO);

numCams=size(s.EO.val,2);
numOPs=size(s.OP.val,2);
numIPs=size(s.IP.val,2);

clear s

% Remove IO parameters
NnoIO=N(nIO+1:end,nIO+1:end);

sp12=sparsity(Part12(NnoIO,nEO,nOP));

if dispProp
    if fix==1
        fprintf('| Name | SC | Cams | OP | IP | sp(B) |\n');
        fprintf('|-\n');
    end
    
    % Name, nIO, nEO, nOP, nIP, sparsity(B)
    fprintf('| %s | %d | %d | %d | %d | %.1f |\n',f,nIO>0,numCams,numOPs,...
            numIPs,sp12*100);
end

fprintf('.');
% Compute post OP cov with classic algorithm.
[classic,C1]=time_classic(NnoIO,0,nEO,nOP,true,false);

[classicSelf,C1S]=time_classic(N,0,nIO+nEO,nOP,true,false);

if 0
    fprintf('.');
    % Compute post OP cov with SI algorithm on IO-EO-OP permutation.
    [siIOfirst,C2]=time_si(NnoIO,0,nEO,nOP,false);
    
    fprintf('.');
    % Compute post OP cov with SI algorithm on OP-EO-IO permutation.
    [siIOlast,C3]=time_si(NnoIO,0,nEO,nOP,true);
else
    siIOfirst=0;
    C2=C1;
    siIOlast=0;
    C3=C1;
end
fprintf('.');
% Compute post OP cov with CIP algorithm on OP-EO-IO permutation.
[icip,C4,spLB,spLC]=time_icip_dense(NnoIO,0,nEO,nOP,false);

sp12*100
spLB*100
spLC*100

fprintf('.');
% Compute post OP cov with CIP algorithm on OP-EO-IO permutation.
%[icipSparse,C4sparse]=time_icip_sparse(NnoIO,0,nEO,nOP,false);
icipSparse=0;
C4sparse=C4;

fprintf('.');
% Compute post OP cov with classic algorithm.
[classicDiag,C1d]=time_classic(NnoIO,0,nEO,nOP,true,true);

fprintf('.');
% ICIP on only the diagonal.
[icipDiag,C4d]=time_icip_dense(NnoIO,0,nEO,nOP,true);

fprintf('\n');

if fix==1
    % Verify result.
    Ni=inv(N);
    C0=Ni(nIO+nEO+1:end,nIO+nEO+1:end);
    C0=mkblkdiag(C0,3);

    % Use no-IO version for classic
    NinoIO=inv(NnoIO);
    C0noIO=NinoIO(nEO+1:end,nEO+1:end);
    C0noIO=mkblkdiag(C0noIO,3);
    C0d=diag(diag(C0noIO));
    
    disp('Classic errors:')
    disp([abserr(C0noIO,C1),relerr(C0noIO,C1)])
    
    disp('SI-IO-EO-OP errors:')
    disp([abserr(C0noIO,C2),relerr(C0noIO,C2)])
    
    disp('SI-OP-EO-IO errors:')
    disp([abserr(C0noIO,C3),relerr(C0noIO,C3)])
    
    disp('ICIP errors:')
    disp([abserr(C0noIO,C4),relerr(C0noIO,C4)])

    disp('ICIP-sparse errors:')
    disp([abserr(C0noIO,C4sparse),relerr(C0noIO,C4sparse)])

    disp('classic-diag errors:')
    disp([abserr(C0d,C1d),relerr(C0d,C1d)])

    disp('ICIP-diag errors:')
    disp([abserr(C0d,C4d),relerr(C0d,C4d)])
end

timeTable(fix,1)=classic(end);
timeTable(fix,2)=siIOfirst(end);
timeTable(fix,2)=classicSelf(end);
timeTable(fix,3)=siIOlast(end);
timeTable(fix,4)=icip(end);
timeTable(fix,5)=icipSparse(end);
timeTable(fix,6)=classicDiag(end);
timeTable(fix,7)=icipDiag(end);
timeTable(fix,[8:11])=[spVis,sp12,spLB,spLC];

timeTable

timeTables{fix,1}=classic;
timeTables{fix,2}=siIOfirst;
timeTables{fix,2}=classicSelf;
timeTables{fix,3}=siIOlast;
timeTables{fix,4}={icip,spVis,sp12,spLB,spLC};
timeTables{fix,5}=icipSparse;
timeTables{fix,6}=classicDiag;
timeTables{fix,7}=icipDiag;

end