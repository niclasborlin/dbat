dataDir=fullfile('/scratch','niclas','dbat_data');

files={'camcaldemo.mat',...
       'stpierre.mat',...
       'romabundledemo-selfcal.mat',...
       'romabundledemo-imagevariant.mat',...
       'vexcel.mat',...
       'mit-2d-3ray.mat',...
       'sewu-filt35.mat',...
       'mit-3d-xhatch2-3ray.mat',...
       'mit-3d-xhatch-3ray.mat',...
       'mit-2lyr-3ray.mat',...
       'mit-3lyr-3ray.mat'};

selfCal=true;

fix=3;
        
f=files{fix}
Z=load(fullfile(dataDir,f));
s=Z.s;
e=Z.E;

sparsity=@(A)nnz(A)/numel(A);

spVis=sparsity(s.IP.vis);

JTJ=e.final.weighted.J'*e.final.weighted.J;
N=JTJ;
        
clear Z

% IO blocks.
[i,j]=ind2sub(size(s.bundle.est.IO),s.bundle.serial.IO.src);
bixIO=full(sparse(i,j,s.bundle.serial.IO.dest,size(s.bundle.est.IO,1),...
                  size(s.bundle.est.IO,2)));
bixIO=bixIO(bixIO~=0);

% EO blocks.
[i,j]=ind2sub(size(s.bundle.est.EO),s.bundle.serial.EO.src);
bixEO=full(sparse(i,j,s.bundle.serial.EO.dest,size(s.bundle.est.EO,1),...
                  size(s.bundle.est.EO,2)));
bixEO=bixEO(bixEO~=0);

% OP blocks.
[i,j]=ind2sub(size(s.bundle.est.OP),s.bundle.serial.OP.src);
bixOP=full(sparse(i,j,s.bundle.serial.OP.dest,size(s.bundle.est.OP,1),...
                  size(s.bundle.est.OP,2)));
bixOP=bixOP(bixOP~=0);
        
% Macros to extract blocks of matrix. Blocksizes are M followed by N
% followed by K.
Blk11=@(A,m,n,k)A(1:m,1:m);
Blk12=@(A,m,n,k)A(1:m,m+1:m+n);
Blk13=@(A,m,n,k)A(1:m,m+n+1:m+n+k);
Blk22=@(A,m,n,k)A(m+1:m+n,m+1:m+n);
Blk23=@(A,m,n,k)A(m+1:m+n,m+n+1:m+n+k);
Blk33=@(A,m,n,k)A(m+n+1:m+n+k,m+n+1:m+n+k);

nIO=nnz(bixIO);
nOP=nnz(bixOP);
nEO=nnz(bixEO);

numCams=size(s.EO.val,2);
numOPs=size(s.OP.val,2);
numIPs=size(s.IP.val,2);

sp12=sparsity(Blk23(N,nIO,nEO,nOP));
    
if selfCal
    Nclassic=N;
else
    % Remove IO parameters
    Nclassic=N(nIO+1:end,nIO+1:end);
    nIO=0;
end        
ix=[nIO+nEO+1:size(Nclassic,2),nIO+1:nIO+nEO,1:nIO];
Nperm=Nclassic(ix,ix);

NN=full(N)~=0;
        
imshow(NN,[255,255,255;204,204,255]/255);
line([0,0],[0,size(NN,1)],'color','k')
line([0,0]+nIO,[0,size(NN,1)],'color','k')
line([0,0]+nIO+nEO,[0,size(NN,1)],'color','k')
line([0,0]+nIO+nEO+nOP,[0,size(NN,1)],'color','k')
line([0,size(NN,1)],[0,0],'color','k')
line([0,size(NN,1)],[0,0]+nIO,'color','k')
line([0,size(NN,1)],[0,0]+nIO+nEO,'color','k')
line([0,size(NN,1)],[0,0]+nIO+nEO+nOP,'color','k')

%matlab2tikz('/tmp/spyfig.tex','imagesAsPng',true);
