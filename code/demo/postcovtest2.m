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

Nperm=JTJ(p,p);
figure(1)
subplot(1,3,1)
spy(Nperm)
lw=4;
pOP=nnz(bixOP);
h1=line([0,sizeN],pOP*[1,1],'color','k','linewidth',lw);
h2=line(pOP*[1,1],[0,sizeN],'color','k','linewidth',lw);
h11=text(pOP*0.5,pOP*0.5,'N_{11}',...
         'horizontalalignment','center','fontsize',24,'fontangle','italic');
h21=text(pOP*0.5,(pOP+sizeN)*0.5,'N_{21}',...
         'horizontalalignment','center','fontsize',24,'fontangle','italic');
h22=text((pOP+sizeN)*0.5,(pOP+sizeN)*0.5,'N_{22}',...
         'horizontalalignment','center','fontsize',24,'fontangle','italic');
set(gca,'xtick',[],'ytick',[])
xlabel('')
%print('-f1','-depsc',strrep(f,'.mat','-nblocks.eps'))

figure(1)
subplot(1,3,2)
L=chol(Nperm)';
spy(L)
lw=4;
pOP=nnz(bixOP);
h1=line([0,sizeN],pOP*[1,1],'color','k','linewidth',lw);
h2=line(pOP*[1,1],[0,sizeN],'color','k','linewidth',lw);
h11=text(pOP*0.5,pOP*0.5,'L_{11}',...
         'horizontalalignment','center','fontsize',24,'fontangle','italic');
h21=text(pOP*0.5,(pOP+sizeN)*0.5,'L_{21}',...
         'horizontalalignment','center','fontsize',24,'fontangle','italic');
h22=text((pOP+sizeN)*0.5,(pOP+sizeN)*0.5,'L_{22}',...
         'horizontalalignment','center','fontsize',24,'fontangle','italic');
set(gca,'xtick',[],'ytick',[])
xlabel('')
%print('-f2','-depsc',strrep(f,'.mat','-lblocks.eps'))

figure(1)
subplot(1,3,3)
U=inv(L);
spy(U)
lw=4;
pOP=nnz(bixOP);
h1=line([0,sizeN],pOP*[1,1],'color','k','linewidth',lw);
h2=line(pOP*[1,1],[0,sizeN],'color','k','linewidth',lw);
h11=text(pOP*0.5,pOP*0.5,'U_{11}',...
         'horizontalalignment','center','fontsize',24,'fontangle','italic');
h21=text(pOP*0.5,(pOP+sizeN)*0.5,'U_{21}',...
         'horizontalalignment','center','fontsize',24,'fontangle','italic');
h22=text((pOP+sizeN)*0.5,(pOP+sizeN)*0.5,'U_{22}',...
         'horizontalalignment','center','fontsize',24,'fontangle','italic');
set(gca,'xtick',[],'ytick',[])
xlabel('')
%print('-f3','-depsc',strrep(f,'.mat','-ublocks.eps'))

print('-f1','-depsc',strrep(f,'.mat','-allblocks.eps'))
