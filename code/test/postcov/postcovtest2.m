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

if false
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
end

Nperm(nOP+1:end,1:nOP)=Nperm(nOP+1:end,1:nOP)+eps;
Nperm(1:nOP,nOP+1:end)=Nperm(1:nOP,nOP+1:end)+eps;

wDir=fullfile(getenv('HOME'),'res','reports','work','nice20');

nIm=double(full(Nperm~=0))+1;
[i,j]=find(nIm);

imwrite(nIm,[1,1,1;0.5,0.5,1],fullfile(wDir,'nblocks.png'));

fid=fopen(fullfile(wDir,'nblocks.tex'),'wt+');
fprintf(fid,'\\begin{tikzpicture}[spy using outlines={circle, magnification=10, connect spies}]\n');
fprintf(fid,'\\begin{axis}[y dir=reverse,width=0.45\\hsize,height=0.45\\hsize,axis line style={draw=none},tick style={draw=none},xtick={\\empty},ytick={\\empty}]\n');
fprintf(fid,'\\addplot graphics[xmin=0,xmax=%d,ymin=0,ymax=%d]{%s};\n',size(Nperm),'nblocks.png');
fprintf(fid,'\\draw[ultra thin] (0,0) -- (%d,0) -- (%d,%d) -- (0,%d) -- cycle;\n',(nOP+nEO)*ones(1,4));
fprintf(fid,'\\draw[ultra thick] (0,%d) -- (%d,%d);\n',nOP,nOP+nEO,nOP);
fprintf(fid,'\\draw[ultra thick] (%d,0) -- (%d,%d);\n',nOP,nOP,nOP+nEO);
for i=90:3:99
    for j=0:3
        fprintf(fid,'\\draw[line width=0.02pt] (%d,%d) -- (%d,%d);\n',i,i+j,i+3,i+j);
        fprintf(fid,'\\draw[line width=0.02pt] (%d,%d) -- (%d,%d);\n',i+j,i,i+j,i+3);
    end
end
fprintf(fid,'\\node at (axis cs: %d,%d) {$N_{11}$};\n',round(nOP)/2,round(nOP)/2);
fprintf(fid,'\\node at (axis cs: %d,%d) {$N_{22}$};\n',nOP+round(nEO)/2,nOP+round(nEO)/2);
fprintf(fid,'\\node at (axis cs: %d,%d) {$N_{21}$};\n',round(nOP)/2,nOP+round(nEO)/2);
fprintf(fid,'\\node at (axis cs: %d,%d) {$N_{21}^T$};\n',nOP+round(nEO)/2,round(nOP)/2);
fprintf(fid,'\\coordinate (spypoint) at (axis cs:96,96);\n');
fprintf(fid,'\\coordinate (magnifyglass) at (axis cs:200,56);\n');
fprintf(fid,'\\end{axis}\n');
fprintf(fid,'\\spy[black,size=0.085\\hsize] on (spypoint) in node[fill=white] at (magnifyglass);\n');
fprintf(fid,'\\end{tikzpicture}\n');
fclose(fid);

L=chol(Nperm)';

if false
    figure(2)
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
    set(gca,'xtick',[],'ytick',[],'xlim',xLim,'ylim',xLim)
    xlabel('')
    print('-f2','-dpng',strrep(f,'.mat','-lblocks.png'))
end

lIm=full(L~=0);
imwrite(lIm,[1,1,1;0.5,0.5,1],fullfile(wDir,'lblocks.png'));

fid=fopen(fullfile(wDir,'lblocks.tex'),'wt+');
fprintf(fid,'\\begin{tikzpicture}[spy using outlines={circle, magnification=10, connect spies}]\n');
fprintf(fid,'\\begin{axis}[y dir=reverse,width=0.45\\hsize,height=0.45\\hsize,axis line style={draw=none},tick style={draw=none},xtick={\\empty},ytick={\\empty}]\n');
fprintf(fid,'\\addplot graphics[xmin=0,xmax=%d,ymin=0,ymax=%d]{%s};\n',size(Nperm),'lblocks.png');
fprintf(fid,'\\draw[ultra thin] (0,0) -- (%d,0) -- (%d,%d) -- (0,%d) -- cycle;\n',(nOP+nEO)*ones(1,4));
fprintf(fid,'\\draw[ultra thick] (0,%d) -- (%d,%d);\n',nOP,nOP+nEO,nOP);
fprintf(fid,'\\draw[ultra thick] (%d,0) -- (%d,%d);\n',nOP,nOP,nOP+nEO);
for i=90:3:99
    for j=0:3
        fprintf(fid,'\\draw[line width=0.02pt] (%d,%d) -- (%d,%d);\n',i,i+j,i+min(3,j+1),i+j);
        fprintf(fid,'\\draw[line width=0.02pt] (%d,%d) -- (%d,%d);\n',i+j,i+min(j-1,2),i+j,i+3);
    end
end
fprintf(fid,'\\node at (axis cs: %d,%d) {$L_{11}$};\n',round(nOP)/2,round(nOP)/2);
fprintf(fid,'\\node at (axis cs: %d,%d) {$L_{22}$};\n',nOP+round(nEO)/2,nOP+round(nEO)/2);
fprintf(fid,'\\node at (axis cs: %d,%d) {$L_{21}$};\n',round(nOP)/2,nOP+round(nEO)/2);
fprintf(fid,'\\node at (axis cs: %d,%d) {$0$};\n',nOP+round(nEO)/2,round(nOP)/2);
fprintf(fid,'\\coordinate (spypoint) at (axis cs:96,96);\n');
fprintf(fid,'\\coordinate (magnifyglass) at (axis cs:200,56);\n');
fprintf(fid,'\\end{axis}\n');
fprintf(fid,'\\spy[black,size=0.085\\hsize] on (spypoint) in node[fill=white] at (magnifyglass);\n');
fprintf(fid,'\\end{tikzpicture}\n');
fclose(fid);

fid=fopen(fullfile(wDir,'ublocks.tex'),'wt+');
fprintf(fid,'\\begin{tikzpicture}[spy using outlines={circle, magnification=10, connect spies}]\n');
fprintf(fid,'\\begin{axis}[y dir=reverse,width=0.45\\hsize,height=0.45\\hsize,axis line style={draw=none},tick style={draw=none},xtick={\\empty},ytick={\\empty}]\n');
fprintf(fid,'\\addplot graphics[xmin=0,xmax=%d,ymin=0,ymax=%d]{%s};\n',size(Nperm),'lblocks.png');
fprintf(fid,'\\draw[ultra thin] (0,0) -- (%d,0) -- (%d,%d) -- (0,%d) -- cycle;\n',(nOP+nEO)*ones(1,4));
fprintf(fid,'\\draw[ultra thick] (0,%d) -- (%d,%d);\n',nOP,nOP+nEO,nOP);
fprintf(fid,'\\draw[ultra thick] (%d,0) -- (%d,%d);\n',nOP,nOP,nOP+nEO);
for i=90:3:99
    for j=0:3
        fprintf(fid,'\\draw[line width=0.02pt] (%d,%d) -- (%d,%d);\n',i,i+j,i+min(3,j+1),i+j);
        fprintf(fid,'\\draw[line width=0.02pt] (%d,%d) -- (%d,%d);\n',i+j,i+min(j-1,2),i+j,i+3);
    end
end
fprintf(fid,'\\node at (axis cs: %d,%d) {$U_{11}$};\n',round(nOP)/2,round(nOP)/2);
fprintf(fid,'\\node at (axis cs: %d,%d) {$U_{22}$};\n',nOP+round(nEO)/2,nOP+round(nEO)/2);
fprintf(fid,'\\node at (axis cs: %d,%d) {$U_{21}$};\n',round(nOP)/2,nOP+round(nEO)/2);
fprintf(fid,'\\node at (axis cs: %d,%d) {$0$};\n',nOP+round(nEO)/2,round(nOP)/2);
fprintf(fid,'\\coordinate (spypoint) at (axis cs:96,96);\n');
fprintf(fid,'\\coordinate (magnifyglass) at (axis cs:200,56);\n');
fprintf(fid,'\\end{axis}\n');
fprintf(fid,'\\spy[black,size=0.085\\hsize] on (spypoint) in node[fill=white] at (magnifyglass);\n');
fprintf(fid,'\\end{tikzpicture}\n');
fclose(fid);

if false
    figure(3)
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
    set(gca,'xtick',[],'ytick',[],'xlim',xLim,'ylim',xLim)
    xlabel('')
    print('-f3','-dpng',strrep(f,'.mat','-ublocks.png'))
    
end
   
%print('-f1','-depsc',strrep(f,'.mat','-allblocks.eps'))

 
