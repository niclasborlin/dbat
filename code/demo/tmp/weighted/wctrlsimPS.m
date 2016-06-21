fprintf('<property name="fixed" value="true"/>\n');
fprintf('<calibration type="frame" class="adjusted">\n');
fprintf('<resolution width="2272" height="1704"/>\n');
fprintf('<fx>%.15e</fx>\n<fy>%.15e</fy>\n',abs(K(1,1)),abs(K(1,1)));
fprintf('<cx>%.15e</cx>\n<cy>%.15e</cy>\n',abs(K(1:2,3)));
for i=1:3
    fprintf('<k%d>0</k%d>\n',i,i);
end
fprintf('</calibration>\n');

fprintf('\n');

for OPid=1001:1004
    OPid
    pts=ss.markPts(:,ss.colPos(ss.OPid==OPid,:));
    i=0:size(pts,2)-1;
    fprintf('<location camera_id="%d" pinned="true" x="%.7e" y="%.7e"/>\n',...
            [i;pts]);
    fprintf('\n');
end

pr=cell(1,20);
for i=1:length(pr)
    fName=fullfile(getenv('HOME'),'photoscan','cptest2',...
                   sprintf('projections%d.ply',i-1));
    [~,~,pr{i}]=ply_read(fName,'tri');
end

fName=fullfile(getenv('HOME'),'photoscan','cptest2','points0.ply');
[~,~,p0]=ply_read(fName,'tri');

fName=fullfile(getenv('HOME'),'photoscan','cptest2','tracks0.ply');
[~,~,t0]=ply_read(fName,'tri');

id0=0:nnz(ss.OPid<1000)-1;

pp0=p0;
pp0.vertex.x=reshape(ss.OP(1,ss.OPid<1000),[],1);
pp0.vertex.y=reshape(ss.OP(2,ss.OPid<1000),[],1);
pp0.vertex.z=reshape(ss.OP(3,ss.OPid<1000),[],1);
pp0.vertex.id=reshape(id0,[],1);

ppr=pr;
for i=1:length(ppr)
    cols=ss.colPos(ss.vis(:,i) & ss.OPid<1000,i);
    id=id0(ss.vis(:,i) & ss.OPid<1000);
    ppr{i}.vertex.x=reshape(ss.markPts(1,cols),[],1);
    ppr{i}.vertex.y=reshape(ss.markPts(2,cols),[],1);
    ppr{i}.vertex.id=reshape(id,[],1);
end

tt0=t0;
tt0.vertex.red=zeros(length(id0),1);
tt0.vertex.green=zeros(length(id0),1);
tt0.vertex.blue=zeros(length(id0),1);

fName=fullfile(getenv('HOME'),'photoscan','cptest2','modified','points0.ply');
ply_write(pp0,fName,'binary_little_endian')

fName=fullfile(getenv('HOME'),'photoscan','cptest2','modified','tracks0.ply');
ply_write(tt0,fName,'binary_little_endian')

for i=1:length(ppr)
    fName=fullfile(getenv('HOME'),'photoscan','cptest2','modified',...
                   sprintf('projections%d.ply',i-1));
    ply_write(ppr{i},fName,'binary_little_endian')
end
