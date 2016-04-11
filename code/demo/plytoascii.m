function plytoascii(plyDir)

asciiDir=fullfile(plyDir,'ascii');
if ~exist(asciiDir)
    mkdir(asciiDir);
end
z=dir(fullfile(plyDir,'*.ply'));
for i=1:length(z)
    %z(i).name
    [~,~,d,~]=ply_read(fullfile(plyDir,z(i).name),'tri');
    ply_write(d,fullfile(asciiDir,z(i).name),'ascii');
end

