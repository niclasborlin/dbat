function plytoascii(plyDir)
%PLYTOASCII Unpack .PLY files to ascii.
%
%   PLYTOASCII(PLYDIR) converts each *.PLY file in PLYDIR into its
%   ASCII version in the subdirectory 'ascii'.

asciiDir=fullfile(plyDir,'ascii');
if ~exist(asciiDir)
    mkdir(asciiDir);
end
z=dir(fullfile(plyDir,'*.ply'));
for i=1:length(z)
    %z(i).name
    [~,~,d,~]=ply_read(fullfile(plyDir,z(i).name),'tri');
    ply_write(d,fullfile(asciiDir,z(i).name),'ascii','double');
end

