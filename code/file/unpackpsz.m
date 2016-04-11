function unpackpsz(psFile)
%UNPACKPSZ Unpack .psz file into /unpacked subdir.
%
%   UNPACKPSZ(PSFILE)

psDir=fileparts(psFile);

psUnpackedDir=fullfile(psDir,'unpacked');

if ~exist(psUnpackedDir)
    mkdir(psUnpackedDir)
end

delete(fullfile(psUnpackedDir,'*'));
delete(fullfile(psUnpackedDir,'ascii','*'));

unzip(psFile,psUnpackedDir)

plytoascii(psUnpackedDir)
