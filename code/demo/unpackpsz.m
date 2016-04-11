function unpackpsz(psFile)

psDir=fileparts(psFile);
psUnpackedDir=fullfile(psDir,'unpacked');

unzip(psFile,psUnpackedDir);

plytoascii(psUnpackedDir);
