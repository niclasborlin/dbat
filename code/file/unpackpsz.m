function unpackpsz(psFile)
%UNPACKPSZ Unpack .PSZ file.
%
%   UNPACKPSZ(PSFILE) unpacks the .PSZ file PSFILE into a
%   subdirectory called 'unpacked'.

psDir=fileparts(psFile);

psUnpackedDir=fullfile(psDir,'unpacked');

% Create unpacked dir if necessary.
if ~exist(psUnpackedDir)
    mkdir(psUnpackedDir)
end

% Delete any old files, including ascii versions of .ply files.
asciiDir=fullfile(psUnpackedDir,'ascii');

if exist(fullfile(psUnpackedDir,'ascii'),'dir')
    delete(fullfile(asciiDir,'*'));
    rmdir(asciiDir);
end
delete(fullfile(psUnpackedDir,'*'));

% Unzip the .PSZ file.
unzip(psFile,psUnpackedDir)

% Unpack any binary .PLY files.
plytoascii(psUnpackedDir)
