function unpackpsz(psFile,keepOld)
%UNPACKPSZ Unpack Photoscan .psz archive
%
%   UNPACKPSZ(PSFILE) unpacks the Photoscan zip (.psz) archive in
%   PSFILE into a subdirectory called unpacked. Furthermore, any
%   unpacked .ply files are converted to ascii into
%   unpacked/ascii. All existing files in those directories are
%   deleted. Use UNPACKPSZ(PSFILE,TRUE) to keep old files. NOTE:
%   Unpacked files are always overwritten by the corresponding
%   files in the archive.

if nargin<2, keepOld=false; end

psDir=fileparts(psFile);
psUnpackedDir=fullfile(psDir,'unpacked');

if ~exist(psUnpackedDir)
    mkdir(psUnpackedDir)
end

if ~keepOld
    asciiDir=fullfile(psUnpackedDir,'ascii');
    if exist(asciiDir)
        rmdir(asciiDir,'s')
    end
    delete(fullfile(psUnpackedDir,'*'))
end

unzip(psFile,psUnpackedDir);

plytoascii(psUnpackedDir);
