function d=unpackpsz(psFile,unpackDir,unpackAscii)
%UNPACKPSZ Unpack .PSZ file.
%
%   UNPACKPSZ(PSFILE,UNPACKDIR) unpacks the .PSZ file PSFILE into
%   the directory UNPACKDIR. UNPACKDIR will be created if it does
%   not exist.
%
%   UNPACKPSZ(PSFILE,UNPACKDIR,TRUE) will furthermore convert each
%   .PLY file into its ascii equivalent in the directory
%   FULLFILE(UNPACKDIR,'ascii').
%
%   DIRS=... returns what directories were created or unpacked
%   into. The returned dirs are sorted such that the last directory
%   should be removed first.

%See also: UNZIP, PLYTOASCII.

% Create unpacked dir if necessary.
if ~exist(unpackDir)
    mkdir(unpackDir)
end

d={unpackDir};

if unpackAscii
    % Delete any old files, including ascii versions of .ply files.
    asciiDir=fullfile(unpackDir,'ascii');

    d{end+1}=asciiDir;
    if exist(asciiDir)
        delete(fullfile(asciiDir,'*'));
        rmdir(asciiDir);
    end
end

delete(fullfile(unpackDir,'*'));

% Unzip the .PSZ file.
unzip(psFile,unpackDir)

if unpackAscii
    % Unpack any binary .PLY files.
    plytoascii(unpackDir)
end
