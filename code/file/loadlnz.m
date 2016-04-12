function s=loadlnz(fName)
%LOADLNZ Load Photoscan LNZ camera calibration file.
%

psDir=fileparts(psFile);
unpackpsz(psFile);
psUnpackedDir=fullfile(psDir,'unpacked');
fName=fullfile(psUnpackedDir,'doc.xml');
s=xml2struct2(fName);
