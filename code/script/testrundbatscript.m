testCase
switch testCase
  case 1
    srcDir=fullfile(fileparts(dbatroot),'/data/script/camcaldemo');

    f=fullfile(srcDir,'camcaldemo.xml');
  case 2
    srcDir=fullfile(fileparts(dbatroot),'/data/script/romabundledemo');

    f=fullfile(srcDir,'romabundledemo.xml');
  case 3
    srcDir=fullfile(fileparts(dbatroot),'/data/script/sxb');

    f=fullfile(srcDir,'sxb.xml');
  otherwise
    error('Bad test case');
end

[s,xml]=rundbatscript(f,true);
