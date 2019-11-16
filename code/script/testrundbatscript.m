camcaldir=fullfile(fileparts(dbatroot),'/data/script/camcaldemo');

f=fullfile(camcaldir,'camcaldemo.xml');

[s,xml]=rundbatscript(f,true);
