camcaldir=fullfile(dbatroot,'/demo/data/script/camcaldemo');

f=fullfile(camcaldir,'camcaldemo.xml');

[s,xml]=rundbatscript(f,true);

%f1='/home/niclas/res/reports/work/lc3d-script/camcaldemo/camcaldemo.xml';
%[s,xml]=rundbatscript(f1,true);
%f2='/home/niclas/res/reports/work/lc3d-script/camcaldemo/camcaldemo2.xml';

c1=fullfile(camcaldir,'reference/camcal-fixed.txt');
c2=fullfile(camcaldir,'reference/camcal-weighted.txt');

pts1=loadctrlpts(c1,'id,label,x,y,z');
pts2=loadctrlpts(c2,'id,label,x,y,z,sx,sy,sz');
pts3=loadctrlpts(c2,'id,label,x,y,z,sxy,dummy,sz');
pts4=loadctrlpts(c2,'id,label,x,y,z,sxyz,dummy,dummy');

i1=fullfile(camcaldir,'images/images.txt');
ims1=loadimagetable(i1,'id,label,path');
ims2=loadimagetable(i1,'id,dummy,path');

imp1=fullfile(camcaldir,'measurements/markpts.txt');
mark1=loadimagepts(imp1,'im,id,x,y,sxy');

%cam1='/home/niclas/res/reports/work/lc3d-script/romabundledemo/cameras/c4040z.xml';
%[cam,camXml]=loadcameras(cam1);
