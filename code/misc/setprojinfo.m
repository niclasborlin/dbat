function s=setprojinfo(s,fileName,title,objUnit,imDir,cptFile,EOfile,x0desc)
%SETPROJINFO Set project info in DBAT struct
%
%   S=SETPROJINFO(S,FILENAME,TITLE,OBJUNIT) sets the project
%   filename, title, and object unit in the DBAT struct S.
%
%   S=SETPROJINFO(S,FILENAME,TITLE,OBJUNIT,IMDIR,CPTFILE,EOFILE,X0DESC)
%   sets all project info fields in the DBAT struct S to the supplied
%   parameters. All parameters from IMDIR onwards default to the blank
%   string.
    
narginchk(4,8)

if nargin<5, imDir=''; end
if nargin<6, cptFile=''; end
if nargin<7, EOfile=''; end
if nargin<8, x0desc=''; end

s.proj.fileName=fileName;
s.proj.title=title;
s.proj.objUnit=objUnit;
s.proj.imDir=imDir;
s.proj.cptFile=cptFile;
s.proj.EOfile=EOfile;
s.proj.x0desc=x0desc;
