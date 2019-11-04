function [dir,rawDir]=getattrbasedir(s,docFile,attrName)
%GETBASEDIR Get base_dir attribute if present.
%
%   GETATTRBASEDIR(S,DOCFILE), where S is a sub-structure loaded from
%   the DBAT XML file DOCFILE, returns the base_dir XML attribute
%   (S.Attributes.base_dir) if present, otherwise the blank string.
%
%   If base_dir starts with any of the special strings $DBAT, $HOME,
%   or $HERE, it is replaced by the DBAT installation directory, the
%   user home directory, or the directory containing DOCFILE,
%   respectively.
% 
%   [DIR,RAWDIR]=GETATTRBASEDIR(S,DOCFILE) also returns the
%   unprocessed base_dir string in RAWDIR.
%
%   Use GETATTRBASEDIR(S,DOCFILE,ATTRNAME) to specify another
%   attribute name ATTRNAME than 'base_dir'.
%
%   NOTE: The user home directory is returned by
%   java.lang.System.getProperty('user.home') and is typically
%   C:\Users\username on Windows, /Users/username on Mac/OS X, and
%   /home/username on Linux, where 'username' is the current user name.
    
narginchk(2,3)

if nargin<3, attrName='base_dir'; end

if isfield(s,'Attributes') && isfield(s.Attributes,attrName)
    rawDir=s.Attributes.(attrName);
    dir=rawDir;
    % Check for  special strings.
    if isprefixed(dir,'$DBAT')
        % The DBAT installation directory is one level up from the
        % DBAT code root.
        root=fileparts(dbatroot);
        dir=fullfile(root,dir(6:end));
    elseif isprefixed(dir,'$HOME')
        home=char(java.lang.System.getProperty('user.home'));
        dir=fullfile(home,dir(6:end));
    elseif isprefixed(dir,'$HERE')
        docDir=fileparts(docFile);
        dir=fullfile(docDir,dir(6:end));
    end
else
    rawDir='';
    dir='';
end

function tf=isprefixed(s,dir)
%Return TRUE if s has dir as its first path component.

if length(s)<length(dir)
    tf=false;
else
    tf=strcmp(s,dir) || startsWith(s,[dir,'/']) || startsWith(s,[dir,'\']);
end

