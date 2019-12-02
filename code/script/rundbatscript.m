function [s,xml]=rundbatscript(f,verbose)
%RUNDBATSCRIPT Run a DBAT XML script
%
%   S=RUNDBATSCRIPT(XMLFILE) loads and executes the DBAT script stored
%   in the xml file XMLFILE. The result is returned in the DBAT
%   structure S. Other files and/or plots may be generated, depending
%   on the content of the script file. See the manual for a further
%   description of the script language.
%
%   Use S=RUNDBATSCRIPT(XMLFILE,TRUE) to turn on verbose output.

if nargin<2, verbose=false; end

if nargin<1
    % Ask for main xml file.
    [f,p]=uigetfile(fullfile(fileparts(dbatroot),'data','script','*.xml'));
    f=fullfile(p,f);
end

if verbose, fprintf('Loading xml file %s...',f); end
xml=dbatxml2struct(f);
if verbose, fprintf('done.\n'); end

% Top level element is just 'document'.
[ok,msg]=checkxmlfields(xml,'document');
if ~ok
    allFields=join(fieldnames(xml),', ');
    error('DBAT XML top level error: %s. Read fields are: %s.', ...
          msg,allFields{1});
end

doc=xml.document;

% Check top level fields. The field 'meta' is optional.
docFields={'input','output','operations','Attributes','meta','c'};
[ok,msg]=checkxmlfields(doc,docFields,[true,true,true,true,false,false]);
if ~ok
    allFields=join(fieldnames(doc),', ');
    error('DBAT XML document field error: %s. Read fields are: %s',...
          msg,allFields{1});
end

% Check that DBAT script version is known.
firstKnownVersion='1.0';
lastKnownVersion='1.0';
[ok,msg]=checkversionattr(doc.Attributes,'dbat_script_version', ...
                          lastKnownVersion,firstKnownVersion);
if ~ok, error('DBAT script error: %s',msg); end

% Parse the meta block.
[name,projUnit]=parsemeta(doc);

% Parse input.
[s,imDir,cptFile,EOfile]=parseinput(doc.input,f);

% Set project info in the DBAT struct.
s=setprojinfo(s,f,name,projUnit,imDir,cptFile,EOfile);

% Parse and execute the loaded operations.
s=parseops(s,doc.operations);

% Parse and execute the loaded output operations.
s=parseoutput(s,doc.output,f);
