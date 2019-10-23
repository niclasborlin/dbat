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

if verbose, fprintf('Loading xml file %s...',f); end
xml=dbatxml2struct(f);
if verbose, fprintf('done.\n'); end

s=[];

% Top level element is just 'document'.
[ok,msg]=checkxmlfields(xml,'document');

if ~ok
    error('DBAT XML top level error: %s',msg);
end

doc=xml.document;

% Document fields are input, output, operations, and Attributes.
docFields={'input','output','operations','Attributes'};

[ok,msg]=checkxmlfields(doc,docFields);
if ~ok
    error('DBAT XML document field error: %s',msg);
end

% Known DBAT script version interval.
firstKnownVersion='1.0';
lastKnownVersion='1.0';

[ok,msg]=checkversionattr(doc.Attributes,'dbat_script_version', ...
                          lastKnownVersion,firstKnownVersion);

if ~ok
    error('DBAT script error: %s',msg);
end

input=doc.input;

% Do stuff with input...

