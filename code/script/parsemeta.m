function [name,projUnit]=parsemeta(doc)
%PARSEMETA Parse meta section of a DBAT XML file.
%
%   [NAME,PROJUNIT]=PARSEMETA(DOC) parses the name and project_unit
%   fields of the meta subfield of the document DOC of a DBAT XML
%   file. Missing fields defaults to the blank string.
%
%See also: SETPROJINFO, PARSEINPUT.

narginchk(1,1)

if isfield(doc,'meta') && isfield(doc.meta,'name') && ...
        isfield(doc.meta.name,'Text')
    name=doc.meta.name.Text;
else
    name='';
end

if isfield(doc,'meta') && isfield(doc.meta,'project_unit') && ...
        isfield(doc.meta.name,'Text')
    projUnit=doc.meta.project_unit.Text;
else
    projUnit='';
end
