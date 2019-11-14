function s=parsesetinitialvalues(s,xml)
%PARSESETINITIALVALUES Parse and set initial values for a bundle operation.
%
%   S=PARSESETINITIALVALUES(S,XML) parses the SET_INITIAL_VALUES block
%   of a DBAT script operation. The DBAT structure S is updated
%   with the new initial values.

narginchk(2,2);

knownFields={'io','eo','op'};
[ok,msg]=checkxmlfields(xml,knownFields,false(size(knownFields)));
if ~ok, error('DBAT XML script set_initial_values error: %s',msg); end

if isfield(xml,'io')
    s=parsesetinitialiovalues(s,xml.io);
end

if isfield(xml,'eo')
    s=parsesetinitialeovalues(s,xml.eo);
end

if isfield(xml,'op')
    s=parsesetinitialopvalues(s,xml.op);
end

