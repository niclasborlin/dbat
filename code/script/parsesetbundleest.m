function s=parsesetbundleest(s,xml)
%PARSESETBUNDLEEST Set what parameters the bundle should estimate.
%
%   S=PARSESETBUNDLEEST(S,XML) parses the SET_BUNDLE_ESTIMATE_PARAMS
%   block of a DBAT script operation and sets the parameters to
%   estimate in the DBAT structure S.

narginchk(2,2);

knownFields={'io','eo','op'};
[ok,msg]=checkxmlfields(xml,knownFields,false(size(knownFields)));
if ~ok, error('DBAT XML script set_initial_values error: %s',msg); end

if isfield(xml,'io')
    %    s=parsesetbundleestio(s,xml.io);
end

if isfield(xml,'eo')
    s=parsesetbundleesteo(s,xml.eo);
end

if isfield(xml,'op')
    s=parsesetbundleestop(s,xml.op);
end

