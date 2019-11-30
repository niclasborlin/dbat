function s=parsesetdatum(s,xml)
%PARSESETDATUM Set what parameters the bundle should estimate.
%
%   S=PARSESETDATUM(S,XML) parses the set_datum block of a DBAT script
%   operation.

narginchk(2,2);

knownFields={'Text','Attributes'};
[ok,msg]=checkxmlfields(xml,knownFields,true(size(knownFields)));
if ~ok, error('DBAT XML script set_datum error: %s',msg); end

switch strip(xml.Text)
  case 'depend'
    knownFields={'ref_cam','ref_base'};
    [ok,msg]=checkxmlfields(xml.Attributes,knownFields);
    if ~ok, error('DBAT XML script set_datum/Attributes error: %s',msg); end

    refCam=sscanf(xml.Attributes.ref_cam,'%d');

    switch xml.Attributes.ref_base
      case 'longest'
        s=seteoest(s,'depend',refCam);
      case {'x','y','z'}
        s=seteoest(s,'depend',refCam,xml.Attributes.ref_base);
      otherwise
        error(['DBAT XML operation set_datum error: Unknown reference ' ...
               'base ''%s'''],xml.Attributes.ref_base);
    end
end
