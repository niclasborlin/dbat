function [ok,msg]=checkxmlfields(s,fields,mandatory)
%CHECKXMLFIELDS Check fields loaded from an XML file.
%
%   OK=CHECKXMLFIELDS(S,FIELDS) checks that the structure S contains
%   only fields listed in the N-cell string array FIELDS. To test for
%   a single field, FIELDS may be a char array. All FIELDS are assumed
%   to be mandatory. Use CHECKXMLFIELDS(S,FIELDS,MAN), where MAN is a
%   logical N-vector. The function returns OK as TRUE if all mandatory
%   fields are present and no unknown field.
%
%   [OK,MSG]=... also returns an error message, indicating the
%   offending field name(s).

if ischar(fields), fields={fields}; end
if nargin<3, mandatory=true(size(fields)); end

msg='';

% Check that all mandatory fields are present.
if any(mandatory)
    mandatoryFields=fields(mandatory);
    tf=isfield(s,mandatoryFields);

    ok=all(tf);
    if ~ok
        msg=sprintf('Missing field ''%s''',mandatoryFields{find(~tf,1)});
        return;
    end
end

% Check that no unknown field is present.
df=setdiff(fieldnames(s),fields);

ok=isempty(df);
if ~ok
    msg=sprintf('Unsupported field ''%s''',df{1});
end
