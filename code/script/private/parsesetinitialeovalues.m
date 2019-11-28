function s=parsesetinitialeovalues(s,xml)
%PARSESETINITIALEOVALUES Parse and set initial EO values for a bundle operation.
%
%   S=PARSESETINITIALEOVALUES(S,XML) parses the EO subblock of the
%   SET_INITIAL_VALUES block of a DBAT script operation. The DBAT
%   structure S is updated with the new initial values.
%
%   The EO subblock can contain the field 'all' that can only contain
%   the string 'loaded'. This sets all loaded EO parameters as initial
%   values.
%   
%   An abbrivated EO subblock with the string 'loaded' is equivalent
%   to 'all'/'loaded'.

narginchk(2,2);

knownFields={'Text','all'};
[ok,msg]=checkxmlfields(xml,knownFields,false(size(knownFields)));
if ~ok, error('DBAT XML script set_initial_values/EO error: %s',msg); end

% Check for abbreviated block
if isfield(xml,'Text')
    switch strip(xml.Text)
      case 'loaded'
        % Translate to <all>loaded</all>
        xml=struct('all',struct('Text','loaded'));
      otherwise
        error('DBAT XML set initial values/EO error: Unknown string ''%s''',...
              xml.Text);
    end
end

% Parse each subblock
fn=fieldnames(xml);
for i=1:length(fn)
    % Every field should have a Text subfield and nothing else
    sub=xml.(fn{i});
    [ok,msg]=checkxmlfields(sub,'Text');
    if ~ok, error('DBAT XML script set_initial_values/EO error: %s',msg); end
    
    switch fn{i}
      case 'all'
        switch strip(sub.Text)
          case 'loaded'
            % Copy all loaded prior EO values.
            s.EO.val=s.prior.EO.val;
          otherwise
            error(['DBAT XML set initial values/EO/all error: Unknown ' ...
                   'string ''%s'''],field.Text);
        end
      otherwise
        error('DBAT XML set initial values/EO error: Unknown field ''%s''',...
              fn{i});
    end        
end

