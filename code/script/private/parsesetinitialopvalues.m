function s=parsesetinitialopvalues(s,xml)
%PARSESETINITIALOPVALUES Parse and set initial OP values for a bundle operation.
%
%   S=PARSESETINITIALOPVALUES(S,XML) parses the OP subblock of the
%   SET_INITIAL_VALUES block of a DBAT script operation. The DBAT
%   structure S is updated with the new initial values.
%
%   The OP subblock can contain the field 'all', which can only
%   contain the string 'loaded'. The field 'all' sets all loaded OP
%   coordinates as initial values.
%   
%   An abbrivated OP subblock with the string 'loaded' is equivalent
%   to 'all'/'loaded'.

narginchk(2,2);

knownFields={'Text','all'};
[ok,msg]=checkxmlfields(xml,knownFields,false(size(knownFields)));
if ~ok, error('DBAT XML script set_initial_values/OP error: %s',msg); end

% Check for abbreviated block
if isfield(xml,'Text')
    switch xml.Text
      case 'loaded'
        % Translate to <all>loaded</all>
        xml=struct('all',struct('Text','loaded'));
      otherwise
        error('DBAT XML set initial values/OP error: Unknown string ''%s''',...
              xml.Text);
    end
end

% Parse each subblock
fn=fieldnames(xml);
for i=1:length(fn)
    field=xml.(fn{i});
    switch fn{i}
      case 'all'
        [ok,msg]=checkxmlfields(field,'Text');
        if ~ok
            error('DBAT XML script set_initial_values/OP/all error: %s',msg);
        end
        switch field.Text
          case 'loaded'
            % Copy all loaded prior OP values.
            s.OP.val=s.prior.OP.val;
          otherwise
            error(['DBAT XML set initial values/OP/all error: Unknown ' ...
                   'string ''%s'''],field.Text);
        end
      otherwise
        error('DBAT XML set initial values/OP error: Unknown field ''%s''',...
              fn{i});
    end        
end

