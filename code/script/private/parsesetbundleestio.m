function s=parsesetbundleestio(s,xml)
%PARSESETBUNDLEESTIO Parse what IO parameters should be estimated by the bundle.
%
%   S=PARSESETBUNDLEESTIO(S,XML) parses the IO subblock of the
%   SET_BUNDLE_ESTIMATE_PARAMS block of a DBAT script operation and
%   sets the parameters to estimate in the DBAT structure S.
%
%   The IO subblock can contain the fields 'all', 'cc', 'pp',
%   'aspect', 'skew', 'K', 'P' plus 'K1', 'K2', etc. and 'P1', 'P2',
%   etc. Each field can contain the strings 'true' or 'false'. A value
%   of 'true' sets all corresponding parameters to be estimated. A
%   value of 'false' sets all corresponding parameters not to be
%   estimated.
%   
%   An abbrivated IO subblock with the string 'true'/'false' affects
%   all parameters.

narginchk(2,2);

knownFields=cat(2,{'Text','all','cc','pp','px','py','K','P', ...
                   'aspect','skew','all'},...
                arrayfun(@(x)sprintf('K%d',x),1:s.IO.model.nK, ...
                         'uniformoutput',false),...
                arrayfun(@(x)sprintf('P%d',x),1:s.IO.model.nP, ...
                         'uniformoutput',false));

[ok,msg]=checkxmlfields(xml,knownFields,false(size(knownFields)));
if ~ok, error('DBAT XML script set_bundle_estimate_params/IO error: %s',msg); end

% Check for abbreviated block
if isfield(xml,'Text')
    switch strip(xml.Text)
      case {'true','false'}
        % Translate to <all>true/false</all>
        xml=struct('all',struct('Text',xml.Text));
      otherwise
        error(['DBAT XML set_bundle_estimate_params/IO error: Unknown ' ...
               'string ''%s'''],xml.Text);
    end
end

% Parse each subblock
fn=fieldnames(xml);
for i=1:length(fn)
    % Every field should have a Text subfield and nothing else
    field=xml.(fn{i});
    [ok,msg]=checkxmlfields(field,'Text');
    if ~ok
        error('DBAT XML script set_bundle_estimate_params/IO error: %s',msg);
    end
    
    switch strip(field.Text)
      case 'true'
        val=true;
      case 'false'
        val=false;
      otherwise
        error(['DBAT XML script set_bundle_estimate_params/IO error: ' ...
               'Bad string %s'],field.Text);
    end        

    switch fn{i}
      case 'aspect'
        fn{i}='as';
      case 'skew'
        fn{i}='sk';
    end
    
    if val
        s=setcamest(s,fn{i});
    else
        s=setcamest(s,'not',fn{i});
    end
end

