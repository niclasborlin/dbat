function s=parsesetbundleestop(s,xml)
%PARSESETBUNDLEESTOP Parse what OP parameters should be estimated by the bundle.
%
%   S=PARSESETBUNDLEESTOP(S,XML) parses the OP subblock of the
%   SET_BUNDLE_ESTIMATE_PARAMS block of a DBAT script operation and
%   sets the parameters to estimate in the DBAT structure S.
%
%   The OP subblock can contain the fields 'all', 'x', 'y', and 'z'.
%   Each field can contain the strings 'true', 'false', or 'default'.
%   A value of 'true' sets all corresponding parameters to be
%   estimated. A value of 'false' sets all corresponding parameters
%   not to be estimated. The value 'default' sets the corresponding
%   parameters to be estimated if any prior observation is non-fixed.
%   
%   An abbrivated OP subblock with the string 'true'/'false'/'default'
%   affects all parameters.

narginchk(2,2);

knownFields={'Text','all','x','y','z'};
[ok,msg]=checkxmlfields(xml,knownFields,false(size(knownFields)));
if ~ok, error('DBAT XML script set_initial_values/OP error: %s',msg); end

% Check for abbreviated block
if isfield(xml,'Text')
    switch xml.Text
      case {'true','false','default'}
        % Translate to <all>loaded</all>
        xml=struct('all',struct('Text',xml.Text));
      otherwise
        error(['DBAT XML set bundle estimate params/OP error: Unknown ' ...
               'string ''%s'''], xml.Text);
    end
end

% Parse each subblock
fn=fieldnames(xml);
for i=1:length(fn)
    field=xml.(fn{i});
    switch fn{i}
      case 'all'
        ix=1:3;
      case {'x','y','z'}
        ix=abs(fn{i})-abs('x')+1;
    end
    [ok,msg]=checkxmlfields(field,'Text');
    if ~ok
        error(['DBAT XML script set_bundle_estimate_params/OP/%s ' ...
               'error: %s'],fn{i},msg);
    end
    switch field.Text
      case {'true','false'}
        % Estimate all or nothing.
        est=strcmp(field.Text,'true');
        s.bundle.est.OP(ix,:)=est;
      case 'default'
        % Estimate all non-control points and all non-fixed control points.
        s.bundle.est.OP(ix,:)=repmat(~s.prior.OP.isCtrl,length(ix),1) | ...
            s.prior.OP.std(ix,:)~=0;
      otherwise
        error(['DBAT XML script set_bundle_estimate_params/OP/%s ' ...
               'error: Unknown string ''%s'''],fn{i},field.Text); 
    end
end

