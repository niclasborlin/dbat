function s=parsesetinitialiovalues(s,xml)
%PARSESETINITIALIOVALUES Parse and set initial IO values for a bundle operation.
%
%   S=PARSESETINITIALIOVALUES(S,XML) parses the IO subblock of the
%   SET_INITIAL_VALUES block of a DBAT script operation. The DBAT
%   structure S is updated with the new initial values.
%
%   The IO subblock can contain the following fields:
%       cc     - camera constant in camera units, or the string
%                'focal' or 'default'. Either of the strings sets
%                the initial value for the camera constant to the
%                focal length.
%       pp     - comma-separated (x,y) with principal point in
%                camera units. The special string 'default' will 
%                set pp as the center of the sensor. 
%       K      - nK-vector with radial lens distortion
%                coefficients. The string 'default' corresponds to
%                all zeros.
%       P      - nP-vector with radial lens distortion
%                coefficients. The string 'default' corresponds to
%                all zeros.
%       aspect - scalar with pixel aspect ratio. The string 'default'
%                sets the aspect ratio to unity.
%       skew   - scalar with skew. The string 'default' sets a zero value.
%       all    - If equal to the string 'loaded', uses the
%                pre-loaded values. If equal to 'default', the
%                default values indicated above are set.
%   
%   An abbrivated IO subblock with the string 'loaded' or 'default' is
%   equivalent to 'all'/'loaded' and 'all'/'default', respectively.

narginchk(2,2);

knownFields={'Text','all','cc','pp','K','P','aspect','skew','all'};
[ok,msg]=checkxmlfields(xml,knownFields,false(size(knownFields)));
if ~ok, error('DBAT XML script set_initial_values/IO error: %s',msg); end

% Check for abbreviated block
if isfield(xml,'Text')
    switch xml.Text
      case {'loaded','default'}
        % Translate to <all>loaded</all>
        xml=struct('all',struct('Text',xml.Text));
      otherwise
        error('DBAT XML set initial values/IO error: Unknown string ''%s''',...
              xml.Text);
    end
end

% Parse each subblock
fn=fieldnames(xml);
for i=1:length(fn)
    switch fn{i}
      case 'all'
        % Copy all loaded prior IO values.
        s.IO.val=s.prior.IO.val;
      otherwise
        error('DBAT XML set initial values/IO error: Unknown field ''%s''',...
              fn{i});
    end        
end

