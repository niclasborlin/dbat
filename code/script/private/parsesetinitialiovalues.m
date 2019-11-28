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
%       aspect - scalar with pixel aspect ratio. The string 'default'
%                sets the aspect ratio to unity.
%       skew   - scalar with skew. The string 'default' sets a zero value.
%       K      - nK-vector with radial lens distortion
%                coefficients. The string 'default' corresponds to
%                all zeros.
%       P      - nP-vector with radial lens distortion
%                coefficients. The string 'default' corresponds to
%                all zeros.
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
    switch strip(xml.Text)
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
    % Every field should have a Text subfield and nothing else
    sub=xml.(fn{i});
    [ok,msg]=checkxmlfields(sub,'Text');
    if ~ok, error('DBAT XML script set_initial_values/IO error: %s',msg); end
    
    switch fn{i}
      case 'all'
        switch strip(sub.Text)
          case 'loaded'
            % Copy all loaded prior IO values.
            s=setcamvals(s,'loaded');
          case 'default'
            % Set all values as default.
            s=setcamvals(s,'default',s.prior.IO.cams{1}.FocalLength);
          otherwise
            error(['DBAT XML set initial values/IO error: Unknown ' ...
                   '''all'' string ''%s'''],sub.Text);
        end
      case 'cc'
        switch strip(sub.Text)
          case {'focal','default'}
            s=setcamvals(s,'cc',s.prior.IO.cams{1}.FocalLength);
          case 'loaded'
            % Copy prior cc values.
            s=setcamvals(s,'cc',getcamvals(s,'prior','cc'));
          otherwise
            cc=sscanf(sub.Text,'%f');
            if isscalar(cc)
                s=setcamvals(s,'cc',cc);
            else
                error(['DBAT XML set initial values/IO/cc error: Bad ' ...
                       'cc value: %s'],sub.Text);
            end
        end
      case 'pp'
        switch strip(sub.Text)
          case 'default'
            s=setcamvals(s,'pp',0.5*diag([1,-1])*s.IO.sensor.ssSize);
          case 'loaded'
            % Copy prior pp values.
            s=setcamvals(s,'pp',getcamvals(s,'prior','pp'));
          otherwise
            pp=sscanf(sub.Text,'%f,');
            if length(pp)==2
                s=setcamvals(s,'pp',repmat(pp,1,size(s.IO.val,2)));
            else
                error(['DBAT XML set initial values/IO/pp error: Bad ' ...
                       'pp values: %s'],sub.Text);
            end
        end
      case 'aspect'
        switch strip(sub.Text)
          case 'default'
            s=setcamvals(s,'as',0);
          case 'loaded'
            % Copy prior aspect value.
            s=setcamvals(s,'as',getcamvals(s,'prior','as'));
          otherwise
            aspect=sscanf(sub.Text,'%f');
            if isscalar(aspect)
                % Use DBATCamera to convert aspect ratios.
                c=DBATCamera;
                c.AspectRatio=aspect;
                s=setcamvals(s,'as',AspectDiff(c));
            else
                error(['DBAT XML set initial values/IO/aspect error: Bad ' ...
                       'aspect value: %s'],sub.Text);
            end
        end
      case 'skew'
        switch strip(sub.Text)
          case 'default'
            s=setcamvals(s,'sk',0);
          case 'loaded'
            % Copy prior skew value.
            s=setcamvals(s,'sk',getcamvals(s,'prior','sk'));
          otherwise
            skew=sscanf(sub.Text,'%f');
            if isscalar(skew)
                s=setcamvals(s,'sk',skew);
            else
                error(['DBAT XML set initial values/IO/skew error: Bad ' ...
                       'skew value: %s'],sub.Text);
            end
        end
      case 'K'
        switch strip(sub.Text)
          case 'default'
            s=setcamvals(s,'K',0);
          case 'loaded'
            % Copy prior K values.
            s=setcamvals(s,'K',getcamvals(s,'prior','K'));
          otherwise
            K=sscanf(sub.Text,'%f,');
            if length(K)==s.IO.model.nK
                s=setcamvals(s,'K',repmat(K,1,size(s.IO.val,2)));
            else
                error(['DBAT XML set initial values/IO/K error: Bad ' ...
                       'K values: %s'],sub.Text);
            end
        end
      case 'P'
        switch strip(sub.Text)
          case 'default'
            s=setcamvals(s,'P',0);
          case 'loaded'
            % Copy prior P values.
            s=setcamvals(s,'P',getcamvals(s,'prior','P'));
          otherwise
            P=sscanf(sub.Text,'%f,');
            if length(P)==s.IO.model.nP
                s=setcamvals(s,'P',repmat(P,1,size(s.IO.val,2)));
            else
                error(['DBAT XML set initial values/IO/P error: Bad ' ...
                       'P values: %s'],sub.Text);
            end
        end
      otherwise
        error('DBAT XML set initial values/IO error: Unknown field ''%s''',...
              fn{i});
    end        
end

