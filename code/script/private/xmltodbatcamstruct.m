function cams=xmltodbatcamstruct(s)
%XMLTODBATCAMSTRUCT Convert DBAT XML camera structure to DBAT camera structure
%
%   CAMS=XMLTODBATCAMSTRUCT(S) converts the DBAT camera XML structure
%   S to a DBAT camera structure CAMS. The structure S should have
%   fields with field names listed below with a Text field holding a
%   string. The string is converted to a same-named field in the
%   corresponding element in CAMS.
%
%   CAMS=XMLTODBATCAMSTRUCT(S), where S is an N-cell array of DBAT
%   DBAT camera XML structures, converts multiple cameras and
%   returns an N-array with camera data.
%
%   The field names are:
%       id     - integer camera id.
%       name   - string with camera name.
%       unit   - string with camera unit.
%       sensor - comma-separated (width, height) in camera units.
%                The special string 'auto' can be used for the
%                width, to have the width computed from the image
%                size and aspect parameter.
%       image  - comma-separated (width, height) in pixel units.
%       aspect - scalar with pixel aspect ratio, or 'auto' (default).
%       focal  - nominal focal length in camera units.
%       cc     - camera constant in camera units.
%       pp     - comma-separated (x,y) with principal point in
%                camera units. The special string 'default' will
%                return pp as the center of the sensor.
%       nK     - number of radial coefficients for lens distortion.
%                Defaults to zero.
%       nP     - number of tangential coefficients for lens distortion.
%                Defaults to zero.
%       K      - nK-vector with radial lens distortion
%                coefficients. Defaults to the empty vector.
%       P      - nP-vector with tangential lens distortion
%                coefficients. Defaults to the empty vector.
%       skew   - scalar with skew. Defaults to 0.
%       model  - scalar with camera model number.
%
%   Any missing field without a listed default above is returned as
%   blank strings or NaN vectors of appropriate sizes, depending on type.

narginchk(1,1),

if ~iscell(s)
    s={s};
end

% Create blank camera struct.
fieldNames={'id','name','unit','sensor','image','aspect','nK','nP',...
            'focal','model','cc','pp',    'K',       'P',       'skew'};
defaults  ={nan, '',    '',    nan(1,2),nan(1,2),nan,    nan, nan,...
            nan,    nan,    nan, nan(1,2),zeros(1,0),zeros(1,0),nan}; 

blankCamera=cell2struct(defaults,fieldNames,2);

% Pre-allocate return array.
cams=repmat(blankCamera,1,length(s));

for i=1:length(s)
    % Interpret each element.
    e=s{i};
    % Check that only expected field names are present.
    [ok,msg]=checkxmlfields(e,fieldNames,false(size(fieldNames)));
    if ~ok
        error('DBAT camera XML error: %s',msg);
    end
    % Check that each field has a 'Text' field and nothing else.
    fn=fieldnames(e);
    for j=1:length(fn)
        if ~checkxmlfields(e.(fn{j}),'Text')
            error(['DBAT camera XML error: Missing ''Text'' field ' ...
                   'in %s'],fn(j));
        end
    end
    % Parse each field.
    cam=blankCamera;
    
    if isfield(e,'id')
        cam.id=sscanf(e.id.Text,'%d');
    end
    if isfield(e,'name')
        cam.name=e.name.Text;
    end
    if isfield(e,'unit')
        cam.unit=e.unit.Text;
    end
    if isfield(e,'sensor')
        ss=strip(split(e.sensor.Text,','));
        if length(ss)~=2
            error(['DBAT camera XML error: Wrong number of sensor ' ...
                   'values: %s'],e.sensor.Text);
        end
        if strcmp(ss{1},'auto')
            cam.sensor(1)=nan;
        else
            cam.sensor(1)=sscanf(ss{1},'%f');
        end
        cam.sensor(2)=sscanf(ss{2},'%f');
    end
    if isfield(e,'image')
        cam.image=sscanf(e.image.Text,'%d,')';
        if length(cam.image)~=2
            error(['DBAT camera XML error: Wrong number of image size ' ...
                   'values: %s'],e.image.Text);
        end
    end
    if isfield(e,'aspect')
        if strcmp(strip(e.aspect.Text),'auto')
            cam.aspect=nan;
        else
            cam.aspect=sscanf(e.aspect.Text,'%f');
        end
    end
    if isfield(e,'focal')
        cam.focal=sscanf(e.focal.Text,'%f');
    end
    if isfield(e,'cc')
        cam.cc=sscanf(e.cc.Text,'%f');
    end
    if isfield(e,'pp')
        if strcmp(strip(e.pp.Text),'default')
            cam.pp=cam.sensor/2;
        else
            cam.pp=sscanf(e.pp.Text,'%f,')';
        end
        if length(cam.pp)~=2
            error(['DBAT camera XML error: Wrong number of principal point ' ...
                   'values: %s'],e.pp.Text);
        end
    end
    if isfield(e,'nK')
        cam.nK=sscanf(e.nK.Text,'%d')';
    end
    if isfield(e,'nP')
        cam.nP=sscanf(e.nP.Text,'%d')';
    end
    if isfield(e,'K')
        cam.K=sscanf(e.K.Text,'%f,')';
    end
    if isfield(e,'P')
        cam.P=sscanf(e.P.Text,'%f,')';
    end
    if isfield(e,'model')
        cam.model=sscanf(e.model.Text,'%d');
    end
    if isfield(e,'skew')
        cam.skew=sscanf(e.skew.Text,'%f');
    end

    if isnan(cam.aspect)
        % Compute aspect
        pixelSize=cam.sensor./cam.image;
        cam.aspect=pixelSize(1)/pixelSize(2);
    end
    
    if isnan(cam.sensor(1))
        % Compute sensor width
        cam.sensor(1)=cam.aspect*cam.sensor(2)*cam.image(1)/cam.image(2);
    end

    % Clean up K and P fields to match nK and nP.
    if ~(cam.nK==length(cam.K))
        if isnan(cam.nK)
            % K specified but not nK => set nK to length of K.
            cam.nK=length(cam.K);
        elseif isempty(cam.K)
            % nK specified but not K => set K to zero vector of length nK.
            cam.K=zeros(1,cam.nK);
        end
        if cam.nK~=length(cam.K)
            error('DBAT camera XML error: K length mismatch');
        end
    end
    if ~(cam.nP==length(cam.P))
        if isnan(cam.nP)
            % P specified but not nP => set nP to length of P.
            cam.nP=length(cam.P);
        elseif isempty(cam.P)
            % nP specified but not P => set P to zero vector of length nP.
            cam.P=zeros(1,cam.nP);
        end
        if cam.nP~=length(cam.P)
            error('DBAT camera XML error: P length mismatch');
        end
        if cam.nP==1
            error('DBAT camera XML error: P cannot be one');
        end
    end

    cams(i)=cam;
end
