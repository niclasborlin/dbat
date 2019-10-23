function cams=xmltodbatcamstruct(s)
%XMLTODBATCAMSTRUCT Convert DBAT XML camera structure to DBAT camera structure
%
%   CAMS=XMLTODBATCAMSTRUCT(S) converts the N-cell array S with XML
%   camera information to a DBAT camera structure N-array CAMS. Each
%   element E in S should contain a structure in DBAT camera XML file
%   format, i.e., with field names listed below with a Text field
%   holding a string. The string is converted to a same-named field in
%   the corresponding element in CAMS.
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
%       K      - nK-vector with radial lens distortion
%                coefficients. Defaults to the empty vector.
%       P      - nP-vector with tangential lens distortion
%                coefficients. Defaults to the empty vector.
%       skew   - scalar with skew. Defaults to 0.
%       model  - scalar with camera model number.
%
%   Any missing field without a listed default above is returned as
%   blank strings or NaN vectors of appropriate sizes, depending on type.
%
%   The returned K and P fields are standardized to have the same
%   length for all elements in CAM.

% Create blank camera struct.
fieldNames={'id','name','unit','sensor','image','aspect','focal', ...
            'model','cc','pp',    'K',       'P',       'skew'};
defaults  ={nan, '',    '',    nan(1,2),nan(1,2),nan,    nan, ...
            nan,    nan, nan(1,2),zeros(1,0),zeros(1,0),nan}; 

blankCamera=cell2struct(defaults,fieldNames,2);

longestK=0;
longestP=0;

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
        cam.pp=sscanf(e.pp.Text,'%f,')';
        if length(cam.pp)~=2
            error(['DBAT camera XML error: Wrong number of principal point ' ...
                   'values: %s'],e.pp.Text);
        end
    end
    if isfield(e,'K')
        cam.K=sscanf(e.K.Text,'%f,')';
        if length(cam.K)>longestK
            longestK=length(cam.K);
        end
    end
    if isfield(e,'P')
        cam.P=sscanf(e.P.Text,'%f,')';
        if length(cam.P)>longestP
            longestP=length(cam.P);
        end
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
    
    cams(i)=cam;
end

% Cannot have only P1.
if longestP==1
    longestP=2;
end

% Normalize lengths of K and P vectors.
for i=1:length(cams)
    if length(cams(i).K)<longestK
        cams(i).K(longestK)=0;
    end
    if length(cams(i).P)<longestP
        cams(i).P(longestP)=0;
    end
end
