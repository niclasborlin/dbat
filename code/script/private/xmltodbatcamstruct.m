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
%       id         - integer camera id.
%       name       - string with camera name.
%       unit       - string with camera unit.
%       sensor     - comma-separated (width, height) in camera
%                    units. The special string 'auto' can be used
%                    for the width, to have the width computed from
%                    the image size and aspect parameter.
%       image      - comma-separated (width, height) in pixel units.
%       aspectDiff - scalar with 1-pixel aspect ratio. Defaults to
%                    zero, i.e., unit aspect ratio.
%       focal      - nominal focal length in camera units.
%       cc         - camera constant in camera units, or the string 'focal'.
%       pp         - comma-separated (x,y) with principal point in
%                    camera units. The special string 'default'
%                    will return pp as the center of the sensor.
%       nK         - number of radial coefficients for lens
%                    distortion. Defaults to zero. 
%       nP         - number of tangential coefficients for lens
%                    distortion. Defaults to zero. 
%       K          - nK-vector with radial lens distortion
%                    coefficients. Defaults to the empty vector. 
%       P          - nP-vector with tangential lens distortion
%                    coefficients. Defaults to the empty vector. 
%       skew       - scalar with skew. Defaults to 0.
%       model      - scalar with camera model number.
%
%   Additionally, the XML field 'all' may be used with the string
%   'default' to set all parameters are set to their default values,
%   i.e., cc to the focal length, pp to the center of the sensor, zero
%   aspectDiff (unit aspect), zero skew, and zero lens distortion.
%
%   Any missing field without a listed default above is returned as
%   blank strings or NaN vectors of appropriate sizes, depending on
%   type.
%
%   If multiple cameras are present, the nK and nP values will be
%   adjusted upwards to match if necessary.

narginchk(1,1),

if ~iscell(s)
    s={s};
end

% Create blank camera struct. Struct field names are the same as
% the XML field names.
fieldNames={'id','name','unit','sensor','image','aspectDiff','nK','nP',...
            'focal','model','cc','pp',    'K',       'P',       'skew'};
defaults  ={nan, '',    '',    nan(1,2),nan(1,2),nan,        nan, nan,...
            nan,    nan,    nan, nan(1,2),zeros(1,0),zeros(1,0),nan}; 
XMLfieldNames={'id','name','unit','sensor','image','aspect','nK','nP',...
               'focal','model','cc','pp',    'K',       'P',       'skew'};

blankCamera=cell2struct(defaults,fieldNames,2);

% Allow 'all' as an extra field name.

extraFields={'all'};
allFields=cat(2,XMLfieldNames,extraFields);

% Pre-allocate return array.
cams=repmat(blankCamera,1,length(s));

% For each XML camera...
for i=1:length(s)
    % Interpret the camera element.
    e=s{i};
    % Check that only expected field names are present.
    [ok,msg]=checkxmlfields(e,allFields,false(size(allFields)));
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

    for j=1:length(fn)
        field=fn{j};
        switch field
            case 'id'
              cam.id=sscanf(e.id.Text,'%d');
          case 'name'
            cam.name=e.name.Text;
          case 'unit'
            cam.unit=e.unit.Text;
          case 'sensor'
            ss=strip(split(e.sensor.Text,','));
            if length(ss)~=2
                error(['DBAT camera XML error: Wrong number of sensor ' ...
                       'values: %s'],e.sensor.Text);
            end
            if strcmp(ss{1},'auto')
                cam.sensor(1)=nan;
                sensorAuto=true;
            else
                cam.sensor(1)=sscanf(ss{1},'%f');
            end
            cam.sensor(2)=sscanf(ss{2},'%f');
          case 'image'
            cam.image=sscanf(e.image.Text,'%d,')';
            if length(cam.image)~=2
                error(['DBAT camera XML error: Wrong number of image size ' ...
                       'values: %s'],e.image.Text);
            end
          case 'aspect'
            if strcmp(strip(e.aspect.Text),'auto')
                cam.aspectDiff=nan;
            else
                aspect=sscanf(e.aspect.Text,'%f');
                cam.aspectDiff=1-aspect;
            end
          case 'focal'
            cam.focal=sscanf(e.focal.Text,'%f');
          case 'cc'
            if strcmp(e.cc.Text,'focal')
                cam.cc=cam.focal;
            else
                cam.cc=sscanf(e.cc.Text,'%f');
            end
          case 'nK'
            cam.nK=sscanf(e.nK.Text,'%d')';
            if isempty(cam.K)
                cam.K=nan(1,cam.nK);
            end
          case 'K'
            cam.K=sscanf(e.K.Text,'%f,')';
            if isnan(cam.nK)
                cam.nK=length(cam.K);
            elseif cam.nK~=length(cam.K)
                error(['DBAT camera XML error: Wrong number of K ' ...
                       'values: %s'],e.image.Text); 
            end
          case 'nP'
            cam.nP=sscanf(e.nP.Text,'%d')';
            if isempty(cam.P)
                cam.P=nan(1,cam.nP);
            end
          case 'P'
            cam.P=sscanf(e.P.Text,'%f,')';
            if isnan(cam.nP)
                cam.nP=length(cam.P);
            elseif cam.nP~=length(cam.P)
                error(['DBAT camera XML error: Wrong number of P ' ...
                       'values: %s'],e.image.Text); 
            end
          case 'model'
            cam.model=sscanf(e.model.Text,'%d');
          case 'skew'
            cam.skew=sscanf(e.skew.Text,'%f');
          case 'pp'
            if strcmp(strip(e.pp.Text),'default')
                cam.pp=evalsensor(cam)/2;
            else
                cam.pp=sscanf(e.pp.Text,'%f,')';
            end
            if length(cam.pp)~=2
                error(['DBAT camera XML error: Wrong number of principal ' ...
                       'point values: %s'],e.pp.Text);
            end
          case 'all'
            if strcmp(strip(e.all.Text),'default')
                % Use default for all parameters
                cam.cc=cam.focal;
                cam.pp=evalsensor(cam)/2;
                cam.aspectDiff=0;
                cam.skew=0;
                cam.K=zeros(1,cam.nK);
                cam.P=zeros(1,cam.nP);
            else
                error(['DBAT camera XML error: Bad string for ''all''' ...
                       'directive: %s'],e.all.Text);
            end
        end
    end
    
    if isnan(cam.aspectDiff)
        cam.aspectDiff=1-evalaspect(cam);
    elseif isnan(cam.sensor(1))
        cam.sensor=evalsensor(cam);
    else
        error(['DBAT camera XML error: Either aspect or sensor ' ...
               'height must be auto']);
    end
    
    cams(i)=cam;
end

nK=cat(1,cams.nK);
nP=cat(1,cams.nP);

if min(nK)~=max(nK)
    % Upgrade all cameras with short K vectors.
    mK=max(nK);
    for i=find(nK<mK)'
        cams(i).nK=mK;
        cams(i).K(mK)=0;
    end
end

if min(nP)~=max(nP)
    % Upgrade all cameras with short P vectors.
    mP=max(nP);
    for i=find(nP<mP)'
        cams(i).nP=mP;
        cams(i).P(mK)=0;
    end
end


function aspect=evalaspect(cam)
%Return the aspect ratio, either as specified or computed from the
%sensor and image sizes.

if isnan(cam.aspectDiff)
    % Compute aspect
    pixelSize=cam.sensor./cam.image;
    aspect=pixelSize(1)/pixelSize(2);
else
    aspect=1-cam.aspectDiff;
end

function sensor=evalsensor(cam)
%Return the sensor size, either as specified or computed from the
%sensor width, aspect ratio, and image size.

if isnan(cam.sensor(1))
    % Compute sensor width
    cam.sensor(1)=(1-cam.aspectDiff)*cam.sensor(2)*cam.image(1)/cam.image(2);
end

sensor=cam.sensor;

