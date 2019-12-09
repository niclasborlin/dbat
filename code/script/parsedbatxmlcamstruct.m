function cams=parsedbatxmlcamstruct(s)
%PARSEDBATXMLCAMSTRUCT Convert DBAT XML camera structure to DBAT camera structure
%
%   CAMS=PARSEDBATXMLCAMSTRUCT(S) parses the DBAT camera XML structure
%   S and returns a 1-cell array CAMS with a DBATCamera object. The
%   structure S should have fields with field names listed below with
%   a Text field holding a string.
%
%   TODO: CAMS=PARSEDBATXMLCAMSTRUCT(S), where S is an N-cell array of
%   DBAT camera XML structures, converts multiple cameras and returns
%   an N-cell array with DBATCameras.
%
%   The XML field names are:
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
%       model      - scalar with projection model number.
%
%   Additionally, the XML field 'all' may be used with the string
%   'default' to set all parameters are set to their default values,
%   i.e., cc to the focal length, pp to the center of the sensor, zero
%   aspectDiff (unit aspect), zero skew, and zero lens distortion.
%
%See also: DBATCamera.

narginchk(1,1);

if ~iscell(s)
    s={s};
end

% Allow 'all' as an extra XML field name.
XMLfieldNames={'id', 'name', 'unit', 'sensor', 'image', 'aspect', ...
               'nK', 'nP', 'focal', 'model', 'cc', 'pp', 'K', 'P', ...
               'skew', 'all','calibrated'};

% Pre-allocate return array.
cams=cell(1,length(s));

% For each XML camera...
for i=1:length(s)
    % Interpret the camera element.
    e=s{i};
    % Check that only expected field names are present.
    [ok,msg]=checkxmlfields(e,XMLfieldNames,false(size(XMLfieldNames)));
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

    % Initialize the camera.
    cam=DBATCamera;

    % Parse each field.
    for j=1:length(fn)
        field=fn{j};
        switch field
          case 'id'
            cam.Id=sscanf(e.id.Text,'%d');
          case 'name'
            cam.Name=e.name.Text;
          case 'unit'
            cam.Unit=e.unit.Text;
          case 'sensor'
            ss=strip(split(e.sensor.Text,','));
            if length(ss)~=2
                error(['DBAT camera XML error: Wrong number of sensor ' ...
                       'values: %s'],e.sensor.Text);
            end
            sensorVal=nan(1,2);
            if strcmp(ss{1},'auto')
                % Do nothing, first value is already nan.
            else
                sensorVal(1)=sscanf(ss{1},'%f');
            end
            sensorVal(2)=sscanf(ss{2},'%f');
            cam.SensorSize=sensorVal;
          case 'image'
            imageVal=sscanf(e.image.Text,'%d,')';
            if length(imageVal)~=2
                error(['DBAT camera XML error: Wrong number of image size ' ...
                       'values: %s'],e.image.Text);
            end
            cam.ImageSize=imageVal;
          case 'aspect'
            if strcmp(strip(e.aspect.Text),'auto')
                aspectVal=nan;
            else
                aspectVal=sscanf(e.aspect.Text,'%f');
            end
            cam.AspectRatio=aspectVal;
          case 'focal'
            cam.FocalLength=sscanf(e.focal.Text,'%f');
          case 'cc'
            if strcmp(e.cc.Text,'focal')
                cam.CameraConstant=cam.FocalLength;
            else
                cam.CameraConstant=sscanf(e.cc.Text,'%f');
            end
          case 'nK'
            nK=sscanf(e.nK.Text,'%d')';
            K=GetStorableK(cam);
            if length(K)>nK
                K=K(1:nK);
            elseif length(K)<nK
                K(end+1:nK)=nan;
            end
            cam=SetStorableK(cam,K);
          case 'K'
            K=sscanf(e.K.Text,'%f,')';
            cam=SetStorableK(cam,K);
          case 'nP'
            nP=sscanf(e.nP.Text,'%d')';
            P=GetStorableP(cam);
            if length(P)>nP
                P=P(1:nP);
            elseif length(P)<nP
                P(end+1:nP)=nan;
            end
            cam=SetStorableP(cam,P);
          case 'P'
            P=sscanf(e.P.Text,'%f,')';
            cam=SetStorableP(cam,P);
          case 'model'
            cam.Model=sscanf(e.model.Text,'%d');
          case 'skew'
            cam.Skew=sscanf(e.skew.Text,'%f');
          case 'pp'
            if strcmp(strip(e.pp.Text),'default')
                cam=SetStorablePrincipalPoint(cam,evalsensor(cam)/2);
            else
                cam=SetStorablePrincipalPoint(cam,sscanf(e.pp.Text,'%f,')');
            end
            if length(GetStorablePrincipalPoint(cam))~=2
                error(['DBAT camera XML error: Wrong number of principal ' ...
                       'point values: %s'],e.pp.Text);
            end
          case 'all'
            if strcmp(strip(e.all.Text),'default')
                % Use default for all parameters
                cam.CameraConstant=cam.FocalLength;
                cam=SetStorablePrincipalPoint(cam,evalsensor(cam)/2);
                cam.AspectRatio=1;
                cam.Skew=0;
                cam=SetStorableK(cam,zeros(1,nK(cam)));
                cam=SetStorableP(cam,zeros(1,nP(cam)));
            else
                error(['DBAT camera XML error: Bad string for ''all''' ...
                       'directive: %s'],e.all.Text);
            end
          case 'calibrated'
            cam.Calibrated=strcmp(strip(e.calibrated.Text),'yes');
        end
    end
    
    if isnan(cam.AspectRatio)
        cam.AspectRatio=evalaspect(cam);
    else
        cam.SensorSize=evalsensor(cam);
    end
    
    cams{i}=cam;
end

nKall=cellfun(@(x)x.nK,cams);
nPall=cellfun(@(x)x.nP,cams);

if min(nKall)~=max(nKall)
    % Upgrade all cameras with short K vectors.
    mK=max(nKall);
    for i=find(nKall<mK)'
        K=GetStorableK(cams{i});
        K(end+1:mK)=nan;
        cams{i}=SetStorableK(cams{i},K);
    end
end

if min(nPall)~=max(nPall)
    % Upgrade all cameras with short P vectors.
    mP=max(nPall);
    for i=find(nPall<mP)'
        P=GetStorableP(cams{i});
        P(end+1:mP)=nan;
        cams{i}=SetStorableP(cams{i},P);
    end
end


function aspect=evalaspect(cam)
%Return the aspect ratio, either as specified or computed from the
%sensor and image sizes.

if isnan(cam.AspectRatio)
    % Compute aspect
    pixelSize=PixelSize(cam);
    aspect=pixelSize(1)/pixelSize(2);
else
    aspect=cam.AspectRatio;
end

function sensor=evalsensor(cam)
%Return the sensor size, either as specified or computed from the
%sensor width, aspect ratio, and image size.

sensor=cam.SensorSize;

if isnan(sensor(1))
    % Compute sensor width
    imageSize=cam.ImageSize;
    sensor(1)=cam.AspectRatio*sensor(2)*imageSize(1)/imageSize(2);
end
