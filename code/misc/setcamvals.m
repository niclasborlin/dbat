function s=setcamvals(s,varargin)
%SETCAMVALS Set initial camera values in DBAT struct.
%
%   S=SETCAMVALS(S,'loaded') sets the initial IO values for all cameras
%   to their loaded prior values.
%
%   S=SETCAMVALS(S,'default',CC) sets the initial IO values to their
%   default values with camera constant CC. The other default values
%   are zero, except the principal point which is at the center of
%   the sensor.
%
%   S=SETCAMVALS(...,<param1>,<val1>,...) specifies invidual initial
%   values for each parameter. For acceptable parameter strings, see
%   BUILDPARAMTYPES.
%
%   S=SETCAMVALS(S,IX,...) sets the camera parameters for the
%   cameras specified in the index vector IX only. IX='all'
%   corresponds to the default: all cameras.
%
%See also: BUILDPARAMTYPES.

ix=1:size(s.IO.val,2);

% Next argument.
i=1;

if length(varargin)>=i
    if isnumeric(varargin{i}) && isvector(varargin{i})
        ix=varargin{i};
        i=i+1;
    elseif ischar(varargin{i}) && strcmp(varargin{i},'all')
        i=i+1;
        ix=1:size(s.IO,val,2);
    end
end

% Check for 'loaded' argument.
if length(varargin)>=i && ischar(varargin{i}) && strcmp(varargin{i},'loaded')
    i=i+1;
    s.IO.val(:,ix)=s.prior.IO.val(:,ix);
end

% Check for 'default' argument.
if length(varargin)>=i && ischar(varargin{i}) && strcmp(varargin{i},'default')
    i=i+1;
    if length(varargin)>=i && isnumeric(varargin{i})
        s.IO.val(1,ix)=varargin{i};
        s.IO.val(2:3,ix)=0.5*diag([1,-1])*s.IO.sensor.ssSize(:,ix);
        s.IO.val(4:end,ix)=0;
        i=i+1;
    else
        error('SETCAMVALS: Bad parameter %d: CC must be numeric',i)
    end
end

% Parse individual parameter assignments
while length(varargin)>=i
    if length(varargin)<i+1
        error('SETCAMVALS: Trailing arguments must come in pairs');
    end
    param=varargin{i};
    arg=varargin{i+1};
    % Verify parameter is char and argument is numeric.
    if ~ischar(param)
        error('SETCAMVALS: Parameter name (arg %d) must be a char array',i);
    end
    if ~isnumeric(arg)
        error('SETCAMVALS: Parameter value (arg %d) must be numeric',i+1);
    end
    % Check for simple arg.
    ii=find(strcmp(param,{'cc','px','py','as','sk'}));
    % Aspect, skew not defined for some models: Argument must be zero.
    mustBeZero=any(ismember(ii,[4,5])) && any(abs(s.IO.model.distModel(ix))<3);
    if isempty(ii) && ~isempty(param)
        % Not found, check Ki, Pi
        switch varargin{i}(1)
          case 'K'
            n=str2double(param(2:end));
            if n<1 || n>s.IO.model.nK
                error('SETCAMVALS: K number out of range');
            end
            ii=5+n;
          case 'P'
            n=str2double(param(2:end));
            if n<1 || n>s.IO.model.nP
                error('SETCAMVALS: P number out of range');
            end
            ii=5+s.IO.model.nK+n;
        end
    end
    if isempty(ii)
        error('SETCAMVALS: Bad parameter %d: ''%s''',i,param);
    end
    if mustBeZero && any(arg~=0)
        models=s.IO.model.distModel(ix);
        badModel=models(abs(models)<3);
        error('SETCAMVALS: Camera model %d only supports zero %s', ...
              badModel(1),param);
    end
    s.IO.val(ii,ix)=arg;
    i=i+2;
end
