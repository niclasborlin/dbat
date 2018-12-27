function s=setcamvals(s,varargin)
%SETCAMVALS Set initial camera values in DBAT struct.
%
%   S=SETCAMVALS(S,'prior') sets the initial IO values for all cameras
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

% Check for 'prior' argument.
if length(varargin)>=i && ischar(varargin{i}) && strcmp(varargin{i},'prior')
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
        error('Bad parameter %d: CC must be numeric',i)
    end
end

% Parse individual parameter assignments
while length(varargin)>=i
    if length(varargin)<i+1
        error('Trailing arguments must come in pairs');
    end
    if ~ischar(varargin{i})
        error('Parameter name (arg %d) must be a char array',i);
    end
    % Check for simple arg.
    ii=find(strcmp(varargin{i},{'cc','px','py','as','sk'}));
    % Not found, check Ki, Pi
    if isempty(ii) && ~isempty(varargin{i})
        switch varargin{i}(1)
          case 'K'
            n=str2double(varargin{i}(2:end));
            if n<1 || n>s.IO.model.nK
                error('K number out of range');
            end
            ii=5+n;
          case 'P'
            n=str2double(varargin{i}(2:end));
            if n<1 || n>s.IO.model.nP
                error('P number out of range');
            end
            ii=5+s.IO.model.nK+n;
        end
    end
    if isempty(ii)
        error('Bad parameter %d: %s',i,varargin{i});
    end
    i=i+1;
    % Verify next parameter is numeric.
    if ~isnumeric(varargin{i})
        error('Parameter value (arg %d) must be numeric',i);
    end
    s.IO.val(ii,ix)=varargin{i};
    i=i+1;
end
