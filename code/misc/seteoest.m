function s=seteoest(s,varargin)
%SETEOEST Set external orientation parameters to be estimated by the bundle.
%
%   S=SETEOEST(S,'all') modifies the S.bundle.est.EO field to indicate that all
%   EO parameters should be estimated by the bundle. Similarly,
%   S=SETEOEST(S,'none') sets all EO parameters to be fixed.
%
%   Individual parameters can be set/unset by S=SETEOEST(S,IX1,<param1>,...),
%   and S=SETEOEST(S,'not',IX1,<param1>,...), respectively. IX1 is a vector
%   with image indices. param1 is an EO parameter name:
%   - 'x', 'y', 'z'    - camera coordinate,
%   - 'pos'            - all camera coordinates,
%   - 'om', 'ph', 'ka' - Euler x-y-z angles,
%   - 'ang'            - all angles.
%   - 'all'            - all parameters.
%   - 'none'           - no parameters.
%   IX1 may be unspecified, in which case it is inherited from the previous
%   argument. Use IX1=[] to indicate all images.
%
%   S=SETEOEST(S,'depend',IMNO) sets up for a dependent relative orientation,
%   i.e. all parameters of the base camera IMNO and the coordinate with largest
%   offset to the base camera is fixed. Use S=SETEOEST(S,'depend',IMNO,<param>)
%   to specify which offset coordinate 'x', 'y', or 'z' to use. The base camera
%   IMNO defaults to 1 if unspecified.

if ~isempty(varargin) && ischar(varargin{1}) && strcmp(varargin{1},'depend')
    % Handle 'depend' argument separately.
    s=setdepend(s,varargin{2:end});
    return;
end

% Estimate parameters until 'not' found.
doEst=true;

% Default to all cameras.
camIx=1:size(s.IO.val,2);

% Parse each argument
while ~isempty(varargin)
    if isnumeric(varargin{1})
        camIx=varargin{1};
        if isempty(camIx)
            camIx=1:size(s.IO.val,2);
        end
    elseif ischar(varargin{1})
        ii=[];
        val=doEst;
        switch varargin{1}
          case 'not'
            if doEst==false
                error('SETEOEST: Cannot specify not twice');
            end
            doEst=false;
          case 'x'
            ii=1;
          case 'y'
            ii=2;
          case 'z'
            ii=3;
          case 'pos'
            ii=1:3;
          case 'om'
            ii=4;
          case 'ph'
            ii=5;
          case 'ka'
            ii=6;
          case 'ang'
            ii=4:6;
          case 'all'
            ii=1:6;
          case 'none'
            ii=1:6;
            val=false;
          otherwise
            error('SETEOEST: Bad argument string %s',varargin{1});
        end
        if ~isempty(ii)
            % Set/clear if argument was specified.
            s.bundle.est.EO(ii,camIx)=val;
        end
    else
        error('SETEOEST: Bad argument type');
    end
    varargin(1)=[];
end


% Handle 'depend' argument.
function s=setdepend(s,varargin)

% Get base camera.
camNo=1;
if ~isempty(varargin) && isnumeric(varargin{1})
    camNo=varargin{1};
    if ~isscalar(camNo)
        error('Base camera index IMNO must be scalar.');
    end
    varargin(1)=[];
end

% Get offset coordinate.
cIx=1:3;
if ~isempty(varargin)
    switch varargin{1}
      case 'x'
        cIx=1;
      case 'y'
        cIx=2;
      case 'z'
        cIx=3;
      case 'pos'
        cIx=1:3;
      otherwise
        error('SETEOEST: Bad position argument for depend.')
    end
end

% Find largest offset coordinate among the allowed coordinates.
basePos=s.EO.val(1:3,camNo);
offset=s.EO.val(cIx,:)-repmat(basePos(cIx),1,size(s.EO.val,2));

[i,j]=find(offset==max(offset(:)));

% 'depend' implies all except base camera and largest offset.
s.bundle.est.EO(:)=true;
s.bundle.est.EO(:,camNo)=false;
s.bundle.est.EO(cIx(i),j)=false;
