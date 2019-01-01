function s=setcamest(s,varargin)
%SETCAMEST Set camera parameters to be estimated by the bundle.
%
%   S=SETCAMEST(S,'all') modifies the S.bundle.est.IO field to indicate that
%   all camera parameters should be estimated by the bundle. Use
%   S=SETCAMEST(S,'all','not',<param1>,<param2>,...) to specify that all
%   parameters except the specified should be estimated. See BUILDPARAMTYPES
%   for valid parameter names. Additionally, 'K' and 'P' means all K and P
%   parameters, respectively. The string 'af' means all affine parameters
%   ('cc', 'px', 'py', and optionally 'as', 'sk'), supported by the current
%   distortion model for the camera. The string 'pp' means 'px' and 'py'.
%
%   S=SETCAMEST(S,'none') modifies the S.bundle.est.IO field to indicate that
%   no camera parameters should be estimated in the bundle. Use
%   S=SETCAMEST(S,'none',<param1>,<param2>,...) to specify that only the
%   listed parameters should be estimated.
%
%   S=SETCAMEST(S,<param1>,...) or S=SETCAMEST(S,'not',<param1>,...) leaves
%   the estimation status for the unlisted parameters unchanged.
%
%   When parameters are specified individually, an error is thrown if an
%   parameter unsupported by the distortion model is to be estimated.
%
%   For the K lens distortion parameters, the first parameters of the sequence
%   K1, ..., KN must be estimated together. Thus, specifying that 'K2'
%   should be estimated implies that 'K1' will also be estimated.
%   Specifying that 'K2' should not be estimated implies that K3, ..., will
%   not be estimated.
%
%   The same restriction applies to the P sequence, except that P1 and P2
%   are always estimated together.
%
%   S=SETCAMEST(S,IX,...) restricts the estimation selection to the cameras
%   specified in the index vector IX only.
%
%See also: BUILDPARAMTYPES.

% Default to all cameras.
camIx=1:size(s.IO.val,2);

% Camera index vector.
if ~isempty(varargin) && isnumeric(varargin{1}) && isvector(varargin{1})
    camIx=varargin{1};
    varargin(1)=[];
end

% Initially, set to estimate parameters.
doEst=true;

% Does the model support aspect and skew?
models=s.IO.model.distModel(camIx);
supportsSkewAspect=abs(models)>=3;

for i=1:length(varargin)
    ii=[];
    val=nan;
    switch varargin{i}
      case 'all'
        ii=1:size(s.bundle.est.IO,1);
        % Mask out skew and aspect for models that do not support it.
        val=repmat(true,length(ii),nnz(camIx));
        val(4:5,:)=val(4:5,:) & repmat(supportsSkewAspect,2,1);
      case 'none'
        ii=1:size(s.bundle.est.IO,1);
        val=false;
      case 'af'
        ii=1:5;
        % Mask out skew and aspect for models that do not support it.
        val=repmat(doEst,length(ii),nnz(camIx));
        val(4:5,:)=val(4:5,:) & repmat(supportsSkewAspect,2,1);
      case 'not'
        if doEst==false
            error('SETCAMEST: Cannot specify not twice');
        end
        doEst=false;
      case 'K'
        ii=5+(1:s.IO.model.nK);
        val=doEst;
      case 'P'
        ii=5+s.IO.model.nK+(1:s.IO.model.nP);
        val=doEst;
      case 'pp'
        ii=2:3;
        val=doEst;
      otherwise
        ii=find(strcmp(varargin{i},{'cc','px','py','as','sk'}));
        if doEst && any(ismember(ii,[4,5])) && any(~supportsSkewAspect)
            badModels=models(~supportsSkewAspect);
            error('SETCAMEST: Model %d does not support skew, aspect.',...
                  badModels(1));
        end
        % Not found, check Ki, Pi
        if isempty(ii) && ~isempty(varargin{i})
            switch varargin{i}(1)
              case 'K'
                n=str2double(varargin{i}(2:end));
                if n<1 || n>s.IO.model.nK
                    error('SETCAMEST: K number out of range');
                end
                if doEst
                    % Set K1..Kn
                    ii=5+(1:n);
                else
                    % Clear from Kn onwards
                    ii=5+(n:s.IO.model.nK);
                end
              case 'P'
                n=str2double(varargin{i}(2:end));
                if n<1 || n>s.IO.model.nP
                    error('SETCAMEST: P number out of range');
                end
                if doEst
                    % Include P2 if P1 was specified.
                    if n==1, n=2; end
                    % Set P1..Pn
                    ii=5+s.IO.model.nK+(1:n);
                else
                    % Include P1 if P2 was specified.
                    if n==2, n=1; end
                    % Clear from Pn onwards
                    ii=5+s.IO.model.nK+(n:s.IO.model.nP);
                end
            end
        end
        if isempty(ii)
            error('SETCAMEST: Bad parameter ''%s''',varargin{i});
        end
        val=doEst;
    end
    if ~isempty(ii)
        s.bundle.est.IO(ii,camIx)=val;
    end
end
