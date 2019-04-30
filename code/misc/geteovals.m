function v=geteovals(s,cams,varargin)
%GETEOVALS Get camera external orientation values by name from DBAT struct.
%
%   V=GETEOVALS(S,CAMS,PARAM) returns the M-by-N array V with the
%   camera external orientation parameters specified by the string
%   PARAM for the camera with indices in the N-vector CAMS. Use
%   CAMS='all' to get the values for all cameras. The number of rows M
%   are different with different PARAM values. Individual elements
%   (M=1) correspond to the following values PARAM strings:
%   - 'X', 'Y', 'Z'           - camera station coordinate.
%   - 'omega', 'phi', 'kappa' - camera station Euler angle.
%
%   For groups of values, the following PARAM strings may be used:
%   - 'pos'    - camera station coordinates [X;Y;Z], M=3.
%   - 'angles' - camera station euler angles [omega;phi;kappa], M=3.
%   
%   V=GETEOVALS(S,CAMS,PARAM1,PARAM2,...) can be used to return
%   arbitrary combination of values. Successive PARAM values will be
%   appended to the bottom of V.
%
%   The special PARAM string 'R' may be used to return a 3-by-3-by-N
%   array with rotation matrices corresponding to the world-to-camera
%   rotation for each camera. Similarly, the special PARAM string 'P'
%   may be used to return a 3-by-4-by-N array with the pose-only
%   camera matrices corresponding to the world-to-camera
%   transformation for each camera. 'R' or 'P' cannot be combined with
%   any other PARAM string.
%
%See also: BUILDPARAMTYPES, GETCAMVALS.

if ischar(cams)
    if strcmp(cams,'all')
        cams=1:size(s.EO.val,2);
    else
        error('Bad CAMS string.');
    end
end

% Check for 'R'
if any(strcmp('R',varargin))
    % 'R' must be by itself.
    if length(varargin)>1
        error('''R'' must appear as solitary PARAM string');
    end
    v=nan(3,3,length(cams));
    for i=1:length(cams)
        v(:,:,i)=rotmat(s.EO.val(4:6,cams(i)));
    end
elseif any(strcmp('P',varargin)) % Check for 'P'
    % 'P' must be by itself.
    if length(varargin)>1
        error('''P'' must appear as solitary PARAM string');
    end
    v=nan(3,4,length(cams));
    for i=1:length(cams)
        v(:,:,i)=rotmat(s.EO.val(4:6,cams(i)))*[eye(3),-s.EO.val(1:3,cams(i))];
    end
else
    ix=zeros(0,1);

    for i=1:length(varargin)
        j=[];
        switch varargin{i}
          case 'X'
            j=1;
          case 'Y'
            j=2;
          case 'Z'
            j=3;
          case 'omega'
            j=4;
          case 'phi'
            j=5;
          case 'kappa'
            j=6;
          case 'pos'
            j=(1:3)';
          case 'angles'
            j=(4:6)';
        end
        if isempty(j)
            error('Bad PARAM ''%s''',varargin{i});
        end
        ix=[ix;j];
    end
    v=s.EO.val(ix,cams);
end
