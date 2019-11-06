function v=getcamvals(s,param,ix)
%GETCAMVALS Get camera values from DBAT struct.
%
%   GETCAMVALS(S,PARAM) returns the camera value specified by the
%   string PARAM for all cameras in S. Use GETCAMVALS(S,PARAM,IX) to
%   return the values for the cameras specified by the index vector
%   IX. A value of IX='all' corresponds to all cameras. The allowed
%   parameter string values are specified in BUILDPARAMTYPES, in
%   addition to 'pp', 'lin', 'K', and 'P', that corresponds to the
%   principal point, the linear parameters, the K values, and the P
%   values, respectively. The linear parameters comprise the camera
%   constant, the principal point, skew, and aspect.
%
%See also: BUILDPARAMTYPES.

if nargin<3, ix='all'; end

if ischar(ix)
    ix=1:size(s.IO.val,2);
end

% Check for simple arg.
ii=find(strcmp(param,{'cc','px','py','as','sk'}));
if isempty(ii)
    % Look for grouped parameters.
    switch param
      case 'lin'
        ii=1:5;
      case 'pp'
        ii=2:3;
      case 'K'
        ii=5+(1:s.IO.model.nK);
      case 'P'
        ii=5+s.IO.model.nK+(1:s.IO.model.nP);
    end
end
if isempty(ii)
    % Not found, check Ki, Pi
    switch param(1)
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
    error('GETCAMVALS: Bad parameter %d: ''%s''',i,param);
end

v=s.IO.val(ii,ix);
