function s=setcammodel(s,model,ix)
%SETCAMMODEL Set the camera model in a DBAT structure.
%
%   S=SETCAMMODEL(S,MODEL) sets the camera model to MODEL for all
%   cameras in S. Model can be a scalar or vector of appropriate
%   length. Use S=SETCAMMODEL(S,MODEL,IX) to set the camera model for
%   the cameras specified in the index vector IX only. A value of
%   IX='all' corresponds to all cameras.
%
%See also: SETCAMVALS.

if nargin<3, ix='all'; end

if ischar(ix)
    ix=1:size(s.IO.val,2);
end

s.IO.model.distModel(ix)=model;
