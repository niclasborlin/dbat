function [s,id,res]=forwintersect(s0,ids)
%FORWINTERSECT Perform forward intersection on points in a project.
%
%   S=FORWINTERSECT(S0,ID) computes the 3D coordinates of all points in
%   the project S0 with IDs listed in ID using forward intersection.
%   Points visible in too few images are given NaN coordinates.
%
%   S=FORWINTERSECT(S0,'all') computes all 3D coordinates.
%
%   [S,ID,RES]=... returns the IDs of each computed points in I and the
%   corresponding rms object space residual in RES.
%
%See also: PM_MULTIFORWINTERSECT.

% $Id$

