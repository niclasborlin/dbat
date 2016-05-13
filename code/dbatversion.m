function [v,date]=dbatversion
%DBATVERSION Version information for the Damped Bundle Adjustment Toolbox.
%
%   DBATVERSION returns a string with the Damped Bundle Adjustment
%   Toolbox version. The string has the format major.minor.patch.revision.
%
%   [V,D]=DBATVERSION also returns a UTC date string D.

% Should always be x.y.z or x.y.z.w.
v='0.5.0';

date='2016-05-13';
