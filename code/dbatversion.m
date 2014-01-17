function [v,date]=dbatversion
%DBATVERSION Version information for the Damped Bundle Adjustment Toolbox.
%
%   DBATVERSION returns a string with the Damped Bundle Adjustment
%   Toolbox version. The string has the format major.minor.patch.revision.
%
%   [V,D]=DBATVERSION also returns a UTC date string D.

% $Id$

% These strings are updated by subversion on commit.
idStr="$Id$";
dateStr="$Date$";
revisionStr="$Revision$";

% Parse date string
d=