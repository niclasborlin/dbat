function [v,date]=dbatversion
%DBATVERSION Version information for the Damped Bundle Adjustment Toolbox.
%
%   DBATVERSION returns a string with the Damped Bundle Adjustment
%   Toolbox version. The string has the format major.minor.patch.revision.
%
%   [V,D]=DBATVERSION also returns a UTC date string D.

% $Id$

% Should always be x.y.z.
toolboxVersion='0.4.1';

% This strings is updated by subversion on commit.
idStr='$Id$';

% Parse string to get revision and UTC date.
z=regexp(idStr,['\$Id:\s*\S+\s+(\d+)\s+([\d-]+\s[\d:]+)Z'],'tokens','once');

if 0
rev=z{1};
date=z{2};
else
    rev='git-1-gc4a63c8';
    date='2016-03-27';
end
v=[toolboxVersion,'.',rev];
