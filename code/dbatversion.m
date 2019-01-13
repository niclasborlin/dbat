function [v,date]=dbatversion(gitToo)
%DBATVERSION Version information for the Damped Bundle Adjustment Toolbox.
%
%   DBATVERSION returns a string with the Damped Bundle Adjustment
%   Toolbox version. The string has the format major.minor.patch.revision.
%
%   [V,D]=DBATVERSION also returns a UTC date string D.
%
%   DBATVERSION(TRUE) appends git branch/describe info. For
%   development purposes only.

if nargin<1, gitToo=false; end

% Should always be x.y.z or x.y.z.w.
v='0.8.5.1';

date='2019-01-13';

if gitToo
    g=gitversion;
    if ~isempty(g)
        v=[v,' (',g,')'];
    end
end

function s=gitversion
% Get git branch/describe info.

% Determine if we need to cd to the dbat working tree.
curDir=fullfile(pwd,filesep);
dbatDir=fullfile(fileparts(mfilename('fullpath')),filesep);
goodDir=strncmp(curDir,dbatDir,length(dbatDir));

if ~goodDir
    cd(dbatDir);
end

% Get branch name.
gitCmd='git rev-parse --abbrev-ref HEAD';
[err,branch]=system(gitCmd);
if err
    s='<git error>';
    disp(['Git error when executing '',gitCmd,''']);
    disp(branch)
    return;
end
branch=trimprompt(strtrim(branch));

% Get description.
gitCmd='git describe';
[err,desc]=system(gitCmd);
if err
    s='<git error>';
    disp(['Git error when executing '',gitCmd,''']);
    disp(desc)
    return;
end
desc=trimprompt(strtrim(desc));

s=['<',branch,'>, ',desc];

if ~goodDir
    cd(curDir);
end

function s=trimprompt(s)

% Sigh. Fix to remove any shell codes, unix prompts, etc.
% Remove everything before last whitespace or ctrl char.
bad=isstrprop(s,'wspace') | isstrprop(s,'cntrl');
lastBad=find(bad,1,'last');
if ~isempty(lastBad)
    s=s(lastBad+1:end);
end
