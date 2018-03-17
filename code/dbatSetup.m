function dbatSetup(dev)
%DBATSETUP Setup paths, etc., for the Damped Bundle Adjustment Toolbox.

if nargin<1, dev=false; end

% Get dir of executing (i.e. this!) file.
baseDir=fileparts(mfilename('fullpath'));

% Add selected subdirectories.
addpath(fullfile(baseDir,'plotting'),'-end')
addpath(fullfile(baseDir,'file'),'-end')
addpath(fullfile(baseDir,'misc'),'-end')
addpath(fullfile(baseDir,'bundle'),'-end')
addpath(fullfile(baseDir,'bundle','cammodel'),'-end')
addpath(fullfile(baseDir,'bundle','cameramodel'),'-end')
addpath(fullfile(baseDir,'bundle','lsa'),'-end')
addpath(fullfile(baseDir,'photogrammetry'),'-end')
addpath(fullfile(baseDir,'demo'),'-end')
addpath(fullfile(baseDir,'xchg','comp_struct'),'-end')
addpath(fullfile(baseDir,'xchg','ply'),'-end')
addpath(fullfile(baseDir,'xchg','struct2xml'),'-end')
addpath(fullfile(baseDir,'xchg','xml2struct'),'-end')
addpath(baseDir,'-end')

fprintf('You can now access DBAT v%s from everywhere.\n',dbatversion(dev))
