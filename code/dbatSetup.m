% Get dir of executing (i.e. this!) file.
baseDir=fileparts(mfilename('fullpath'));

% Add selected subdirectories.
addpath(fullfile(baseDir,'plotting'),'-end')
addpath(fullfile(baseDir,'file'),'-end')
addpath(fullfile(baseDir,'misc'),'-end')
addpath(fullfile(baseDir,'demo'),'-end')

disp('You can now access DBAT from everywhere.')

clear baseDir
