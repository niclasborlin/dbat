clear

genfigures=true;

% Extract name of current directory.
curDir=fileparts(mfilename('fullpath'));

% Defult to Olympus Camedia C4040Z camera calibration dataset if no data
% file is specified.
if ~exist('fName','var')
    fName=fullfile(curDir,'data','C4040Z-2272x1704.txt');
    fprintf('No data file specified, using ''%s''.\n',fName);
    disp(['Set variable ''fName'' to name of Photomodeler Export file if ' ...
          'you wish to use another file.']);
    disp(' ')
end
    
% Load data...
fprintf('Loading data file %s...',fName);
prob=loadpm(fName);
disp('done.')

% ...convert to useful struct...
s=prob2dbatstruct(prob);

% ...plot it.
h=plotnetwork(s,'title','Camera calibration data set (computed by Photomodeler)','camerasize',0.2);

if genfigures
    h=get(h,'parent');
    figDir=fullfile('..','doc','manual','ill');
    files={'ccam.eps'};
    for i=1:length(h)
        print(h(i),'-depsc2',fullfile(figDir,files{i}));
    end
end
