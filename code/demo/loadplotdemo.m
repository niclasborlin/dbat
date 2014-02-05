clear

genfigures=true;

% Extract name of current directory.
curDir=fileparts(mfilename('fullpath'));

% Defult to Roma dataset if no data file is specified.
if ~exist('fName','var')
    fName=fullfile(curDir,'data','roma.txt');
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
disp(['Plotting the camera network aligned with camera 1 and with +Z being ' ...
      '''up'''])
h=plotnetwork(s,'title','Roma data (computed by Photomodeler)','align',1,'trans','up');

if genfigures
    h=get(h,'parent');
    figDir=fullfile('..','doc','manual','ill');
    files={'roma.eps'};
    for i=1:length(h)
        print(h(i),'-depsc2',fullfile(figDir,files{i}));
    end
end
