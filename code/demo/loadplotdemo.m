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
    
fprintf('Loading data file %s...',fName);
prob=loadpm(fName);
disp('done.')

s=prob2dbatstruct(prob);

% Rotate to have +Z up.
T0=blkdiag(1,[0,-1;1,0],1);

plotnetwork(s,'title','Roma','align',1,'trans',T0);
