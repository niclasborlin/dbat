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

if ~exist('prob','var')
    fprintf('Loading data file %s...',fName);
    prob=loadpm(fName);
    disp('done.')
else
    disp('Using pre-loaded data. Do ''clear prob'' to reload.');
end
s0=prob2dbatstruct(prob);

% Solve the datum problem by fixing camera 1...
s0.cEO(:,1)=false;
% ...and the largest absolute camera coordinate.
camDiff=abs(s0.EO(1:3,:)-repmat(s0.EO(1:3,1),1,size(s0.EO,2)));
[i,j]=find(camDiff==max(camDiff(:)));
s0.cEO(i,j)=false;

% Run the bundle.
s1=bundle(s0,'none');
