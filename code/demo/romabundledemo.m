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

% Fix the datum by fixing camera 1...
s0.cEO(:,1)=false;
% ...and the largest other absolute camera coordinate.
camDiff=abs(s0.EO(1:3,:)-repmat(s0.EO(1:3,1),1,size(s0.EO,2)));
[i,j]=find(camDiff==max(camDiff(:)));
s0.cEO(i,j)=false;

disp('Running the bundle');
% Run the bundle.
[s1,ok,iters,s0,E]=bundle(s0,'none','trace');

if ok
    fprintf('Bundle ok after %d iterations with sigma0=%.2f pixels\n',iters,s0);
else
    fprintf('Bundle failed after %d iterations. Last sigma0 estimate=%.2f pixels\n',iters, s0);
end

% Rotate to have +Z up.
T0=blkdiag(1,[0,-1;1,0],1);

plotnetwork(s1,E,'trans',T0,'align',1,'title','Iteration %d of %d', ...
            'pause','on');