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

dampings={'none','gna','lm','lmp'};

dampings=dampings(2);

result=cell(size(dampings));
ok=nan(size(dampings));
iters=nan(size(dampings));
sigma0=nan(size(dampings));
E=cell(size(dampings));

for i=1:length(dampings)
    fprintf('Running the bundle with damping %s...\n',dampings{i});

    % Run the bundle.
    [result{i},ok(i),iters(i),sigma0(i),E{i}]=bundle(s0,dampings{i},'trace');

    if ok(i)
        fprintf('Bundle ok after %d iterations with sigma0=%.2f pixels\n', ...
                iters(i),sigma0(i));
    else
        fprintf(['Bundle failed after %d iterations. Last sigma0 estimate=%.2f ' ...
                 'pixels\n'],iters(i),sigma0(i));
    end
end

% Rotate to have +Z up.
T0=blkdiag(1,[0,-1;1,0],1);

for i=1:length(E)
    plotnetwork(result{i},E{i},'trans',T0,'align',1,'title',...
                ['Damping: ',dampings{i},'. Iteration %d of %d'], ...
                'axes',figure(i),'pause','on');
end
