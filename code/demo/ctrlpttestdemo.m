% Extract name of current directory.
curDir=fileparts(mfilename('fullpath'));

% Defult to Roma dataset if no data file is specified.
if ~exist('fName','var')
    fName=fullfile(curDir,'data','ctrlpttest.txt');
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

% Don't estimate IO data (treat it as exact).
s0.IO=s0.IOobs; % No really necessary...
s0.estIO=false(size(s0.IO));
s0.useIOobs=false(size(s0.IO));

% Use supplied EO data as initial values. Treat EO data as free.
s0.EO=s0.EOobs; % No really necessary...
s0.estEO(1:6,:)=true;
s0.useEOobs=false(size(s0.EO));

% Use supplied OP data as initial values. Treat control points as
% exact.
s0.OP=s0.OPobs; % No really necessary...
s0.estOP=repmat(~s0.isCtrl(:)',3,1);
s0.useOPobs=repmat(s0.isCtrl(:)',3,1);

% Adjust EO and OP to not be exactly what was given by PM.
s0.EO(1:3,:)=round(s0.EO(1:3,:)*10)/10;
s0.EO(4:6,:)=round(s0.EO(4:6,:)*10)/10;
s0.OP=round(s0.OP*10)/10;

% Use sigma0=1 as first approximation.
s0.markStd(:)=1;

% The datum by fixed already via the control points...
if 0
    s0.estEO(:,1)=false;
    % ...and the largest other absolute camera coordinate.
    camDiff=abs(s0.EO(1:3,:)-repmat(s0.EO(1:3,1),1,size(s0.EO,2)));
    [i,j]=find(camDiff==max(camDiff(:)));
    s0.estEO(i,j)=false;
end

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

resFile=strrep(fName,'.txt','_result_file.txt');

COP=bundle_result_file(result{1},E{1},resFile);

fprintf('\nBundle result file %s generated.\n',resFile);

plotparams(result{1},E{1},'noop');

plotcoverage(result{1},true);

plotimagestats(result{1},E{1});

plotopstats(result{1},E{1},COP);

for i=1:length(E)
    fig=tagfigure(sprintf('network%d',i));
    fprintf('Displaying bundle iteration playback for method %s in figure %d.\n',E{i}.damping.name,double(fig));
    plotnetwork(result{i},E{i},'trans','up','align',1,'title',...
                ['Damping: ',dampings{i},'. Iteration %d of %d'], ...
                'axes',fig,'pause','on');
end
