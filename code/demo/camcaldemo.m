% Extract name of current directory.
curDir=fileparts(mfilename('fullpath'));

% Defult to Olympus Camedia C4040Z dataset if no data file is specified.
if ~exist('fName','var')
    fName=fullfile(curDir,'data','C4040Z-2272x1704.txt');
    fprintf('No data file specified, using ''%s''.\n',fName);
    disp(['Set variable ''fName'' to name of Photomodeler Export file if ' ...
          'you wish to use another file.']);
    disp(' ')
end

if ~exist('prob','var')
    fprintf('Loading data file %s...',fName);
    prob=loadpm(fName);
    if any(isnan(cat(2,prob.images.imSz)))
        error('Image sizes unknown!');
    end
    disp('done.')
else
    disp('Using pre-loaded data. Do ''clear prob'' to reload.');
end
s0=prob2dbatstruct(prob);

% Set CP 1-4 to nominal coordinates. Initial values for the EO and OP
% parameters are computed based on these points.
s0.OP(:)=nan;
s0.OP(:,ismember(s0.OPid,1001))=[0,1,0]';
s0.OP(:,ismember(s0.OPid,1002))=[1,1,0]';
s0.OP(:,ismember(s0.OPid,1003))=[0,0,0]';
s0.OP(:,ismember(s0.OPid,1004))=[1,0,0]';

% Fix the bundle datum by fixing CP 1, 2 + the Z coordinate of OP 3.
s0.cOP(:,ismember(s0.OPid,1001:1002))=false;
s0.cOP(3,ismember(s0.OPid,1003))=false;

% Estimate px,py,c,K1-K3,P1-P2.
s0.cIO(1:8,:)=true;

% Set initial IO parameters.
s0.IO(1)=s0.IO(11)/2;  % px = center of sensor
s0.IO(2)=-s0.IO(12)/2; % py = center of sensor (sign is due to camera model)
s0.IO(3)=s0.IO(3);     % c = half of true
s0.IO(4:8)=0;          % K1-K3, P1-P2 = 0.

cpId=1001:1004;
s1=resect(s0,'all',cpId,cpId);
s2=forwintersect(s1,'all',true);

s2.x0desc='Camera calibration';

s=s2;
plotnetwork(s,'title','Initial network',...
            'axes',tagfigure(sprintf('network%d',i)),'pause','on',...
            'camerasize',0.1);

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
    [result{i},ok(i),iters(i),sigma0(i),E{i}]=bundle(s,dampings{i},'trace');
    
    if ok(i)
        fprintf('Bundle ok after %d iterations with sigma0=%.2f pixels\n', ...
                iters(i),sigma0(i));
    else
        fprintf(['Bundle failed after %d iterations. Last sigma0 estimate=%.2f ' ...
                 'pixels\n'],iters(i),sigma0(i));
    end
end

resFile=strrep(fName,'.txt','_result_file.txt');

bundle_result_file(result{1},E{1},resFile);

fprintf('\nBundle result file %s generated.\n',resFile);

plotparams(result{1},E{1});

plotcoverage(result{1},true);

plotimagestats(result{1},E{1});

plotopstats(result{1},E{1});

for i=1:length(E)
    plotnetwork(result{i},E{i},'title',...
                ['Damping: ',dampings{i},'. Iteration %d of %d'], ...
                'axes',tagfigure(sprintf('network%d',i)),'pause','on','camerasize',0.1);
end
