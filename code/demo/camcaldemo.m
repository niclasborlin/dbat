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

% Estimate px,py,c,K1-K3,P1-P2.
s0.cIO(1:8,:)=true;

% Set initial IO parameters.
s0.IO(1)=s0.IO(11)/2;  % px = center of sensor
s0.IO(2)=-s0.IO(12)/2; % py = center of sensor (sign is due to camera model)
s0.IO(3)=7.3;          % c = EXIF value.
s0.IO(4:8)=0;          % K1-K3, P1-P2 = 0.

cpId=1001:1004;
s0.isCtrl=ismember(s0.OPid,cpId);

% Fix the bundle datum by fixing all control points.
s0.cOP(:,ismember(s0.OPid,cpId))=false;

s1=resect(s0,'all',cpId,1,0,cpId);
s2=forwintersect(s1,'all',true);

s2.x0desc='Camera calibration from EXIF value';

s=s2;
h=plotnetwork(s,'title','Initial network',...
              'axes',tagfigure(sprintf('network%d',i)),...
              'camsize',0.1);

if printdemofigures
    h=get(h,'parent');
    figDir=fullfile('..','docsrc','manual','ill');
    files={'ccamx0.eps'};
    for i=1:length(h)
        print(h(i),'-depsc2',fullfile(figDir,files{i}));
    end
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

COP=bundle_result_file(result{1},E{1},resFile);

fprintf('\nBundle result file %s generated.\n',resFile);

h=plotparams(result{1},E{1});

if printdemofigures
    figDir=fullfile('..','docsrc','manual','ill');
    files={'ccamiotrace.eps','ccameotrace.eps','ccamoptrace.eps', ...
           'ccamgnatrace.eps'};
    for i=1:length(h)
        print(h(i),'-depsc2',fullfile(figDir,files{i}));
    end
end

h=plotcoverage(result{1},true);

if printdemofigures
    figDir=fullfile('..','docsrc','manual','ill');
    files={'ccamcoverage.eps'};
    for i=1:length(h)
        print(h(i),'-depsc2',fullfile(figDir,files{i}));
    end
end

h=plotimagestats(result{1},E{1});

if printdemofigures
    figDir=fullfile('..','docsrc','manual','ill');
    files={'ccamimstats.eps'};
    for i=1:length(h)
        print(h(i),'-depsc2',fullfile(figDir,files{i}));
    end
end

h=plotopstats(result{1},E{1},COP);

if printdemofigures
    figDir=fullfile('..','docsrc','manual','ill');
    files={'ccamopstats.eps'};
    for i=1:length(h)
        print(h(i),'-depsc2',fullfile(figDir,files{i}));
    end
end

if printdemofigures, doPause=0; else doPause='on'; end

for i=1:length(E)
    h=plotnetwork(result{i},E{i},'title',...
                  ['Damping: ',dampings{i},'. Iteration %d of %d'], ...
                  'axes',tagfigure(sprintf('network%d',i)),...
                  'pause',doPause,'camsize',0.1); 
end

if printdemofigures
    h=get(h,'parent');
    figDir=fullfile('..','docsrc','manual','ill');
    files={'ccamxfinal.eps'};
    for i=1:length(h)
        print(h(i),'-depsc2',fullfile(figDir,files{i}));
    end
end
